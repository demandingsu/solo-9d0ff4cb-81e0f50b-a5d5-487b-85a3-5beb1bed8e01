"""两运行读段内容重合比对（后台任务）。

用途：怀疑两个运行的 FASTQ 被重复投递时，按"读段内容"定位重合读段，
交给人工复核后再决定是否进入下游比对。

分组规则：忽略读段名，仅当 R1、R2 的序列与质量串四条字节串分别
逐字节相同时才视为同组；同组在两边各出现 c_a、c_b 次时，按
min(c_a, c_b) 计入重合对数。

设计约束（对应交付要求）：
* 快照一致：对主库只开一个只读事务（WAL 快照），两个运行的状态、
  已提交对数与全部样例都来自同一次已提交快照；采集继续提交新配对
  不会影响本次比对结果；
* 不阻塞采集：比对在独立后台线程执行，主库走只读快照事务
  （WAL 模式下读者不阻塞写者），分组数据写入独立的临时 SQLite 库；
* 内存有界：pairs 表逐对流式读取，分组键落临时库磁盘表，
  任何时刻只持有一个批次，大运行不会整份载入内存；
* 失败安全：主库连接置 query_only，从机制上无法写主库；临时库
  独立存放、结束即清理；任务失败只影响任务本身的状态记录，
  原始读段与采集进度不受影响。

注意：分组键是四条字节串的 128-bit SHA-256 截断摘要（带长度前缀，
无拼接歧义）。128-bit 下生日碰撞需要约 2^64 条记录，工程上可视为
与逐字节比较等价；如需严格逐字节可改为在临时库存放完整键。
"""
from __future__ import annotations

import collections
import hashlib
import os
import shutil
import sqlite3
import tempfile
import threading
import time
import uuid
from typing import Any

from . import db

# 结果中样例分组数与每组每侧序号数的上限默认值
SAMPLE_GROUPS = 10
SAMPLE_INDICES = 32
# 批量写入临时库的行数（内存占用上界 ≈ 批大小 × 行宽）
INSERT_BATCH = 10_000
# 序列预览截断长度（仅供人工复核快速辨认分组）
PREVIEW_LEN = 80


class CompareError(RuntimeError):
    """比对过程中的可预期失败（数据损坏、运行不存在等）。"""


def _seq_qual(raw: bytes) -> tuple[bytes, bytes]:
    """从一条完整四行 FASTQ 记录字节中取出序列行与质量行。

    记录落库前已通过 parser.validate 的结构校验（严格四行、序质等长），
    这里只做防御性检查：损坏即抛错，任务失败但不影响主库数据。
    """
    parts = raw.split(b"\n")
    if len(parts) < 5 or not parts[0].startswith(b"@"):
        raise CompareError(f"pairs 表中的记录字节损坏: {raw[:40]!r}")
    return parts[1], parts[3]


def group_key(r1_record: bytes, r2_record: bytes) -> bytes:
    """一对读段的分组键：忽略读段名，只看 R1/R2 的序列与质量串。

    每个字段带 4 字节长度前缀后再哈希，任意字节内容都不会产生拼接歧义。
    """
    s1, q1 = _seq_qual(r1_record)
    s2, q2 = _seq_qual(r2_record)
    h = hashlib.sha256()
    for field in (s1, q1, s2, q2):
        h.update(len(field).to_bytes(4, "little"))
        h.update(field)
    return h.digest()[:16]


def _load_side(
    src: sqlite3.Connection, tmp: sqlite3.Connection, run_id: str, side: int
) -> int:
    """把一侧运行的全部配对流式读出，分组键写入临时库。返回读取条数。"""
    cur = src.execute(
        "SELECT seq, r1_record, r2_record FROM pairs WHERE run_id=? ORDER BY seq",
        (run_id,),
    )
    batch: list[tuple[int, bytes, int]] = []
    n = 0
    tmp.execute("BEGIN")
    for row in cur:
        batch.append((side, group_key(row["r1_record"], row["r2_record"]), row["seq"]))
        n += 1
        if len(batch) >= INSERT_BATCH:
            tmp.executemany("INSERT INTO g(side, h, seq) VALUES (?,?,?)", batch)
            tmp.commit()
            tmp.execute("BEGIN")
            batch.clear()
    if batch:
        tmp.executemany("INSERT INTO g(side, h, seq) VALUES (?,?,?)", batch)
    tmp.commit()
    return n


def _indices(
    tmp: sqlite3.Connection, side: int, h: bytes, cap: int
) -> tuple[list[int], bool]:
    """某分组在一侧的配对序号（升序、封顶 cap 条），并报告是否被截断。"""
    rows = tmp.execute(
        "SELECT seq FROM g WHERE side=? AND h=? ORDER BY seq LIMIT ?",
        (side, h, cap + 1),
    ).fetchall()
    seqs = [r[0] for r in rows]
    return seqs[:cap], len(seqs) > cap


def _preview(src: sqlite3.Connection, run_id: str, seq: int) -> dict[str, Any] | None:
    """取某分组一条代表记录的序列预览（截断），供人工复核快速辨认。"""
    row = src.execute(
        "SELECT r1_record, r2_record FROM pairs WHERE run_id=? AND seq=?",
        (run_id, seq),
    ).fetchone()
    if row is None:
        return None
    s1, _ = _seq_qual(row["r1_record"])
    s2, _ = _seq_qual(row["r2_record"])
    return {
        "r1_seq": s1[:PREVIEW_LEN].decode("ascii", "replace"),
        "r1_seq_len": len(s1),
        "r2_seq": s2[:PREVIEW_LEN].decode("ascii", "replace"),
        "r2_seq_len": len(s2),
    }


def compare_runs(
    db_path: str,
    run_a: str,
    run_b: str,
    *,
    tmp_root: str | None = None,
    sample_groups: int = SAMPLE_GROUPS,
    sample_indices: int = SAMPLE_INDICES,
) -> dict[str, Any]:
    """比对两个运行的读段内容重合情况。全部读操作在同一个主库快照事务内完成。

    返回结果字典（见 CompareJob.public 中的 result 部分）。
    任何异常都不会改动主库数据：主库连接为 query_only，分组数据只写临时库。
    """
    if tmp_root is not None:
        os.makedirs(tmp_root, exist_ok=True)
    tmp_dir = tempfile.mkdtemp(prefix="compare-", dir=tmp_root)
    src: sqlite3.Connection | None = None
    tmp: sqlite3.Connection | None = None
    try:
        src = db.connect(db_path)
        # 机制性防写：即使后续代码有缺陷也无法改动主库
        src.execute("PRAGMA query_only=ON")
        # 固定一次已提交快照：之后的状态、计数、配对字节全部来自该快照；
        # WAL 模式下此读事务不阻塞采集写事务
        src.execute("BEGIN")

        runs: dict[str, dict[str, Any]] = {}
        for rid in (run_a, run_b):
            row = src.execute("SELECT * FROM runs WHERE run_id=?", (rid,)).fetchone()
            if row is None:
                raise CompareError(f"运行 {rid} 不存在")
            runs[rid] = {
                "run_id": row["run_id"],
                "status": row["status"],
                "pairs": row["pairs"],
                "r1_offset": row["r1_offset"],
                "r2_offset": row["r2_offset"],
            }

        # 独立临时库：分组数据落盘，内存只持有一个批次
        tmp = sqlite3.connect(os.path.join(tmp_dir, "groups.db"), isolation_level=None)
        tmp.execute("PRAGMA journal_mode=OFF")
        tmp.execute("PRAGMA synchronous=OFF")
        tmp.execute(
            "CREATE TABLE g (side INTEGER NOT NULL, h BLOB NOT NULL, seq INTEGER NOT NULL)"
        )
        loaded = {
            run_a: _load_side(src, tmp, run_a, 0),
            run_b: _load_side(src, tmp, run_b, 1),
        }
        # 与快照中的已提交对数互相印证（同一事务内必然一致）
        for rid, info in runs.items():
            if info["pairs"] != loaded[rid]:
                raise CompareError(
                    f"运行 {rid} 快照内计数不一致: runs={info['pairs']} pairs={loaded[rid]}"
                )
        tmp.execute("CREATE INDEX idx_g_side_h ON g(side, h)")
        tmp.execute("CREATE TABLE ca (h BLOB PRIMARY KEY, c INTEGER NOT NULL)")
        tmp.execute("INSERT INTO ca SELECT h, COUNT(*) FROM g WHERE side=0 GROUP BY h")
        tmp.execute("CREATE TABLE cb (h BLOB PRIMARY KEY, c INTEGER NOT NULL)")
        tmp.execute("INSERT INTO cb SELECT h, COUNT(*) FROM g WHERE side=1 GROUP BY h")

        overlap_pairs, overlap_groups = tmp.execute(
            "SELECT COALESCE(SUM(MIN(ca.c, cb.c)), 0), COUNT(*)"
            " FROM ca JOIN cb ON ca.h = cb.h"
        ).fetchone()

        samples: list[dict[str, Any]] = []
        for (h,) in tmp.execute(
            "SELECT ca.h FROM ca JOIN cb ON ca.h = cb.h"
            " ORDER BY MIN(ca.c, cb.c) DESC, ca.h LIMIT ?",
            (sample_groups,),
        ).fetchall():
            count_a = tmp.execute("SELECT c FROM ca WHERE h=?", (h,)).fetchone()[0]
            count_b = tmp.execute("SELECT c FROM cb WHERE h=?", (h,)).fetchone()[0]
            idx_a, trunc_a = _indices(tmp, 0, h, sample_indices)
            idx_b, trunc_b = _indices(tmp, 1, h, sample_indices)
            samples.append(
                {
                    "count_a": count_a,
                    "count_b": count_b,
                    "overlap_pairs": min(count_a, count_b),
                    "pair_indices_a": idx_a,
                    "pair_indices_b": idx_b,
                    "truncated_a": trunc_a,
                    "truncated_b": trunc_b,
                    "preview": _preview(src, run_a, idx_a[0]),
                }
            )

        src.rollback()  # 释放快照（只读事务，无内容可提交）
        return {
            "run_a": runs[run_a],
            "run_b": runs[run_b],
            "overlap_pairs": overlap_pairs,
            "overlap_groups": overlap_groups,
            "sample_groups": len(samples),
            "samples": samples,
        }
    finally:
        if src is not None:
            src.close()
        if tmp is not None:
            tmp.close()
        shutil.rmtree(tmp_dir, ignore_errors=True)


class CompareJob:
    """一次比对任务的状态记录（仅内存，进程重启即消失）。"""

    def __init__(self, run_a: str, run_b: str):
        self.job_id = uuid.uuid4().hex
        self.run_a = run_a
        self.run_b = run_b
        self.status = "pending"  # pending|running|done|failed
        self.created_at = time.time()
        self.started_at: float | None = None
        self.finished_at: float | None = None
        self.result: dict[str, Any] | None = None
        self.error: str | None = None

    def public(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "job_id": self.job_id,
            "run_a": self.run_a,
            "run_b": self.run_b,
            "status": self.status,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }
        if self.status == "done" and self.result is not None:
            out["result"] = self.result
            out["duration_sec"] = round(self.finished_at - self.started_at, 3)
        if self.status == "failed":
            out["error"] = self.error
        return out


class CompareManager:
    """比对任务队列：有界后台线程执行，避免慢比对拖垮采集进程。"""

    def __init__(
        self,
        db_path: str,
        *,
        max_workers: int = 2,
        keep_jobs: int = 200,
        tmp_root: str | None = None,
        sample_groups: int = SAMPLE_GROUPS,
        sample_indices: int = SAMPLE_INDICES,
    ):
        self.db_path = db_path
        self.max_workers = max(1, max_workers)
        self.keep_jobs = max(10, keep_jobs)
        # 临时库默认放在主库所在目录（通常是持久卷），避免小容量 tmpfs
        self.tmp_root = tmp_root or os.path.join(
            os.path.dirname(os.path.abspath(db_path)), "compare-tmp"
        )
        self.sample_groups = sample_groups
        self.sample_indices = sample_indices
        self._lock = threading.Lock()
        self._jobs: dict[str, CompareJob] = {}
        self._order: collections.deque[str] = collections.deque()
        self._pending: collections.deque[str] = collections.deque()
        self._active = 0

    # ------------------------------------------------------------------ 查询
    def get(self, job_id: str) -> CompareJob | None:
        with self._lock:
            return self._jobs.get(job_id)

    def list(self) -> list[CompareJob]:
        with self._lock:
            return [self._jobs[jid] for jid in self._order if jid in self._jobs]

    # ------------------------------------------------------------------ 提交
    def submit(self, run_a: str, run_b: str) -> CompareJob:
        job = CompareJob(run_a, run_b)
        with self._lock:
            self._jobs[job.job_id] = job
            self._order.append(job.job_id)
            self._pending.append(job.job_id)
            self._evict_locked()
            self._pump_locked()
        return job

    def _evict_locked(self) -> None:
        """任务记录只保留最近 keep_jobs 条，优先淘汰已结束的。"""
        while len(self._jobs) > self.keep_jobs:
            victim = None
            for jid in self._order:
                job = self._jobs.get(jid)
                if job is not None and job.status in ("done", "failed"):
                    victim = jid
                    break
            if victim is None:
                break  # 全是 pending/running，不淘汰
            del self._jobs[victim]
        if len(self._order) > len(self._jobs):
            self._order = collections.deque(
                jid for jid in self._order if jid in self._jobs
            )

    def _pump_locked(self) -> None:
        while self._active < self.max_workers and self._pending:
            jid = self._pending.popleft()
            job = self._jobs.get(jid)
            if job is None:
                continue
            self._active += 1
            job.status = "running"
            job.started_at = time.time()
            t = threading.Thread(
                target=self._run_guarded,
                args=(job,),
                name=f"compare-{jid[:8]}",
                daemon=True,
            )
            t.start()

    # ------------------------------------------------------------------ 执行
    def _run_guarded(self, job: CompareJob) -> None:
        try:
            job.result = compare_runs(
                self.db_path,
                job.run_a,
                job.run_b,
                tmp_root=self.tmp_root,
                sample_groups=self.sample_groups,
                sample_indices=self.sample_indices,
            )
            job.status = "done"
        except Exception as exc:  # 失败只记录到任务，不触碰主库数据
            job.status = "failed"
            job.error = f"{type(exc).__name__}: {exc}"
        finally:
            job.finished_at = time.time()
            with self._lock:
                self._active -= 1
                self._pump_locked()
