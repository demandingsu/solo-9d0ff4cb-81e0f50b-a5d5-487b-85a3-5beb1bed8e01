"""读段内容重合比对：分组键单元测试 + 失败安全性 + 真实子进程集成测试。"""
import time

import pytest

from collector import compare, db
from tests.test_integration import _write_run, make_scenario, poll_collect


# ------------------------------------------------------------------ 单元测试
def test_group_key_ignores_read_names():
    r1a = b"@alpha/1\nACGT\n+\nIIII\n"
    r2a = b"@alpha/2\nTGCA\n+\nJJJJ\n"
    # 读段名不同、序列与质量串完全相同 -> 同组
    r1b = b"@beta/1\nACGT\n+\nIIII\n"
    r2b = b"@beta/2\nTGCA\n+\nJJJJ\n"
    assert compare.group_key(r1a, r2a) == compare.group_key(r1b, r2b)
    # R2 质量串不同 -> 不同组
    assert compare.group_key(r1a, r2a) != compare.group_key(r1b, b"@beta/2\nTGCA\n+\nJJJI\n")
    # R1 序列不同 -> 不同组
    assert compare.group_key(r1a, r2a) != compare.group_key(b"@beta/1\nACGA\n+\nIIII\n", r2b)
    # R1 质量串不同 -> 不同组
    assert compare.group_key(r1a, r2a) != compare.group_key(b"@beta/1\nACGT\n+\nIIIJ\n", r2b)


def test_group_key_rejects_corrupted_record():
    with pytest.raises(compare.CompareError):
        compare.group_key(b"garbage", b"@x/2\nACGT\n+\nIIII\n")


def test_compare_failure_preserves_data_and_cleans_tmp(tmp_path):
    """比对失败（pairs 中记录字节损坏）不得改动原始读段与进度，临时目录须清理。"""
    db_path = str(tmp_path / "c.db")
    conn = db.connect(db_path)
    db.init_db(conn)
    now = time.time()
    db.create_run(conn, "ra", "http://x/ra", now)
    db.create_run(conn, "rb", "http://x/rb", now)
    good1 = b"@a0/1\nACGT\n+\nIIII\n"
    good2 = b"@a0/2\nTGCA\n+\nJJJJ\n"
    db.commit_pair(conn, "ra", 0, "a0", good1, good2, len(good1), len(good2), now)
    db.commit_pair(conn, "rb", 0, "b0", b"corrupted", good2, 9, len(good2), now)
    conn.close()

    with pytest.raises(compare.CompareError):
        compare.compare_runs(db_path, "ra", "rb", tmp_root=str(tmp_path))

    # 原始读段与进度原样保留
    conn2 = db.connect(db_path)
    ra, rb = db.get_run(conn2, "ra"), db.get_run(conn2, "rb")
    assert ra["pairs"] == 1 and ra["r1_offset"] == len(good1)
    assert rb["pairs"] == 1 and rb["r1_offset"] == 9
    assert list(db.iter_records(conn2, "rb", 1)) == [b"corrupted"]
    conn2.close()
    # 临时目录已清理
    assert list(tmp_path.glob("compare-*")) == []


def test_compare_counts_and_samples_directly(tmp_path):
    """不经过 HTTP，直接验证 min 计数与样例序号。"""
    db_path = str(tmp_path / "c.db")
    conn = db.connect(db_path)
    db.init_db(conn)
    now = time.time()

    def pair(name, s1, q1, s2, q2):
        r1 = f"@{name}/1\n{s1}\n+\n{q1}\n".encode()
        r2 = f"@{name}/2\n{s2}\n+\n{q2}\n".encode()
        return r1, r2

    db.create_run(conn, "ra", "http://x/ra", now)
    db.create_run(conn, "rb", "http://x/rb", now)
    # ra: X, Y, X, Z（X 出现 2 次）
    xa = pair("a0", "ACGT", "IIII", "TGCA", "JJJJ")
    ya = pair("a1", "GGGG", "HHHH", "CCCC", "FFFF")
    za = pair("a3", "AAAA", "####", "TTTT", "$$$$")
    for i, (r1, r2) in enumerate([xa, ya, xa, za]):
        db.commit_pair(conn, "ra", i, f"a{i}", r1, r2, 0, 0, now)
    # rb: X, X, X, W（X 出现 3 次，读段名完全不同）
    xb = pair("b0", "ACGT", "IIII", "TGCA", "JJJJ")
    wb = pair("b3", "GCGC", "%%%%", "CGCG", "&&&&")
    for i, (r1, r2) in enumerate([xb, xb, xb, wb]):
        db.commit_pair(conn, "rb", i, f"b{i}", r1, r2, 0, 0, now)
    conn.close()

    res = compare.compare_runs(db_path, "ra", "rb", tmp_root=str(tmp_path))
    assert res["run_a"]["pairs"] == 4 and res["run_b"]["pairs"] == 4
    # X: min(2, 3) = 2；Y/Z/W 不重合
    assert res["overlap_pairs"] == 2
    assert res["overlap_groups"] == 1
    (s,) = res["samples"]
    assert s["count_a"] == 2 and s["count_b"] == 3
    assert s["overlap_pairs"] == 2
    assert s["pair_indices_a"] == [0, 2]
    assert s["pair_indices_b"] == [0, 1, 2]
    assert s["truncated_a"] is False and s["truncated_b"] is False
    assert s["preview"]["r1_seq"] == "ACGT" and s["preview"]["r2_seq"] == "TGCA"


# ------------------------------------------------------------------ 集成测试
def wait_job(client, collector, job_id, timeout=60):
    deadline = time.time() + timeout
    while time.time() < deadline:
        r = client.get(f"{collector}/compare/jobs/{job_id}")
        assert r.status_code == 200
        j = r.json()
        if j["status"] in ("done", "failed"):
            assert j["status"] == "done", j
            return j
        time.sleep(0.2)
    raise AssertionError(f"比对任务未在 {timeout}s 内完成")


def _fq(name, s1, q1, s2, q2):
    return (
        f"@{name}/1\n{s1}\n+\n{q1}\n".encode(),
        f"@{name}/2\n{s2}\n+\n{q2}\n".encode(),
    )


def _collect_done(client, collector, relay, meta, run_id, pairs):
    r1 = b"".join(p[0] for p in pairs)
    r2 = b"".join(p[1] for p in pairs)
    _write_run(meta, run_id, r1, r2)
    snap, _ = poll_collect(client, collector, run_id, f"{relay}/runs/{run_id}")
    assert snap["pairs"] == len(pairs)
    return snap


def test_compare_api_validation(servers, client):
    _relay, collector, _meta = servers
    r = client.post(f"{collector}/compare/jobs",
                    json={"run_a": "nope1", "run_b": "nope2"})
    assert r.status_code == 404
    r = client.post(f"{collector}/compare/jobs", json={"run_a": "x", "run_b": "x"})
    assert r.status_code == 422
    assert client.get(f"{collector}/compare/jobs/{'0' * 32}").status_code == 404


def test_compare_partial_overlap_min_count_and_samples(servers, client):
    relay, collector, meta = servers
    # 运行 A：X, Y, X, Z；运行 B：X, X, X, W；读段名全部不同
    pairs_a = [
        _fq("a0", "ACGT", "IIII", "TGCA", "JJJJ"),
        _fq("a1", "GGGG", "HHHH", "CCCC", "FFFF"),
        _fq("a2", "ACGT", "IIII", "TGCA", "JJJJ"),
        _fq("a3", "AAAA", "####", "TTTT", "$$$$"),
    ]
    pairs_b = [
        _fq("b0", "ACGT", "IIII", "TGCA", "JJJJ"),
        _fq("b1", "ACGT", "IIII", "TGCA", "JJJJ"),
        _fq("b2", "ACGT", "IIII", "TGCA", "JJJJ"),
        _fq("b3", "GCGC", "%%%%", "CGCG", "&&&&"),
    ]
    _collect_done(client, collector, relay, meta, "cmpA1", pairs_a)
    _collect_done(client, collector, relay, meta, "cmpB1", pairs_b)

    r = client.post(f"{collector}/compare/jobs",
                    json={"run_a": "cmpA1", "run_b": "cmpB1"})
    assert r.status_code == 202
    job = wait_job(client, collector, r.json()["job_id"])
    res = job["result"]
    assert res["run_a"]["run_id"] == "cmpA1"
    assert res["run_a"]["status"] == "completed"
    assert res["run_a"]["pairs"] == 4
    assert res["run_a"]["r1_offset"] > 0 and res["run_a"]["r2_offset"] > 0
    assert res["run_b"]["pairs"] == 4
    # X 组按两边次数较小值计入：min(2, 3) = 2
    assert res["overlap_pairs"] == 2
    assert res["overlap_groups"] == 1
    assert res["sample_groups"] == 1
    (s,) = res["samples"]
    assert s["count_a"] == 2 and s["count_b"] == 3
    assert s["pair_indices_a"] == [0, 2]
    assert s["pair_indices_b"] == [0, 1, 2]
    assert s["preview"]["r1_seq"] == "ACGT" and s["preview"]["r2_seq"] == "TGCA"

    # 任务列表可查
    jobs = client.get(f"{collector}/compare/jobs").json()["jobs"]
    assert any(j["job_id"] == job["job_id"] for j in jobs)


def test_compare_identical_content_different_names_full_overlap(servers, client):
    relay, collector, meta = servers
    # 两个运行每对内容完全相同、读段名完全不同 -> 全部重合
    pairs_a = [_fq(f"ra{i}", "ACGT", "IIII", "TGCA", "JJJJ") for i in range(3)]
    pairs_b = [_fq(f"rb{i}", "ACGT", "IIII", "TGCA", "JJJJ") for i in range(3)]
    _collect_done(client, collector, relay, meta, "cmpA2", pairs_a)
    _collect_done(client, collector, relay, meta, "cmpB2", pairs_b)

    r = client.post(f"{collector}/compare/jobs",
                    json={"run_a": "cmpA2", "run_b": "cmpB2"})
    job = wait_job(client, collector, r.json()["job_id"])
    res = job["result"]
    assert res["overlap_pairs"] == 3
    assert res["overlap_groups"] == 1
    (s,) = res["samples"]
    assert s["pair_indices_a"] == [0, 1, 2]
    assert s["pair_indices_b"] == [0, 1, 2]


def test_compare_no_overlap(servers, client):
    relay, collector, meta = servers
    pairs_a = [_fq("na0", "AAAA", "IIII", "TTTT", "JJJJ")]
    pairs_b = [_fq("nb0", "CCCC", "IIII", "GGGG", "JJJJ")]
    _collect_done(client, collector, relay, meta, "cmpA3", pairs_a)
    _collect_done(client, collector, relay, meta, "cmpB3", pairs_b)
    r = client.post(f"{collector}/compare/jobs",
                    json={"run_a": "cmpA3", "run_b": "cmpB3"})
    res = wait_job(client, collector, r.json()["job_id"])["result"]
    assert res["overlap_pairs"] == 0
    assert res["overlap_groups"] == 0
    assert res["samples"] == []


def test_compare_while_collecting_uses_one_snapshot_and_does_not_block(servers, client):
    """运行仍在采集时：计数/状态/样例来自同一次已提交快照，且采集不被阻塞。"""
    relay, collector, meta = servers
    # 参照运行：快速完成
    pairs_a = [_fq(f"sa{i}", "ACGT", "IIII", "TGCA", "JJJJ") for i in range(2)]
    _collect_done(client, collector, relay, meta, "cmpA4", pairs_a)

    # 目标运行：缓慢追加，采集期间保持 running
    run_b = "cmpB4"
    client.post(f"{relay}/runs",
                json=make_scenario(run_b, reads=6, append_delay=1.0,
                                   r2_extra_delay=0.2, chunk=1))
    base_url = f"{relay}/runs/{run_b}"
    # 采到至少 2 对后停手（采集是拉取式的，不再调用则进度冻结）
    deadline = time.time() + 30
    snap = None
    while time.time() < deadline:
        snap = client.post(f"{collector}/collect/runs",
                           json={"run_id": run_b, "base_url": base_url}).json()
        if snap["pairs"] >= 2:
            break
        time.sleep(0.2)
    assert snap and snap["status"] == "running" and snap["pairs"] >= 2
    frozen_pairs = snap["pairs"]

    r = client.post(f"{collector}/compare/jobs",
                    json={"run_a": "cmpA4", "run_b": run_b})
    assert r.status_code == 202
    job = wait_job(client, collector, r.json()["job_id"])
    res = job["result"]
    # 状态与计数来自同一次已提交快照：仍是 running，且对数与冻结时一致
    assert res["run_b"]["status"] == "running"
    assert res["run_b"]["pairs"] == frozen_pairs
    assert res["run_a"]["status"] == "completed" and res["run_a"]["pairs"] == 2
    # 计数与样例自洽：重合不超过两边较小值，样例序号都在已提交对数内
    assert 0 <= res["overlap_pairs"] <= min(res["run_a"]["pairs"], res["run_b"]["pairs"])
    for s in res["samples"]:
        assert all(i < res["run_a"]["pairs"] for i in s["pair_indices_a"])
        assert all(i < res["run_b"]["pairs"] for i in s["pair_indices_b"])

    # 慢比对结束后采集照常推进到完成，下载也不受影响
    snap, _ = poll_collect(client, collector, run_b, base_url, timeout=60)
    assert snap["pairs"] == 6
    got = client.get(f"{collector}/collect/runs/{run_b}/R1")
    assert got.status_code == 200
    src = client.get(f"{relay}/runs/{run_b}/R1.fastq")
    assert got.content == src.content
