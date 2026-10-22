"""Additional intermediate-frequency validity boundary and far-limit checks."""

from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time
from pair_asymptotics import work
from common import resumable_run, finish, atomic_json


def main():
    start = time.perf_counter()
    out = resumable_run(
        "pair_asymptotics",
        dict(extension="intermediate f plus formal far-limit", workers=6),
        ["Zhang2011 Eqs3-11; DoinikovBouakaz2015 Eqs8-13,39-43,67-68"],
    )
    parent = json.loads(
        (
            out / "PARENT_STATE.json"
            if (out / "PARENT_STATE.json").exists()
            else Path("results/stage2_5/pair_asymptotics/STATE.json")
        ).read_text()
    )
    atomic_json(out / "PARENT_STATE.json", parent)
    jobs = [
        (k, q, f, p, 1482.0, str(out))
        for k in ["zhang", "doinikov"]
        for q in [10, 20, 40, 80, 160]
        for f in [1e6, 1.5e6]
        for p in [100, 500]
    ]
    jobs += [
        (k, q, f, 100, 1e9, str(out))
        for k in ["zhang", "doinikov"]
        for q in [320, 640]
        for f in [2.15e6, 2.95e6]
    ]
    with ProcessPoolExecutor(max_workers=6) as pool:
        rows = list(pool.map(work, jobs))
    summary = dict(
        status="COMPARISON_COMPLETE",
        physics_status="partially validated",
        rows=parent["rows"] + rows,
        parent_run=parent["run"],
        scaling=parent["scaling"],
        extension_runtime_s=time.perf_counter() - start,
        original_runtime_s=parent.get("runtime_s"),
        additional_cases=len(rows),
    )
    atomic_json(out / "PARENT_STATE.json", parent)
    finish(out, summary)


if __name__ == "__main__":
    main()
