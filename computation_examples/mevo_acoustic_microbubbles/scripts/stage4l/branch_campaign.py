"""Finite primary and dt/2 jobs; an explicit completed S0 gate is required."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from multiprocessing import get_context
from pathlib import Path

from scripts.stage4l.branch_runs import run_case
from scripts.stage4l.dense_probe import BASE, write


def worker(task):
    legacy, case, dt = task
    return str(run_case(legacy, case, dt, max_wall_s=7200))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--gate", type=Path, required=True)
    parser.add_argument(
        "--cases", nargs="+", choices=["B0", "B1", "B2", "B3"], default=["B0", "B1", "B2", "B3"]
    )
    parser.add_argument("--workers", type=int, choices=[1, 2, 3, 4], default=4)
    args = parser.parse_args()
    gate = json.loads(args.gate.read_text())
    if gate.get("dense_probe_gate") != "PASS_WITH_SCOPE":
        raise ValueError("Completed short-probe gate required")
    out = BASE / "primary_baseline" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_campaign")
    out.mkdir(parents=True, exist_ok=False)
    tasks = [(args.legacy_root, c, 1.0) for c in args.cases]
    tasks += [(args.legacy_root, c, 0.5) for c in args.cases if c in ["B2", "B3"]]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=get_context("spawn")) as pool:
        futures = {pool.submit(worker, t): t for t in tasks}
        for future in as_completed(futures):
            _, case, dt = futures[future]
            try:
                row = dict(case=case, dt_factor=dt, status="COMPLETE", run=future.result())
            except Exception as exc:
                row = dict(case=case, dt_factor=dt, status="FAILED", error=repr(exc))
            rows.append(row)
            write(out / "ROWS.json", rows)
            print(json.dumps(row), flush=True)
    write(
        out / "SUMMARY.json",
        dict(
            status="COMPLETE" if all(r["status"] == "COMPLETE" for r in rows) else "INCOMPLETE",
            rows=rows,
            gate_source=str(args.gate),
        ),
    )


if __name__ == "__main__":
    main()
