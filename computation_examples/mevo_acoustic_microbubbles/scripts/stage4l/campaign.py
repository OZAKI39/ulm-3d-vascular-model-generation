"""Finite dense probes only. No implicit launch of branch or full secondary runs."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from multiprocessing import get_context
from pathlib import Path

from scripts.stage4l.dense_probe import BASE, run, write


def worker(task):
    legacy, case, backend, duration, budget = task
    out, summary = run(legacy, case, backend, duration, max_wall_s=budget, max_steps=100000)
    return dict(case=case, backend=backend, run=str(out), status=summary["status"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--duration", type=float, choices=[0.001, 0.005], default=0.001)
    parser.add_argument(
        "--backends",
        nargs="+",
        default=["lubrication_coupled", "legacy_stokes"],
        choices=["lubrication_coupled", "legacy_stokes"],
    )
    parser.add_argument("--workers", type=int, choices=[1, 2], default=1)
    parser.add_argument("--max-wall-s", type=float, default=600.0)
    parser.add_argument("--gate", type=Path)
    args = parser.parse_args()
    if args.duration == 0.005:
        if (
            args.gate is None
            or json.loads(args.gate.read_text()).get("dense_probe_gate") != "PASS_WITH_SCOPE"
        ):
            raise ValueError("The paired 1 ms probe must pass before 5 ms")
    out = BASE / "dense_probe" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_campaign")
    out.mkdir(parents=True, exist_ok=False)
    tasks = [
        (args.legacy_root, case, backend, args.duration, args.max_wall_s)
        for backend in args.backends
        for case in ["B2", "B3"]
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=get_context("spawn")) as pool:
        futures = [pool.submit(worker, task) for task in tasks]
        for future in as_completed(futures):
            rows.append(future.result())
            write(out / "ROWS.json", rows)
    write(
        out / "SUMMARY.json",
        dict(status="COMPLETE" if all(r["status"] == "COMPLETE" for r in rows) else "INCOMPLETE", rows=rows),
    )


if __name__ == "__main__":
    main()
