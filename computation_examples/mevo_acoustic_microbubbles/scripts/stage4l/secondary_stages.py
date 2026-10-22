"""Finite secondary stages, never a full 6 s submission."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from scripts.stage4l.branch_runs import run_case
from scripts.stage4l.dense_probe import BASE, write


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--legacy-root", type=Path, required=True)
    p.add_argument("--gate", type=Path, required=True)
    p.add_argument("--duration", type=float, choices=[0.1, 0.5, 1.0], required=True)
    args = p.parse_args()
    gate = json.loads(args.gate.read_text())
    if gate.get("dense_probe_gate") != "PASS_WITH_SCOPE":
        raise ValueError("Dense probe gate required")
    if args.duration > 0.1 and gate.get("secondary_0p1_gate") != "PASS_WITH_SCOPE":
        raise ValueError("The 0.1 s secondary stage must pass before this stage")
    out = BASE / "secondary_short" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_campaign")
    out.mkdir(exist_ok=False, parents=True)
    rows = []
    for case in ["B2", "B3"]:
        try:
            folder = run_case(args.legacy_root, case, secondary=True, duration=args.duration, max_wall_s=1800)
            row = dict(case=case, run=str(folder), status="COMPLETE")
        except Exception as exc:
            row = dict(case=case, status="FAILED", error=repr(exc))
        rows.append(row)
        write(out / "ROWS.json", rows)
        print(json.dumps(row), flush=True)
    write(
        out / "SUMMARY.json",
        dict(
            status="COMPLETE" if all(r["status"] == "COMPLETE" for r in rows) else "INCOMPLETE",
            duration_s=args.duration,
            rows=rows,
        ),
    )


if __name__ == "__main__":
    main()
