"""Verify only the frozen inputs and source boundary of the current M0 runs."""

import json
from pathlib import Path
from scripts.stage4m.common import BASE, refs, write


def audit(legacy, stage4l):
    frozen = json.loads((BASE / "existing_reference/BACKEND_IDENTITY.json").read_text())[
        "stage4l_current_source_sha256"
    ]
    rows = []
    for identity_path in sorted((BASE / "dense_probe").glob("*/*/IDENTITY.json")):
        d = json.loads(identity_path.read_text())
        a0 = refs(legacy, stage4l, d["case"])["dense_A0"]
        old = json.loads((a0 / "IDENTITY.json").read_text())
        mismatch = {
            p: {"A0": v, "new": d["input_sha256"].get(p)}
            for p, v in old["input_sha256"].items()
            if d["input_sha256"].get(p) != v
        }
        changed = {
            p: {"stage4l": v, "new": d["source_sha256"].get(p)}
            for p, v in frozen.items()
            if d["source_sha256"].get(p) != v
        }
        prohibited = set(changed) - {"src/mevo_mb/bubbles/forces.py", "src/mevo_mb/bubbles/tracker.py"}
        rows.append(
            dict(
                run=str(identity_path.parent),
                policy=d["policy"],
                case=d["case"],
                existing_A0_input_hashes_all_match=not mismatch,
                input_mismatches=mismatch,
                changed_existing_production_files=changed,
                prohibited_changes=sorted(prohibited),
                passed=not mismatch and not prohibited,
            )
        )
    result = dict(
        passed=bool(rows) and all(r["passed"] for r in rows),
        rows=rows,
        scope="Read-only input hashes and frozen source boundary; initial checkpoint arrays were separately checked bitwise by every M0 runner.",
    )
    write(BASE / "existing_reference/M0_INPUT_AND_SOURCE_AUDIT.json", result)
    print(json.dumps(dict(passed=result["passed"], checked_runs=len(rows))), flush=True)
    if not result["passed"]:
        raise RuntimeError("Frozen M0 input/source mismatch")
    return result


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--legacy-root", type=Path, required=True)
    p.add_argument("--stage4l-root", type=Path, required=True)
    a = p.parse_args()
    audit(a.legacy_root, a.stage4l_root)
