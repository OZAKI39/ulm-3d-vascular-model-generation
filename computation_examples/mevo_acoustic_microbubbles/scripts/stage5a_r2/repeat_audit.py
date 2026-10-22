"""Post-timing scientific repeatability with unchanged registered metrics.

Adaptive grids can differ across repeated float64 runs. Never infer scientific
equivalence from the same command, accepted count, or similar wall time alone.
All repeats remain in the timing ledger; this audit cannot upgrade campaign gates.
"""

import json
from pathlib import Path
from scripts.stage5a_r2.seal import BASE
from scripts.stage5a_r2.compare import compare
from scripts.stage5a_r2.partial_comparisons import common_duration_and_identity
from scripts.stage4l.dense_probe import write, sha


def main():
    source = BASE / "final/CAMPAIGN_V4A/RESULT.json"
    campaign = json.loads(source.read_text())
    assert campaign["complete"], "Run only after isolated timing ends"
    entries = [
        e
        for e in campaign["ledger"]
        if e.get("summary") and e["summary"]["status"] == "COMPLETE" and e["kind"] in ["nominal", "reference"]
    ]
    rows = []
    for e in entries:
        if e["repeat"] == 0:
            continue
        kinds = [e["kind"]] + (["reference"] if e["kind"] == "nominal" else [])
        for kind in kinds:
            anchor = next(
                q
                for q in entries
                if q["case"] == e["case"]
                and q["duration"] == e["duration"]
                and q["kind"] == kind
                and q["repeat"] == 0
            )
            paths = [Path(q["run"]) / "trajectories.h5" for q in [anchor, e]]
            assert common_duration_and_identity(paths) >= e["duration"] - 1e-12
            label = f"{e['case']}_{e['duration']:g}_v4a_repeat{e['repeat']}_{e['kind']}_vs_{kind}0"
            r = compare(*paths, BASE / "comparisons", label)
            r.update(
                case=e["case"],
                repeat=e["repeat"],
                kind=e["kind"],
                anchor_kind=kind,
                candidate_accepted_steps=e["summary"]["accepted_steps"],
                anchor_accepted_steps=anchor["summary"]["accepted_steps"],
                role="Additional post-timing audit under unchanged registered scientific thresholds; no parameter tuning or timing-sample removal",
            )
            write(BASE / "comparisons" / f"{label}.json", r)
            rows.append(r)
    status = (
        (
            "PASS"
            if all(r["pass_scientific"] for r in rows)
            else "FAIL"
            if any(r["status"] == "FAIL" for r in rows)
            else "INCOMPLETE"
        )
        if rows
        else "NOT_RUN"
    )
    write(
        BASE / "final/REPEAT_AUDIT.json",
        dict(
            status=status,
            rows=rows,
            campaign_sha256=sha(source),
            script_sha256=sha(__file__),
            timing_samples_removed=0,
            scope="Completed equal-duration repetitions, identical original IDs/initial coordinates, actual accepted states and sampling uncertainty. Raw campaign results remain unchanged.",
        ),
    )
    print([(r["label"], r["status"], r["position_endpoint_max_difference_m"]) for r in rows], flush=True)


if __name__ == "__main__":
    main()
