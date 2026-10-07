"""Describe actual overlapping prefixes without upgrading full-duration gates."""

import json
from pathlib import Path
import h5py
import numpy as np
from scripts.stage5a_r2.seal import BASE
from scripts.stage5a_r2.compare import compare
from scripts.stage4l.dense_probe import write, sha


def common_duration_and_identity(paths):
    with h5py.File(paths[0]) as h, h5py.File(paths[1]) as k:
        assert h["time_s"][0] == k["time_s"][0]
        for key in ["positions_m", "active"]:
            assert np.array_equal(h[key][0], k[key][0], equal_nan=True)
        for key in ["particle_id", "radii_m"]:
            assert np.array_equal(h[key][:], k[key][:])
        return float(min(h["time_s"][-1], k["time_s"][-1]) - h["time_s"][0])


def archived_science(campaign, source):
    """Existing complete frozen data can inform science, never fresh timing."""
    selection = json.loads((BASE / "fixtures/HOLDOUT_SELECTION.json").read_text())
    rows = []
    for case in ["B2", "B3"]:
        entries = {
            e["kind"]: e
            for e in campaign["ledger"]
            if e["case"] == case and e["duration"] == 0.005 and e.get("repeat") == 0 and e.get("summary")
        }
        if entries.get("reference", {}).get("summary", {}).get("status") == "COMPLETE":
            continue
        original = next(r for r in selection if r["case"] == case)
        reference = Path(original["source"])
        assert sha(reference) == original["source_sha256"]
        for kind in ["nominal", "refined"]:
            e = entries.get(kind)
            if not e or e["summary"]["status"] != "COMPLETE":
                continue
            paths = [reference, Path(e["run"]) / "trajectories.h5"]
            assert common_duration_and_identity(paths) >= 0.005 - 1e-12
            suffix = "legacy" if kind == "nominal" else "refined_legacy"
            r = compare(*paths, BASE / "comparisons", f"{case}_0.005_v4a_archived_{suffix}")
            r.update(
                reference_kind="ARCHIVED_STAGE4L",
                scope="Full 5 ms science against immutable archived Stage4L; not a fresh same-environment timing run",
                eligible_for_formal_campaign_gate=False,
                eligible_for_performance_gate=False,
                interpretation="Supplementary full-duration scientific comparison only. Original formal campaign gates remain unchanged; no historical wall time enters completion speedup.",
            )
            write(BASE / "comparisons" / f"{r['label']}.json", r)
            rows.append(r)
    write(
        BASE / "final/ARCHIVED_SCIENCE.json",
        dict(
            campaign_sha256=sha(source),
            rows=rows,
            affects_formal_gates=False,
            affects_performance=False,
            script_sha256=sha(__file__),
        ),
    )


def main():
    source = BASE / "final/CAMPAIGN_V4A/RESULT.json"
    campaign = json.loads(source.read_text())
    assert campaign["complete"]
    rows = []
    for case in ["B2", "B3"]:
        entries = {
            e["kind"]: e
            for e in campaign["ledger"]
            if e["case"] == case and e["duration"] == 0.005 and e.get("repeat") == 0 and e.get("summary")
        }
        for name, left, right in [
            ("legacy", "reference", "nominal"),
            ("self", "refined", "nominal"),
            ("refined_legacy", "reference", "refined"),
        ]:
            if left not in entries or right not in entries:
                continue
            a, b = entries[left], entries[right]
            if all(e["summary"]["status"] == "COMPLETE" for e in [a, b]):
                continue
            paths = [Path(e["run"]) / "trajectories.h5" for e in [a, b]]
            if common_duration_and_identity(paths) < 0.001 - 1e-12:
                continue
            r = compare(*paths, BASE / "comparisons", f"{case}_0.005_v4a_partial_{name}")
            r.update(
                scope="PARTIAL OVERLAP ONLY: requested 5 ms task not completed by both methods",
                target_duration_s=0.005,
                eligible_for_target_gate=False,
                reference_run_status=a["summary"]["status"],
                candidate_run_status=b["summary"]["status"],
                interpretation="PASS/FAIL below describes registered metrics over the actual common prefix only. It does not modify the campaign or establish complete 5 ms scientific validation.",
            )
            write(BASE / "comparisons" / f"{r['label']}.json", r)
            rows.append(r)
    write(
        BASE / "final/PARTIAL_COMPARISONS.json",
        dict(campaign_sha256=sha(source), rows=rows, affects_formal_gates=False, script_sha256=sha(__file__)),
    )
    print([(r["label"], r["compared_duration_s"], r["status"]) for r in rows], flush=True)
    archived_science(campaign, source)


if __name__ == "__main__":
    main()
