"""Describe actual overlapping prefixes without upgrading full-duration gates."""

import json
from pathlib import Path
import h5py
import numpy as np
from scripts.stage5a_r2.seal import BASE
from scripts.stage5a_r2.compare import compare
from scripts.stage4l.dense_probe import write, sha


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
            with h5py.File(paths[0]) as h, h5py.File(paths[1]) as k:
                common = min(h["time_s"][-1], k["time_s"][-1]) - max(h["time_s"][0], k["time_s"][0])
                if common < 0.001 - 1e-12:
                    continue
                assert h["time_s"][0] == k["time_s"][0]
                for key in ["positions_m", "active"]:
                    assert np.array_equal(h[key][0], k[key][0], equal_nan=True)
                for key in ["particle_id", "radii_m"]:
                    assert np.array_equal(h[key][:], k[key][:])
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


if __name__ == "__main__":
    main()
