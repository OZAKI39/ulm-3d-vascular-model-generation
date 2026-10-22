"""Recompute archived primary nominal/dt-half comparisons without new transport."""

from pathlib import Path
import json
from scripts.stage4l.compare_branches import compare
from scripts.stage4l.dense_probe import sha, write

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"


def main():
    source = Path(
        "/workspace/mevo_acoustic_microbubbles_stage4l/results/stage4/lubrication_upgrade/comparison_legacy/20261005T232125840331Z_branches/SUMMARY.json"
    )
    expected = "92ce93854a392ccfc4beb878fa963303ade9a00f23ba8745383cf8906ad88c1e"
    if sha(source) != expected:
        raise ValueError("Archived branch report digest mismatch")
    saved = json.load(open(source))
    results = {}
    manifest = json.load(open(ROOT / "SYNC_MANIFEST.json"))
    included = {v["path"]: v for v in manifest["included_files"]}
    for case in ["B2", "B3"]:
        old = saved["timestep_sensitivity"][case]["comparison"]
        nominal = Path(old["reference"])
        half = Path(old["candidate"])
        hashes = []
        for folder in [nominal, half]:
            for name in ["RUN_IDENTITY.json", "BRANCH_METRICS.json", "metrics.json"]:
                p = folder / name
                relative = str(p).removeprefix("/workspace/mevo_acoustic_microbubbles_stage4l/")
                digest = sha(p)
                expected = included[relative]["sha256"]
                if digest != expected:
                    raise ValueError("Archived primary input digest mismatch")
                hashes.append(dict(path=str(p), sha256=digest))
        new = compare(nominal, half)
        results[case] = dict(status="PASS" if new == old else "FAIL", recomputed=new, input_hashes=hashes)
    write(
        BASE / "audit_v2/PRIMARY_TRANSPORT_REGRESSION.json",
        dict(
            status="PASS" if all(v["status"] == "PASS" for v in results.values()) else "FAIL",
            new_trajectories_run=False,
            scope="Archived Stage4L transport regression only; does not validate new Stage5A long transport",
            cases=results,
        ),
    )


if __name__ == "__main__":
    main()
