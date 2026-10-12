"""Final engineering/data-integrity gates, distinct from scientific success."""

from pathlib import Path
import json
import subprocess
import hashlib
from datetime import datetime, timezone
import h5py
import numpy as np
from .common import ROOT, BASE, F10A, js, write, sha


def main():
    baseline = "adeceb1d3b8ce82312e338bfe5a8b7d15f20b87f"
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    repo = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip())
    diff = subprocess.check_output(["git", "diff", baseline, "--name-status"], cwd=ROOT, text=True)
    lines = diff.splitlines()
    assert all(line.startswith("A\t") for line in lines), lines
    for rel, meta in js(F10A / "ARTIFACT_MANIFEST.json")["files"].items():
        assert sha(F10A / rel) == meta["sha256"], rel
    write(
        BASE / "HISTORY_PRESERVATION.json",
        dict(
            status="PASS",
            baseline_commit=baseline,
            tracked_changes=lines,
            existing_files_modified=0,
            F10A_artifacts_verified=228,
        ),
    )
    checked = set()
    for path in BASE.glob("*/IDENTITY.json"):
        data = js(path)
        commit = data["execution_commit"]
        subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=ROOT, check=True)
        for rel, digest in data["source_sha256"].items():
            key = (commit, rel, digest)
            if key in checked:
                continue
            blob = subprocess.check_output(["git", "show", commit + ":" + prefix + rel], cwd=repo)
            assert hashlib.sha256(blob).hexdigest() == digest, (commit, rel)
            checked.add(key)
    write(
        BASE / "SOURCE_VERIFICATION.json",
        dict(
            status="PASS",
            unique_commit_file_hashes=len(checked),
            science_jobs=[str(p.relative_to(BASE)) for p in BASE.glob("*/IDENTITY.json")],
        ),
    )
    remote = js(BASE / "REMOTE_MANIFEST.json")
    for rel, meta in remote["files"].items():
        assert sha(BASE / rel) == meta["sha256"], rel
    result = js(BASE / "analysis/RESULTS.json")
    assert result["models"]["G1"]["surrogate_numerical_admission"] == "PASS"
    assert result["models"]["G0"]["surrogate_numerical_admission"] == "FAIL"
    assert result["models"]["G1"]["performance"]["body_gate"] == "FAIL"
    checks = js(BASE / "checks_v2/RESULTS.json")
    assert checks["status"] == "PASS"
    assert all(s["status"] == "SURROGATE_NOT_ADMITTED" for s in checks["scope"][3:])
    assert js(BASE / "support_v1/RESULTS.json")["status"] == "PASS"
    for mode in ["G0", "G1"]:
        for key in ["OFF", "S_1000", "O_1000"]:
            path = BASE / "rollout_v1" / (mode + "_" + key)
            cfg = js(path / "CONFIG.json")
            original = js(F10A / "extended_10_v1" / key / "CONFIG.json")
            np.testing.assert_array_equal(cfg["initial_positions_m"], original["environment"]["positions_m"])
            assert cfg["event"] == original["events"][0]
            with h5py.File(path / "centroid.h5") as h:
                assert set(h.keys()) == {"time_s", "centroid_m"}
                assert h["centroid_m"].shape == (201, 3) and np.isfinite(h["centroid_m"][:]).all()
                assert abs(h["time_s"][-1] - 0.01) < 1e-14
                meta = json.loads(h.attrs["metadata_json"])
                assert meta["training_ready"] is False and meta["quality"]["error_bits"] == 0
    review = js(BASE / "figures/VISUAL_REVIEW.json")
    assert review["status"] == "PASS"
    for rel, digest in review["files"].items():
        assert sha(BASE / "figures" / rel) == digest, rel
    doc = js(BASE / "progress_document/VERIFICATION.json")
    assert doc["status"] == "PASS"
    assert sha(doc["backup"]) == doc["before_sha256"] and sha(doc["original_path"]) == doc["after_sha256"]
    gate = dict(
        status="COMPLETE_WITH_SCIENTIFIC_LIMITATIONS",
        timestamp_UTC=datetime.now(timezone.utc).isoformat(),
        implementation_and_checks="PASS",
        contact_diagnosis="COMPLETE_A0_A1_A2_REJECT_A3_SHORT_PASS",
        surrogate_numerical_admission={"G0": "FAIL", "G1": "PASS"},
        standalone_2x_engineering_value="FAIL",
        physics_experimentally_calibrated=False,
        particle_morphology_validated=False,
        training_ready=False,
        production_bidirectional_switching="NOT_IMPLEMENTED_NOT_ADMITTED",
        archive="PASS",
        source_identity="PASS",
        historical_preservation="PASS",
        figures="PASS",
        progress_document="PASS",
        scope="One mature32-bubble template,OFF/S/O,10ms,originalCFD,M0; all requested bounded tasks complete, neither physical calibration nor2x performance established",
    )
    write(BASE / "analysis/FINAL_GATE.json", gate)
    files = {
        str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
        for p in sorted(BASE.rglob("*"))
        if p.is_file() and p.name != "ARTIFACT_MANIFEST.json"
    }
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(files=files, count=len(files), bytes=sum(m["bytes"] for m in files.values())),
    )
    print(
        json.dumps(
            dict(
                status=gate["status"],
                files=len(files),
                source_hashes=len(checked),
                bytes=sum(m["bytes"] for m in files.values()),
            )
        )
    )


if __name__ == "__main__":
    main()
