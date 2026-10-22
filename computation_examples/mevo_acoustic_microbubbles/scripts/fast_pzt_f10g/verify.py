"""Offline delivery checks. Never run scientific trajectories during verification."""

import gzip
import json
import subprocess
import numpy as np
import h5py
from .common import ROOT, BASE, BASELINE, CASES, js, sha, write


def main():
    checks = {}
    audit = js(BASE / "BASELINE_AUDIT.json")
    for p, h in audit["historical_sha256"].items():
        assert sha(ROOT / p) == h, p
    for p, v in audit["inputs"].items():
        assert sha(ROOT / "results/fast_pzt/f5" / p) == v["sha256"], p
    changes = subprocess.check_output(
        ["git", "diff", BASELINE, "--name-only"], cwd=ROOT, text=True
    ).splitlines()
    assert all("f10g" in p.lower() for p in changes), changes
    checks["historical_data_sources_and_inputs"] = "PASS"
    rm = js(BASE / "REMOTE_MANIFEST.json")
    for p, v in rm["files"].items():
        assert sha(BASE / p) == v["sha256"], p
    checks["remote_file_hashes"] = "PASS"
    small = js(BASE / "prototype_v1/SMALL_SYSTEMS.json")
    assert len(small["two"]) == 90 and len(small["four"]) == 18
    for a in small["two"] + small["four"]:
        result = a["result"]
        ref = a["reference"]
        assert (
            result["converged"]
            and np.max(np.abs(np.array(result["positions_m"]) - ref["positions_m"])) < 1e-12
        )
        assert min(result["pair_force_N"]) >= 0 and np.linalg.norm(result["net_contact_force_N"]) < 1e-20
    assert js(BASE / "checks_v1/RESULTS.json")["status"] == "PASS"
    checks["independent_small_systems_and_adapter"] = "PASS"
    static = js(BASE / "integration_v1/STATIC_FULL_GATE.json")
    assert static["SOFT_CONTACT_32_STATIC_PASS"] and all(
        r["quality"]["error_bits"] == 0 for r in static["rows"]
    )
    for case, expected in zip(CASES, [500, 235]):
        path = BASE / f"integration_v1/{case}_B2_K2"
        with gzip.open(path / "SOFT_STEPS.json.gz", "rt") as f:
            steps = json.load(f)
        assert len(steps) == expected and all(s["accepted"] and s["error_bits"] == 0 for s in steps)
        for i, s in enumerate(steps):
            r = s["result"]
            assert max(r["geometric_overlap_m"]) <= 17.8e-9 and r["max_single_iteration_m"] <= 178e-9
            assert r["converged"] and r["force_position_consistency_m"] < 1e-13
            assert np.allclose(s["actual_accepted_positions_m"], r["positions_m"], rtol=0, atol=1e-15)
            if i:
                assert np.array_equal(s["initial_positions_m"], steps[i - 1]["actual_accepted_positions_m"])
        with h5py.File(path / "trajectory.h5") as h:
            assert np.array_equal(h["particle_id"][:], np.arange(32))
            assert np.all(h["radii_m"][:] == 1.78e-6)
            assert len(h["substeps/accepted_end_s"]) == expected
            for i, time in enumerate(h["time_s"][1:], 1):
                np.testing.assert_array_equal(
                    h["positions_m"][i, 0], steps[round(time / 1e-5) - 1]["actual_accepted_positions_m"]
                )
    checks["735_accepted_steps_and_recording"] = "PASS"
    with gzip.open(BASE / "checks_v1/ENGINEERING_EVIDENCE.json.gz", "rt") as f:
        r = json.load(f)
    assert r["accepted_clock"] == 0 and r["rejected_step"]["error_bits"] != 0
    np.testing.assert_array_equal(
        r["rejected_step"]["initial_positions_m"], r["rejected_step"]["actual_accepted_positions_m"]
    )
    checks["retained_raw_rejection"] = "PASS"
    gate = js(BASE / "analysis/FINAL_GATE.json")
    assert (
        gate["PRODUCTION_NOT_ADMITTED"]
        and not gate["training_ready"]
        and gate["GEOMETRIC_SPRING_LAW_EXACT"] is False
    )
    assert all(not js(p)["exit_code"] for p in BASE.glob("*_PROCESS.json"))
    if (BASE / "figures/VISUAL_REVIEW.json").exists():
        review = js(BASE / "figures/VISUAL_REVIEW.json")
        assert review["status"] == "PASS"
        for p, h in review["files_sha256"].items():
            assert sha(BASE / "figures" / p) == h, p
        checks["reviewed_figures"] = "PASS"
    else:
        checks["reviewed_figures"] = "PENDING"
    if (BASE / "progress_document/VERIFICATION.json").exists():
        v = js(BASE / "progress_document/VERIFICATION.json")
        assert v["status"] == "PASS"
        assert sha("/home/lzy/projects/microbubble_aggregation.md") == v["after_sha256"]
        assert sha(v["backup"]) == v["before_sha256"]
        checks["progress_document"] = "PASS"
    else:
        checks["progress_document"] = "PENDING"
    write(
        BASE / "analysis/DELIVERY_VERIFICATION.json",
        dict(status="PASS" if all(v == "PASS" for v in checks.values()) else "PENDING", checks=checks),
    )
    print(checks)


if __name__ == "__main__":
    main()
