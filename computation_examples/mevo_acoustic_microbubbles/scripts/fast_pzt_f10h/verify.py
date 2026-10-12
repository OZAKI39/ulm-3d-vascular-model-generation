"""Read-only offline scientific verification plus separate derived admissions."""

import subprocess
import numpy as np
import h5py
from .common import ROOT, BASE, BASELINE, js, sha, write
from .analyze import gz, read_steps


def main():
    checks = {}
    old = ROOT / "results/fast_pzt/f10g"
    for name, v in js(old / "ARTIFACT_MANIFEST.json")["files"].items():
        assert sha(old / name) == v["sha256"], name
    for name, v in js(BASE / "BASELINE_AUDIT.json")["inputs"].items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == v["sha256"], name
    changes = subprocess.check_output(
        ["git", "diff", BASELINE, "--name-status"], cwd=ROOT, text=True
    ).splitlines()
    assert all(s.startswith("A\t") and "f10h" in s.lower() for s in changes), changes
    checks["F10G_historical_bytes_inputs_default_preserved"] = "PASS"
    for name, v in js(BASE / "REMOTE_MANIFEST.json")["files"].items():
        assert sha(BASE / name) == v["sha256"], name
    assert js(BASE / "REMOTE_VERIFICATION.json")["status"] == "PASS"
    for p in BASE.glob("*_PROCESS.json"):
        assert js(p)["exit_code"] == 0, p
    for label in ["checks_v1", "checks_v2"]:
        assert js(BASE / label / "RESULTS.json")["status"] == "PASS"
    checks["fixed_sources_remote_hashes_and_engineering"] = "PASS"
    small = js(BASE / "geometry_v1/SMALL_SYSTEMS.json")
    assert len(small["rows"]) == 18 and small["two_pass"] and small["four_pass"]
    maxforce = 0.0
    for r in small["rows"]:
        z = r["result"]
        ref = r["reference"]
        assert z["converged"] and ref["reference_pass"] and r["reference_error_m"] < 1e-12
        error = float(np.max(abs(np.array(z["contact_force_N"]) - ref["contact_force_N"])))
        maxforce = max(maxforce, error)
        assert error < 1e-15
        assert np.min(z["pair_force_N"]) >= 0 and np.linalg.norm(z["net_contact_force_N"]) < 1e-20
    assert js(BASE / "geometry_v1/DENSE_STATIC.json")["both_pass"]
    assert js(BASE / "b2_short_v1/STATIC_FULL_GATE.json")["status"] == "PASS"
    checks["independent_two_four_static32_and_force_reference"] = "PASS"
    trials = (
        js(BASE / "analysis/LUBRICATION_TRAJECTORIES.json")["trials"]
        + js(BASE / "analysis/B2_GEOMETRY_TRAJECTORIES.json")["trials"]
    )
    total_accepted = 0
    failures = 0
    admissions = BASE / "analysis/admissions"
    admissions.mkdir(exist_ok=True)
    for r in trials:
        path = ROOT / r["path"]
        times, x, steps = read_steps(path)
        ledger = gz(path / "SUBSTEP_LEDGER.json.gz")["steps"]
        mode = r["mode"]
        dt = 5e-6 if mode.endswith("DT5") else 1e-5
        n = sum(s["accepted"] for s in ledger)
        total_accepted += n
        assert n == r["accepted_steps"] == len(times) - 1
        assert np.allclose(times, np.arange(n + 1) * dt, rtol=0, atol=1e-14)
        for i, s in enumerate(steps):
            np.testing.assert_array_equal(s["initial_positions_m"], x[min(i, n)])
            if s["accepted"]:
                assert not s["error_bits"]
                np.testing.assert_array_equal(s["positions_after_step_m"], x[i + 1])
                d = np.linalg.norm(x[i + 1, :, None, :] - x[i + 1, None, :, :], axis=-1)
                np.fill_diagonal(d, np.inf)
                assert max(0, 3.56e-6 - d.min()) <= 17.8e-9
            else:
                failures += 1
                assert i == len(steps) - 1
                np.testing.assert_array_equal(s["positions_after_step_m"], s["initial_positions_m"])
            if mode.startswith("M0_"):
                assert len(s["lubrication"]) == 4
        for s in ledger:
            if s["accepted"]:
                assert (
                    s["final_overlap_m"] <= 17.8e-9 and max(s["round_max_correction_m"], default=0) <= 178e-9
                )
        with h5py.File(path / "trajectory.h5") as h:
            np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
            assert np.all(h["radii_m"][:] == 1.78e-6)
            assert not h.attrs["training_ready"] and len(h["substeps/accepted_end_s"]) == n
            for t, v in zip(h["time_s"][:], h["positions_m"][:, 0]):
                np.testing.assert_array_equal(v, x[round(t / dt)])
        if mode.startswith("B2"):
            soft = gz(path / "SOFT_STEPS.json.gz")
            assert len(soft) == n
            for s in soft:
                z = s["result"]
                assert (
                    z["converged"]
                    and z["force_position_consistency_m"] <= 1.01e-13
                    and len(z["iterations"]) <= 24
                )
                assert z["kc_N_m"] == 0.04641319518378658 and z["dt_s"] == dt
                gap = np.array(z["geometric_gap_m"])
                reaction = np.array(z["pair_force_N"])
                assert np.all(reaction[gap >= 0] == 0)
                np.testing.assert_array_equal(reaction, z["kc_N_m"] * np.maximum(-gap, 0))
            shape_status = "POSITIONS_SHAPE_PASS_EDGE_TIMING_FAIL"
        else:
            shape_status = "MORPHOLOGY_NOT_CONVERGED"
        raw = js(path / "NUMERICAL_IDENTITY.json")
        assert raw["numerical_scheme_id"] == mode and not raw["complete_trajectory"]
        parents = js(ROOT / "results/fast_pzt/f10e/SELECTION.json")
        parent = parents["branch"] if r["case"] == "O_1000_1102" else parents["preparation"]
        assert raw["split_group"] == js(ROOT / parent["relative_path"] / "SUMMARY.json")["split_group"]
        write(
            admissions / (path.parent.name + "__" + path.name + ".json"),
            dict(
                source_path=r["path"],
                source_sha256=sha(path / "trajectory.h5"),
                numerical_scheme_id=mode,
                morphology_admission=shape_status,
                timestep_convergence="FAIL_FULL_GATE",
                timestep_tested_end_s=0.005 if r["case"] == "O_1000_1102" else 0.00235,
                complete_trajectory=False,
                contact_quality="SHORT_WINDOW_PASS_LONG_ATTEMPT_REJECTED"
                if r["status"] == "FAILED"
                else "REQUESTED_SHORT_WINDOW_PASS",
                failure_reason=raw["failure_reason"],
                training_ready=False,
                scope="Derived from actual comparisons; original run-time PENDING fields and raw HDF5 are preserved. No convergence inference beyond tested window.",
            ),
        )
    assert total_accepted == 4341 and failures == 1
    checks["actual_physical_steps_budgets_clocks_ids_split_groups"] = "PASS"
    failure = js(BASE / "b2_short_v1/FAIL_STOP_CHECK.json")
    assert failure["first_rejection"]["error_bits"] and not failure["step"]["accepted"]
    assert js(BASE / "analysis/EDGE_AUDIT.json")["status"] == "PASS_EXACT_GPU_EVENT_RECONSTRUCTION"
    checks["independent_all_accepted_edge_reconstruction"] = "PASS"
    gate = js(BASE / "analysis/FINAL_GATE.json")
    assert gate["PRODUCTION_NOT_ADMITTED"] and not gate["training_ready"]
    review = BASE / "figures/VISUAL_REVIEW.json"
    if review.exists():
        v = js(review)
        assert v["status"] == "PASS"
        for name, h in v["files_sha256"].items():
            assert sha(BASE / "figures" / name) == h, name
        checks["actual_PNG_PDF_visual_review"] = "PASS"
    else:
        checks["actual_PNG_PDF_visual_review"] = "PENDING"
    progress = BASE / "progress_document/VERIFICATION.json"
    if progress.exists():
        v = js(progress)
        assert v["status"] == "PASS" and v["original_bytes_preserved"] and v["headings_preserved"]
        assert sha("/home/lzy/projects/microbubble_aggregation.md") == v["after_sha256"]
        assert sha(v["backup"]) == v["before_sha256"]
        checks["progress_backup_bytes_headings"] = "PASS"
    else:
        checks["progress_backup_bytes_headings"] = "PENDING"
    write(
        BASE / "analysis/DELIVERY_VERIFICATION.json",
        dict(
            status="PASS" if all(v == "PASS" for v in checks.values()) else "PENDING",
            checks=checks,
            science_accepted_steps=total_accepted,
            science_rejections=failures,
            small_force_reference_max_error_N=maxforce,
        ),
    )
    print(checks)


if __name__ == "__main__":
    main()
