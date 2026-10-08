"""Offline closure checks on actual F3 evidence; never advances or edits trajectories."""

import hashlib
import subprocess
import numpy as np
import h5py
from mevo_mb.fast_sim.kinematics import affine_velocity
from scripts.fast_pzt_f2.analyze import read
from .common import ROOT, BASE, F2, PROTOCOL, js, write, state


def main():
    rng = np.random.default_rng(20261008)
    x = rng.normal(size=(32, 3)) * 4e-6
    L = np.array([[2, -3, 1], [5, -1, 2], [0.5, 3, 4]])
    u = x @ L.T + np.array([0.001, -0.002, 0.003])
    fit = affine_velocity(x, u)
    np.testing.assert_allclose(fit["L_per_s"], L, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(fit["E_per_s"], (L + L.T) / 2, atol=1e-12)
    np.testing.assert_allclose(fit["W_per_s"], (L - L.T) / 2, atol=1e-12)
    moved = affine_velocity(x + [0.001, -0.002, 0.003], u + [0.01, 0.01, -0.01])
    np.testing.assert_allclose(moved["L_per_s"], L, rtol=1e-11, atol=1e-11)
    degenerate = affine_velocity(np.column_stack([np.arange(32), np.zeros((32, 2))]), u)
    assert degenerate["status"] == "DEGENERATE" and "L_per_s" not in degenerate

    rows = {}
    mature = state()
    reference = read(BASE / "repeats_v1/stable_reduction_0")
    for hpath in sorted(BASE.glob("*/*/trajectory.h5")):
        d = read(hpath.parent)
        s = js(hpath.parent / "SUMMARY.json")
        cfg = js(hpath.parent / "CONFIG.json")
        assert s["status"] == "COMPLETE" and s["numerical_quality"] == "PASS"
        assert s["actual_substep_count"] == s["macro_step_count"]
        assert s["nominal_dt_s"] == PROTOCOL["dt_s"] and s["substeps_per_macro"] == 1
        assert s["quality"]["error_bits"] == 0
        assert s["quality"]["max_accepted_penetration_m"] <= 0.01 * 1.78e-6
        assert s["quality"]["max_contact_pass_correction_m"] <= 0.1 * 1.78e-6
        assert np.isfinite(d["positions_m"]).all()
        np.testing.assert_array_equal(d["particle_id"], np.arange(d["positions_m"].shape[2]))
        np.testing.assert_array_equal(d["radii_m"], np.full(d["positions_m"].shape[2], 1.78e-6))
        np.testing.assert_allclose(np.diff(d["time_s"]), cfg["sample_interval_s"], atol=2e-14, rtol=0)
        np.testing.assert_allclose(d["time_s"][-1], cfg["duration_s"], atol=2e-14, rtol=0)
        if d["positions_m"].shape[2] == 32:
            assert d["positions_m"].shape == (601, 1, 32, 3)
            if not hpath.parent.name.startswith("initial_perturb"):
                np.testing.assert_array_equal(d["positions_m"][0], mature)
            if hpath.parent.name in ["stable_reduction_1", "stable_reduction_2"]:
                for key in ["positions_m", "time_s", "edge_times", "edge_ids"]:
                    np.testing.assert_array_equal(d[key], reference[key])
        rows[str(hpath.parent.relative_to(BASE))] = dict(
            frames=len(d["time_s"]),
            count=d["positions_m"].shape[2],
            quality=s["numerical_quality"],
            steps=s["actual_substep_count"],
            max_accepted_overlap_nm=s["quality"]["max_accepted_penetration_m"] * 1e9,
            max_contact_pass_correction_nm=s["quality"]["max_contact_pass_correction_m"] * 1e9,
        )
    assert len(rows) == 12
    driving = js(BASE / "driving_v3/RESULT.json")
    for row in driving["candidate_rows"].values():
        assert row["field"]["derivative_step_halving_relative_L2"] < 1e-8
        assert row["coupled_pressure_quality"]["error_bits"] == 0
    with h5py.File(BASE / "driving_v3/FORCE_FIELD_DATA.h5") as h:
        np.testing.assert_array_equal(h["positions_m"][:], mature[0])
        for name in driving["candidate_rows"]:
            np.testing.assert_allclose(
                h[name + "/total_N"][:],
                h[name + "/primary_N"][:] + h[name + "/secondary_N"][:],
                atol=1e-24,
                rtol=1e-12,
            )
            assert np.linalg.norm(h[name + "/secondary_N"][:].sum(0)) < 1e-23
    # Compare every F1/F2 artifact listed in the earlier manifest, not only Git-tracked files.
    preserved = {}
    for stage in ["f1", "f2"]:
        folder = F2.parent / stage
        manifest = js(folder / "ARTIFACT_MANIFEST.json")
        checked = 0
        absent = []
        for entry in manifest["artifacts"]:
            path = folder / entry["relative_path"]
            if not path.is_file():
                absent.append(entry["relative_path"])
                continue
            assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"], str(path)
            checked += 1
        preserved[stage] = dict(
            local_checked=checked,
            local_not_present=absent,
            note="Server originals are checked separately, including excluded large files",
        )
    old_paths = [
        "src",
        "scripts/fast_pzt",
        "scripts/fast_pzt_f2",
        "configs/fast_pzt/f1",
        "configs/fast_pzt/f2",
        "results/fast_pzt/f1",
        "results/fast_pzt/f2",
    ]
    changes = subprocess.check_output(
        ["git", "diff", "--name-status", PROTOCOL["f2_archive"], "--", *old_paths], cwd=ROOT, text=True
    )
    unexpected = [
        line
        for line in changes.splitlines()
        if not (
            line.startswith("A\t")
            and "/fast_sim/" in line
            and any(n in line for n in ["f3_engine.py", "f3_kernels.py", "kinematics.py"])
        )
    ]
    assert not unexpected, unexpected
    write(
        BASE / "ANALYSIS_CHECKS.json",
        dict(
            status="PASS",
            trajectories=rows,
            checks=[
                "Known affine velocity gradient and E/W decomposition",
                "Translation-invariant affine fit",
                "Rank-degenerate fit withheld",
                "All 12 saved trajectories satisfy clocks/IDs/radii/quality",
                "Stable saved raw trajectories and events repeat exactly",
                "Same mature positions before control branches",
                "Finite-difference step-halving convergence",
                "Force decomposition and pair action/reaction",
            ],
            preservation=preserved,
            inherited_source_diff=changes,
        ),
    )
    print("PASS: F3 offline evidence, kinematics and preserved-source checks")


if __name__ == "__main__":
    main()
