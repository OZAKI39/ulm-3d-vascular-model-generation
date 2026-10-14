"""Independent archive, geometry, force mapping, controls and admission checks."""

import time
import subprocess
import numpy as np
import h5py
from scipy.spatial.distance import pdist
from mevo_mb.fast_sim.f10g_contact import resistance
from .common import BASE, OLD, ROOT, CONFIG, read, write, sha, identity


def verify(include_visual=True):
    start = time.perf_counter()
    p = read(CONFIG)
    checks = {}
    evidence = {}
    protected = read(BASE / "cleanup/AUDIT.json")["protected_sha256"]
    checks["protected_F10P_B2_M0_R3_fields"] = all(sha(ROOT / k) == v for k, v in protected.items())
    oldp = read(ROOT / "configs/straight_tube_exp1/f10p/protocol.json")
    thresholds = [
        "contact_overlap_m",
        "hard_geometry_m",
        "contact_correction_m",
        "mobility_residual_m",
        "primal_violation_m",
        "dual_negative_N",
        "force_residual_atol_N",
        "force_residual_rtol",
        "complementarity_equivalent_gap_m",
        "actual_geometry_complementarity_m",
        "group_COM_tolerance_m",
        "group_Rg_tolerance_m",
        "group_Q_relative",
        "group_xy_tolerance",
        "group_aspect_relative",
        "principal_axis_angle_deg",
        "distribution_RMS_over_R",
        "fixed_id_tolerance_m",
        "XY_targets_m",
        "XY_tolerance_m",
        "persistence_s",
        "instant_edge_time_tolerance_s",
    ]
    checks["all_inherited_thresholds_unchanged"] = all(p[k] == oldp[k] for k in thresholds)
    checks["input_hashes"] = all(
        sha(OLD / k) == v for k, v in read(BASE / "PREREGISTRATION.json")["input_sha256"].items()
    )
    manifest = read(BASE / "REMOTE_MANIFEST.json")
    checks["all_remote_sha256"] = all(sha(BASE / k) == r["sha256"] for k, r in manifest["files"].items())
    evidence["returned_files"] = len(manifest["files"])
    evidence["returned_bytes"] = manifest["bytes"]
    checks["frozen_and_small_independent_refs"] = read(BASE / "local/FOUNDATIONS_v2.json")["status"] == "PASS"
    checks["first_three_QP_failure_retained"] = read(BASE / "local/FOUNDATIONS_v1.json")["status"] == "FAIL"
    checks["sequential_cap_five"] = (
        read(CONFIG.parent / "algorithm_revision_1.json")["max_sequential_qp"] == 5
    )
    audit = read(BASE / "analysis/TRUE_GAP_FAILURE_DIAGNOSIS.json")
    checks["weighted_comp_reproduces"] = (
        abs(audit["recomputed_true_comp_m"] - audit["saved_true_comp_m"]) < 1e-22
    )
    violations = []
    max_recomputed_force_error = 0.0
    max_saved_force_residual_difference = 0.0
    count = 0
    for rp in sorted(BASE.glob("*/runs/*/RESULT.json")):
        r = read(rp)
        inp = read(rp.parent / "INPUT.json")
        dt = inp["dt_s"]
        R = p["radius_m"]
        zeta = 6 * np.pi * 0.001 * R
        with h5py.File(rp.parent / "trajectory.h5") as h:
            x = h["positions_m"][:]
            t = h["time_s"][:]
            ids = h["particle_id"][:]
            rr = h["step_metrics"][:]
            accepted = rr[:, 0].astype(bool)
            ee = h["extra_metrics"][:]
            lam = h["pair_force_N"][:]
            mask = rr[accepted]
            assert np.array_equal(ids, np.arange(32))
            assert np.array_equal(x[0], np.load(OLD / "initial_states.npz")["aggregate_m"])
            assert h.attrs["flow"] == "STRICT_ZERO"
            assert len(x) == 1 + int(accepted.sum())
            assert np.allclose(t, np.arange(len(t)) * dt, rtol=0, atol=1e-15)
            overlaps = np.array([max(0, -(pdist(y) - 2 * R).min()) for y in x[1:]])
            assert np.max(abs(overlaps - mask[:, 2]), initial=0) < 1e-20
            assert np.all(overlaps <= p["hard_geometry_m"]) and np.all(overlaps <= p["contact_overlap_m"])
            assert (
                np.all(mask[:, 1] <= p["mobility_residual_m"])
                and np.all(mask[:, 4] <= p["contact_correction_m"])
                and np.all(mask[:, 5] == 0)
            )
            assert (
                np.all(mask[:, 12] <= p["primal_violation_m"])
                and np.all(mask[:, 13] <= p["dual_negative_N"])
                and np.all(mask[:, 15] <= p["complementarity_equivalent_gap_m"])
                and np.all(mask[:, 17] <= p["actual_geometry_complementarity_m"])
            )
            force_scale = np.maximum.reduce([ee[:, 5], abs(lam).max(1), np.full(len(rr), 1e-18)])
            assert np.all(
                mask[:, 14] <= p["force_residual_atol_N"] + p["force_residual_rtol"] * force_scale[accepted]
            )
            assert np.all((ee[:, 3] == 0) | (ee[:, 3] >= 0.1)) and np.all(ee[:, 4] <= 40000)
            assert np.max(rr[:, 9]) <= 5
            for field_time, F in zip(h["force_samples/time_s"][:], h["force_samples/total_force_N"][:]):
                j = round(field_time / dt)
                if j + 1 >= len(x):
                    continue
                i, k = np.triu_indices(32, 1)
                dy = x[j + 1, k] - x[j + 1, i]
                normal = dy / np.linalg.norm(dy, axis=1)[:, None]
                contact = np.zeros((32, 3))
                np.add.at(contact, i, -normal * lam[j, :, None])
                np.add.at(contact, k, normal * lam[j, :, None])
                A = resistance(x[j], R)["A"]
                v = (x[j + 1] - x[j]) / dt
                err = float(np.linalg.norm((zeta * A @ v.ravel()).reshape(-1, 3) - F - contact, axis=1).max())
                max_recomputed_force_error = max(max_recomputed_force_error, err)
                max_saved_force_residual_difference = max(
                    max_saved_force_residual_difference, abs(err - rr[j, 14])
                )
                count += 1
            if r["status"] == "FAIL":
                assert "failure_proposal_m" in h and t[-1] < r["requested_s"]
                violations.append(dict(name=r["name"], legal_s=t[-1], reasons=r["reasons"]))
    checks["accepted_states_geometry_pressure_correction_balance"] = True
    checks["independent_true_normal_reaction_reconstruction"] = max_saved_force_residual_difference < 1e-19
    evidence.update(
        force_samples_rechecked=count,
        max_recomputed_force_residual_N=max_recomputed_force_error,
        max_residual_difference_from_saved_N=max_saved_force_residual_difference,
        failures=violations,
    )
    gates = read(BASE / "analysis/FINAL_GATE.json")
    g = gates["groups"]
    checks["A_failure_does_not_block_B"] = g["STRONG"]["BASIC_NUMERICS_PASS"] == "FAIL" and all(
        g[k]["BASIC_NUMERICS_PASS"] == "PASS" for k in ["NAV_X", "NAV_Y"]
    )
    checks["actual_group_vs_ID_separation"] = (
        g["NAV_Y"]["GROUP_POSITION_NUMERICALLY_ADMITTED"] == "PASS"
        and g["NAV_Y"]["PARTICLE_ID_NUMERICALLY_ADMITTED"] == "FAIL"
        and g["NAV_Y"]["GROUP_SHAPE_NUMERICALLY_ADMITTED"] == "FAIL"
    )
    checks["targets_not_moved"] = (
        g["NAV_X"]["TARGET_REACHED"]["X"]
        and not g["NAV_Y"]["TARGET_REACHED"]["Y"]
        and p["XY_tolerance_m"] == 7.5e-8
    )
    costs = read(BASE / "analysis/COSTS.json")
    checks["budget_limits"] = costs["total_numerical_s"] < 600 and all(
        costs["routes_wall_s"][k] <= v for k, v in p["route_budget_s"].items()
    )
    checks["untrained_uncalibrated"] = (
        not gates["TRAINING_READY"] and gates["EXPERIMENTALLY_CALIBRATED"] == "NOT_VERIFIED"
    )
    for file in BASE.glob("execution/*_v1.json"):
        d = read(file)
        assert sha(BASE / "execution" / (d["sha"] + ".tar.gz")) == d["archive_sha256"]
        subprocess.run(["git", "cat-file", "-e", d["sha"] + "^{commit}"], cwd=ROOT, check=True)
    checks["immutable_remote_source_packages"] = True
    if include_visual:
        v = read(BASE / "figures/VISUAL_REVIEW.json")
        pdf = read(BASE / "figures/PDF_CHECKS.json")
        checks["all_8_actual_images_reviewed"] = (
            v["status"] == "PASS" and len(v["opened_and_reviewed"]) == 8 and pdf["status"] == "PASS"
        )
        for rel, digest in read(BASE / "figures/FIGURE_PROVENANCE.json")["figures_sha256"].items():
            assert sha(BASE / "figures" / rel) == digest
        pr = read(BASE / "progress/UPDATE_AUDIT.json")
        checks["progress_history_preserved"] = (
            pr["all_checks_pass"]
            and sha("/home/lzy/projects/microbubble_aggregation.md") == pr["after_sha256"]
        )
    result = dict(
        status="PASS" if all(checks.values()) else "FAIL",
        checks=checks,
        evidence=evidence,
        wall_s=time.perf_counter() - start,
        identity=identity(),
    )
    write(BASE / "VERIFICATION.json" if include_visual else BASE / "PRELIMINARY_VERIFICATION.json", result)
    assert all(checks.values()), checks
    print(result)
    return result


if __name__ == "__main__":
    import sys

    verify("--preliminary" not in sys.argv)
