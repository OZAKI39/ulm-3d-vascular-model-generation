"""Separate quality, fixed-ID sensitivity, target success and model calibration."""

import numpy as np
import h5py
from .common import BASE, CONFIG, read, write, sha


def main():
    p = read(CONFIG)
    runs = read(BASE / "analysis/RUN_ANALYSIS.json")
    comp = read(BASE / "analysis/COMPARISONS.json")
    costs = read(BASE / "analysis/COSTS.json")
    strong = [runs["validation_v1/STRONG_" + s] for s in ["DT1", "DT05", "FINE"]]
    weak = runs["validation_v1/WEAK_DT1"]
    checks = read(BASE / "pilot_v1/CHECKS.json")
    failure = []
    for r in strong:
        path = BASE / "validation_v1/runs" / r["name"]
        row = read(path / "LAST_STEP.json")
        b = row["result"]
        overlap = np.array(b["geometric_overlap_m"])
        perforce = np.array(b["pair_force_N"])
        pairs = np.array(b["pairs"])
        j = int(np.argmax(overlap))
        failure.append(
            dict(
                name=r["name"],
                last_accepted_time_s=r["completed_s"],
                attempted_time_s=r["completed_s"] + r["dt_s"],
                residual_m=b["mobility_equation_residual_m"],
                mathematical_pass=b["converged"],
                overlap_m=float(overlap[j]),
                overlap_budget_m=0.01 * p["radius_m"],
                overlap_excess_m=float(overlap[j] - 0.01 * p["radius_m"]),
                pair_ids=pairs[j],
                max_pair_contact_force_N=float(perforce.max()),
                max_pair_force_at_budget_N=p["kc_N_m"] * 0.01 * p["radius_m"],
                max_iteration_m=b["max_single_iteration_m"],
                wall_correction_m=row["wall_correction_m"],
                reason=r["reasons"],
                data_sha256=sha(path / "LAST_STEP.json"),
            )
        )
    solver_s = resistance_s = step_s = 0
    for rel, r in runs.items():
        with h5py.File(BASE / rel.split("/")[0] / "runs" / rel.split("/")[1] / "trajectory.h5") as h:
            m = h["step_metrics"][:]
            step_s += m[:, 6].sum()
            solver_s += m[:, 7].sum()
            resistance_s += m[:, 8].sum()
    weak_dt = comp["WEAK_DT05"]
    weak_cache = comp["WEAK_FINE"]
    quality = all(r["status"] == "PASS" for r in strong)
    particle = quality and all(comp["STRONG_" + s]["fixed_id_pass"] for s in ["DT05", "FINE"])
    gate = dict(
        cleanup_status="PASS",
        cleanup_actual_deleted_code_lines=5,
        cleanup_helper="vessel_kernels.f64",
        cleanup_commit=p["cleanup_commit"],
        bridge_commit="3f5a117c1",
        engineering_status=checks["status"],
        original_CPU_B2_solver_reused=True,
        kc_N_m=p["kc_N_m"],
        strong_B2_mathematical_pass=all(f["mathematical_pass"] for f in failure),
        strong_B2_2ms_quality_pass=quality,
        strong_B2_2ms_fixed_id_pass=particle,
        strong_short_prefix_dt=comp["STRONG_DT05"],
        strong_short_prefix_cache=comp["STRONG_FINE"],
        strong_failures=failure,
        XY_targets_m=p["XY_targets_m"],
        XY_tolerance_m=p["XY_tolerance_m"],
        navigation_status="NOT_RUN_STRONG_B2_QUALITY_GATE_FAILED",
        X_target_reached=False,
        Y_target_reached=False,
        XY_independent_control_demonstrated=False,
        Z_independent_control_demonstrated=False,
        longest_completed_s=max(r["completed_s"] for r in runs.values()),
        longest_strong_valid_prefix_s=max(r["completed_s"] for r in strong),
        windows={
            "2ms": "OFF_AND_WEAK_COMPLETE; STRONG_REJECTED",
            "5ms": "NOT_RUN",
            "10ms": "NOT_RUN",
            "30ms": "NOT_RUN",
        },
        weak_recheck_status="PASS"
        if weak["status"] == "PASS" and weak_dt["fixed_id_pass"] and weak_cache["fixed_id_pass"]
        else "FAIL",
        weak_Rg_change_m=weak["rg_initial_final_m"][1] - weak["rg_initial_final_m"][0],
        weak_xy_change=weak["xy_initial_final"][1] - weak["xy_initial_final"][0],
        weak_target_reached=weak["xy_initial_final"][1] >= p["shape_target_xy"],
        weak_effect_exceeds_dt_cache=bool(
            abs(weak["xy_initial_final"][1] - 1)
            > max(weak_dt["directional_xy"], weak_cache["directional_xy"])
        ),
        numerical_process_total_s=costs["total_process_s"],
        integration_s=costs["integration_s"],
        measured_step_sum_s=step_s,
        CPU_solver_sum_s=solver_s,
        CPU_resistance_sum_s=resistance_s,
        GPU_kernel_only_time_s=None,
        GPU_time_note="Existing CUDA forces and walls used; host step timing includes transfers. No isolated GPU kernel timing claim.",
        peak_process_rss_MiB=max(r["rss_peak_MiB"] for r in runs.values()),
        scope="EXPLORATORY / PHYSICS_UNCALIBRATED / HOMOGENEOUS_WATER_CONTROL_PROXY",
        physical_calibration_pass=False,
        training_ready=False,
        stop_reason="Geometric soft overlap exceeds unchanged 0.01R before 0.175ms in all three strong tests. No navigation search or extension admitted. Weak-control shape analysis reuses required numerical comparison; no extra shape optimization.",
    )
    assert not quality and all("RESIDUAL_OVERLAP_BUDGET" in f["reason"] for f in failure)
    write(BASE / "analysis/FINAL_GATE.json", gate)
    print(gate)


if __name__ == "__main__":
    main()
