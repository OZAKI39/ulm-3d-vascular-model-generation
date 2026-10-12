"""Read-only synthesis of archived science, including the rejected GPU candidate."""

import numpy as np
from mevo_mb.fast_sim.f10g_contact import resistance, normals
from .common import ROOT, BASE, js, write, identity, sha


def main():
    a = BASE / "analysis"
    p = js(ROOT / "configs/fast_pzt/f10i/protocol.json")
    raw = js(BASE / "gpu_v2/RESULTS.json")
    row = raw["rows"][-1]
    assert row["name"] == "saved32_O_1102" and not row["pass_"]
    r, kc, dt = p["radius_m"], p["kc_N_m"], row["dt_s"]
    x = np.array(row["initial_positions_m"])
    y = np.array(row["gpu"]["positions_m"])
    free = np.array(row["free_velocity_m_s"])
    A = resistance(x, r)["A"]
    zeta = 6 * np.pi * p["viscosity_Pa_s"] * r
    i, j, distance, normal = normals(y)
    pair_force = kc * np.maximum(2 * r - distance, 0)
    force = np.zeros_like(y)
    np.add.at(force, i, -pair_force[:, None] * normal)
    np.add.at(force, j, pair_force[:, None] * normal)
    pred = x + dt * np.linalg.solve(A, free.ravel()).reshape(-1, 3)
    residual = y - pred - dt / zeta * np.linalg.solve(A, force.ravel()).reshape(-1, 3)
    # Offline exact CPU resistance gives an independent force/position residual.
    gradient = A @ ((y - pred) / r).ravel() - dt / (zeta * r) * force.ravel()
    assert np.all(pair_force >= 0) and np.all(pair_force[distance >= 2 * r] == 0)
    assert np.linalg.eigvalsh(A)[0] > 0.999999
    assert np.linalg.norm(residual, axis=1).max() < p["gpu_comparison"]["equation_residual_m"]
    assert raw["status"] == "GPU_B2_NOT_ADMITTED"
    small = raw["rows"][:-1]
    assert len(small) == 18 and all(q["pass_"] for q in small)
    negative = [q for q in small if q["contact_quality"] == "EXPECTED_BUDGET_REJECTION"]
    assert len(negative) == 2
    result = dict(
        source=identity(),
        raw_sha256=sha(BASE / "gpu_v2/RESULTS.json"),
        engineering_v1=js(BASE / "gpu_v1/ENGINEERING_FAILURE.json"),
        small_reference_agreement_count=18,
        small_positive_quality_count=16,
        negative_budget_fixture_count=2,
        saved32_metrics=row["metrics"],
        saved32_solver_converged=False,
        saved32_stop_reason=row["gpu"]["reason"],
        stored_conservative_gradient_bound_m=row["gpu"]["conservative_equation_bound_m"],
        offline_reconstructed_gradient_bound_m=float(r * np.linalg.norm(gradient)),
        offline_equation_residual_m=float(np.linalg.norm(residual, axis=1).max()),
        resistance_eigenvalue_range=np.linalg.eigvalsh(A)[[0, -1]],
        nonnegative_pair_force=True,
        separated_pair_force_zero=True,
        accepted_outer_updates=len(row["gpu"]["iterations"]),
        base_PCG_tolerance_pass=False,
        final_metric_PCG_tolerance_pass=True,
        stored_linear_solves=row["gpu"]["linear_solves"],
        cpu_timing_s={k: row["cpu"][k] for k in ["initialization_s", "contact_s", "total_s"]},
        gpu_timing=row["gpu"]["timing"],
        construction=raw["constructions"][-1],
        benchmark="NOT_RUN_AFTER_SOLVER_GATE_FAIL",
        memory_scope="Explicit 78376 bytes excludes PCG private arrays, CUDA context and graph storage. Process samples cannot establish a high-water memory bound.",
        diagnosis="Reference errors and original overlap/correction budgets pass, but the registered conservative gradient bound did not reach 1e-14m; Armijo exhausted 20 trials. The frozen-resistance PCG recursive residual also exceeded its own target after 128 iterations. Objective cancellation or base-solve error are possible contributors, not isolated by this bounded experiment. No tolerance or physics was changed to relabel the candidate.",
        decision="GPU_B2_NOT_ADMITTED",
        production=False,
        training_ready=False,
    )
    write(a / "GPU_FAILURE_ANALYSIS.json", result)
    e = js(a / "EVENT_MATCHES.json")
    g = js(a / "GROUP_TOPOLOGY.json")
    s = js(a / "SOLVER_TOLERANCE.json")
    t = js(a / "TOLERANCE_PROPAGATION.json")
    checks = {
        "centroid_norm_m": "centroid_m",
        "max_fixed_ID_m": "fixed_ID_m",
        "Kabsch_RMSD_m": "Kabsch_m",
        "Rg_relative": "Rg_relative",
        "aspect_relative": "aspect_relative",
        "tensor_relative": "tensor_relative",
    }
    assert all(g["maxima"][k] <= p["historical_morphology_limits"][v] for k, v in checks.items())
    assert len(e["events"]) == 29 and e["timing_fail_count"] == 12
    assert e["edge_id_sequence_pass"] and not e["cross_pair_order_pass"]
    assert all(v["members_equal"] for v in g["comparison"])
    assert all(set(v["component_count"]) == {1} for v in g["modes"].values())
    counts = {
        str(tol): sum(
            v["result"]["converged"] for row in s["rows"] for v in row["variants"] if v["tolerance_m"] == tol
        )
        for tol in p["cpu_tolerances_m"]
    }
    assert list(counts.values()) == [14, 14, 13]
    assert t["trajectory"]["accepted_steps"] == 200 and t["trajectory"]["first_rejection"] is None
    assert t["same_dt_tolerance_comparison"]["matched_edge_max_time_difference_s"] == 0
    write(
        a / "FINAL_GATE.json",
        dict(
            source=identity(),
            raw_edge_timing="RAW_EDGE_TIMING_FAIL",
            raw_edge_failures=12,
            matched_events=29,
            max_event_difference_s=e["max_event_difference_s"],
            edge_identity="EDGE_ID_SEQUENCE_PASS",
            cross_pair_order="FAIL",
            group_position="GROUP_POSITION_PASS",
            group_shape="GROUP_SHAPE_PASS",
            group_member_state="IDENTICAL_SINGLE_COMPONENT_32_IDS_ALL_ACCEPTED_FRAMES",
            task_level_evidence="GROUP_LEVEL_NUMERICAL_EVIDENCE",
            persistent_observation="NO_PERSISTENT_SPLIT_OR_MERGE_OBSERVED",
            persistent_event_accuracy="PERSISTENT_TOPOLOGY_NOT_EVALUATED",
            right_censored_events=0,
            cpu_converged_counts=counts,
            cpu_states=14,
            tolerance_control="ONE_ORIGINAL_INITIAL_STATE_2MS_DT10_TOL14_COMPLETE",
            finer_2p5us="NOT_RUN",
            gpu="GPU_B2_NOT_ADMITTED",
            gpu_static_numerical_agreement="18_SMALL_CASES_MATCH_CPU; SAVED32_COMPARISON_PASS_BUT_SOLVER_NOT_CONVERGED",
            physical_status="EXPLORATORY / PHYSICS_UNCALIBRATED",
            inherited_full_morphology_gate="MORPHOLOGY_NOT_CONVERGED",
            production_default_changed=False,
            training_ready=False,
            new_10ms=False,
            new_seeds=False,
            B64=False,
            next="F10J: independently enforce true base/outer residuals on the same fixed state and reduce host synchronization before any production proposal; then seek bounded task-level persistent-topology validation under a separately registered protocol. Do not chase instantaneous 130us alone.",
        ),
    )
    print("Final layered gates written; GPU rejected, original strict event failure retained.")


if __name__ == "__main__":
    main()
