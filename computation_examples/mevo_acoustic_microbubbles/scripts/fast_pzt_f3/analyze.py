"""Offline F3 inference; distinguish exact repeats, initial sensitivity, and control effects."""

import numpy as np
import h5py
from mevo_mb.fast_sim.kinematics import rms_vectors
from scripts.fast_pzt_f2.analyze import time_windows
from .common import BASE, F2, PROTOCOL, write, js

REF = BASE / "repeats_v1/stable_reduction_0"
BRANCHES = dict(unchanged=REF, phase=BASE / "controls_v1/phase", amplitude=BASE / "controls_v1/amplitude")
KEYS = ["rg_relative", "eigenvalues_relative", "aligned_internal_rmsd_m"]


def main():
    driving = js(BASE / "driving_v3/RESULT.json")
    repeat = js(BASE / "repeats_v1/RESULT.json")
    controls = js(BASE / "controls_v1/RESULT.json")
    with h5py.File(BASE / "driving_v3/FORCE_FIELD_DATA.h5") as h:
        x = h["positions_m"][:]
        r = x - x.mean(0)
        zeta = 6 * np.pi * 0.001 * 1.78e-6
        s0 = h["unchanged/secondary_N"][:]
        for name, row in driving["candidate_rows"].items():
            s = h[name + "/secondary_N"][:]
            u = h[name + "/mobility_velocity_m_s"][:]
            f = h[name + "/total_N"][:]
            fd = h[name + "/primary_from_pressure_derivative_N"][:]
            alpha = float(np.sum(s * s0) / np.sum(s0 * s0))
            delta = s - s0
            residual = s - alpha * s0
            row["secondary_scale_decomposition"] = dict(
                scale_relative_to_baseline=alpha,
                nonuniform_residual_RMS_N=rms_vectors(residual),
                control_delta_RMS_N=rms_vectors(delta),
                fraction_of_delta_energy_in_global_secondary_scaling=float(
                    1 - np.sum(residual**2) / max(np.sum(delta**2), 1e-300)
                )
                if name != "unchanged"
                else None,
                interpretation="Least-squares decomposition on fixed positions; global scaling of the existing nonuniform Secondary pattern is distinct from direct differential Primary actuation",
            )
            torque = np.sum(np.cross(r, f), axis=0)
            drag_torque = zeta * np.sum(np.cross(r, u), axis=0)
            row["angular_balance_diagnostic"] = dict(
                acoustic_torque_N_m=torque,
                drag_velocity_torque_N_m=drag_torque,
                residual_N_m=drag_torque - torque,
                interpretation="Exact symmetric normal-resistance solve would balance total torque; truncated local iteration corrects mean force only. This audit does not change the frozen solver or prove a sole cause of rearrangement.",
            )
            row["field"]["differential_primary_from_pressure_derivative_RMS_N"] = rms_vectors(fd - fd.mean(0))
    write(BASE / "CONTROL_FORCE_DIAGNOSTICS.json", driving)
    write(
        BASE / "PZT_CONTROL_SCREENING.json",
        dict(
            selected=js(BASE / "driving_v3/SELECTED.json"),
            status=driving["screening_status"],
            differential_primary_status=driving["differential_driving_status"],
            candidates=driving["candidate_rows"],
            criteria=PROTOCOL["screening"],
            no_post_trajectory_selection=True,
        ),
    )
    repeat["fixed_state_checks"] = js(BASE / "checks_v2/RESULT.json")
    repeat["initial_condition_sensitivity"] = dict(
        controls["perturbations"],
        endpoint_comparisons={k: v for k, v in controls["analysis"].items() if k.startswith("initial_")},
    )
    repeat["interpretation"] = (
        "On tested hardware, pair atomic order is an experimentally observed source of initial force variation; both ordered reductions yield bitwise repeated trajectories. Component-only stabilization does not remove pair variation; component atomics alone did not vary in the fixed-state sample. Contact/mobility are fixed-order; graph/direct short replay agrees. Actual 1 pm/1 nm perturbations still amplify: repeatability is not robust individual-trajectory accuracy."
    )
    repeat["morphology_sensitivity_status"] = "MORPHOLOGY_NUMERICALLY_SENSITIVE"
    write(BASE / "NUMERICAL_REPRODUCIBILITY.json", repeat)
    noise = np.load(BASE / "repeats_v1/stable_reduction_envelope.npz")
    old_models = [np.load(F2 / "mature_v1" / n / "CONTRAST.npz") for n in ["F1_near_field", "mobility4"]]
    perturbs = [np.load(BASE / "controls_v1" / f"initial_perturb_{i}" / "CONTRAST.npz") for i in range(2)]
    thresholds = dict(
        rg_relative=PROTOCOL["observation"]["rg_relative"],
        eigenvalues_relative=PROTOCOL["observation"]["eigenvalues_relative"],
        aligned_internal_rmsd_m=PROTOCOL["observation"]["aligned_rmsd_over_radius"] * 1.78e-6,
    )
    envelopes = {
        k: np.maximum.reduce([a[k] for a in old_models + perturbs]) for k in KEYS + ["centroid_distance_m"]
    }
    np.savez_compressed(
        BASE / "SHAPE_UNCERTAINTY.npz",
        time_s=noise["time_s"],
        **{k: envelopes[k] for k in KEYS},
        **{"repeat_" + k: noise[k] for k in KEYS},
    )
    effects = {}
    for name in ["phase", "amplitude"]:
        c = np.load(BRANCHES[name] / "CONTRAST.npz")
        t = c["time_s"]
        metrics = {}
        robust_last = False
        transient = False
        for k, limit in thresholds.items():
            detected = (c[k] > limit) & (c[k] > 3 * noise[k])
            robust = detected & (c[k] > envelopes[k])
            windows = time_windows(t, robust)
            ends = [w for w in windows["intervals"] if abs(w["end_s"] - 0.12) < 1e-10]
            sustained_final = bool(
                ends and ends[0]["observed_span_s"] >= PROTOCOL["observation"]["sustained_window_s"] - 1e-14
            )
            robust_last |= sustained_final
            transient |= bool(windows["longest_sampled_span_s"] > 0)
            metrics[k] = dict(
                endpoint=float(c[k][-1]),
                maximum=float(c[k].max()),
                repeat_endpoint=float(noise[k][-1]),
                model_and_initial_envelope_endpoint=float(envelopes[k][-1]),
                above_repeat_windows=time_windows(t, detected),
                above_repeat_model_initial_windows=windows,
                endpoint_exceeds_all=bool(robust[-1]),
                sustained_through_endpoint=sustained_final,
                short_window_maximum=float(
                    c[k][t <= PROTOCOL["observation"]["short_window_s"] + 1e-14].max()
                ),
            )
        effects[name] = dict(
            status="ROBUST_SHAPE_RESPONSE"
            if robust_last
            else "TRANSIENT_SHAPE_RESPONSE"
            if transient
            else "NO_ROBUST_SHAPE_RESPONSE",
            endpoint_status="ROBUST_ENDPOINT_SHAPE_RESPONSE"
            if robust_last
            else "NO_ROBUST_ENDPOINT_SHAPE_RESPONSE",
            metrics=metrics,
            translation_endpoint_m=float(c["centroid_distance_m"][-1]),
            translation_resolved=bool(
                c["centroid_distance_m"][-1]
                > 3 * max(noise["centroid_distance_m"][-1], envelopes["centroid_distance_m"][-1])
            ),
            endpoint_details=controls["analysis"][name],
            scope="Fixed thresholds; same-mode max over three repeat pairs; descriptive historical F2 model and two initial-perturbation envelopes, not statistical confidence or physical calibration",
        )
    quality = {}
    for path in list(BRANCHES.values()) + [BASE / "controls_v1" / f"initial_perturb_{i}" for i in range(2)]:
        s = js(path / "SUMMARY.json")
        quality[str(path.relative_to(BASE))] = dict(
            status=s["numerical_quality"],
            max_accepted_overlap_nm=s["quality"]["max_accepted_penetration_m"] * 1e9,
            max_pass_correction_nm=s["quality"]["max_contact_pass_correction_m"] * 1e9,
            correction_over_mobility_path=s["quality"]["correction_fraction"],
            correction_over_free_path=s["quality"]["correction_over_free_Stokes_path"],
            wall_s=s["wall_s"],
            actual_steps=s["actual_substep_count"],
        )
    result = dict(
        effects=effects,
        physics_scope="EXPLORATORY",
        maturity_source="F2 102.2 ms, exact saved 32 IDs",
        controls=controls["selected_controls"],
        fixed_criteria=PROTOCOL["observation"],
        quality=quality,
        unchanged=controls["analysis"]["unchanged"],
        initial_sensitivity=controls["perturbations"],
        envelope_scope="F2 historical near-field/mobility variations at the same initial state plus separate F3 1 pm/1 nm initial perturbations. No F3 physical model change. This is a conservative descriptive comparison, not a validated model-error bound.",
        topology_status="NO_PERSISTENT_SPLIT"
        if all(
            controls["analysis"][n]["topology"]["counts"]["split_persistent"] == 0
            for n in ["unchanged", "phase", "amplitude"]
        )
        else "PERSISTENT_GRAPH_SPLIT_OBSERVED",
        inference_limit="Field gradients are separately interpolated nodal derivatives and disagree with local derivatives of interpolated pressure; no spatial-grid calibration. Kabsch angles with large residual deformation do not establish controlled rigid rotation.",
    )
    write(BASE / "MORPHODYNAMICS_COMPARISON.json", result)
    perfpath = BASE / "performance_v1/RESULT.json"
    if perfpath.exists():
        write(BASE / "PERFORMANCE.json", js(perfpath))
    print(
        {
            name: dict(
                status=e["status"],
                endpoint=e["endpoint_status"],
                translation_um=e["translation_endpoint_m"] * 1e6,
                windows={
                    k: v["above_repeat_model_initial_windows"]["longest_sampled_span_s"]
                    for k, v in e["metrics"].items()
                },
            )
            for name, e in effects.items()
        }
    )


if __name__ == "__main__":
    main()
