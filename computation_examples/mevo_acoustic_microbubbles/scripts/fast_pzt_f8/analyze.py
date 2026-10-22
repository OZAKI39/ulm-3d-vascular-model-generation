"""Closed accepted-data analysis; numerical, task and physical labels stay separate."""

import json
import h5py
import numpy as np
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import ROOT, BASE, js, sha, write


def load(path):
    s = js(path / "SUMMARY.json")
    cfg = js(path / "CONFIG.json")
    with h5py.File(path / "trajectory.h5") as h:
        d = dict(
            t=h["time_s"][:],
            x=h["positions_m"][:, 0],
            c=h["metrics/global_centroid_m"][:, 0],
            rg=h["metrics/global_rg_m"][:, 0],
            fraction=h["metrics/largest_size"][:, 0] / 32,
            quality=h["quality/raw_F1_diagonal"][:, 0],
            wall=h["wall_stats"][:, 0],
            active=h["active"][:, 0],
            error=h["trajectory/target_error_m"][:],
            flow=h["flow/velocity_m_s"][:, 0],
            control_index=h["control_index"][:],
            labels=h["metrics/labels"][:, 0],
        )
        assert np.isfinite(d["x"]).all()
        np.testing.assert_allclose(d["c"], d["x"].mean(1), atol=1e-17, rtol=0)
        target = np.asarray(cfg["environment"]["target_m"])
        expected = np.linalg.norm(
            (d["c"] - target)[:, 1:] if cfg["task_type"] == "flow_guided_transverse" else d["c"] - target,
            axis=1,
        )
        np.testing.assert_allclose(d["error"], expected, atol=1e-17, rtol=1e-10)
        assert json.loads(h["control"].attrs["events_json"]) == s["controls"]
        if s["status"] == "COMPLETE":
            np.testing.assert_allclose(
                h["trajectory/physical_dt_s"][:].sum(), cfg["duration_s"], atol=1e-12, rtol=0
            )
        assert list(h["particle_id"][:]) == list(range(32))
        assert all(
            k in h
            for k in ["source", "medium", "flow", "bubbles", "control", "trajectory", "morphology", "quality"]
        )
    for k, v in dict(
        control_t=[e["time_s"] for e in s["controls"]],
        amplitudes=[e["amplitude_Pa"] for e in s["controls"]],
        phases=[e["phase_rad"] for e in s["controls"]],
    ).items():
        d[k] = np.array(v)
    gate = js(path / "GATE_CROSSING.json") if (path / "GATE_CROSSING.json").exists() else None
    if gate and gate["status"] == "CROSSED":
        lo, hi = np.asarray(gate["accepted_before_m"]), np.asarray(gate["accepted_after_m"])
        cross = lo + gate["fraction"] * (hi - lo)
        np.testing.assert_allclose(cross, gate["derived_positions_m"], rtol=0, atol=1e-18)
        assert abs(cross.mean(0)[0] - cfg["environment"]["x_gate_m"]) < 1e-14
        assert gate["bracket_width_s"] <= cfg["dt_s"] * 1.000001
        d["cross_x"] = cross
        d["cross_c"] = cross.mean(0)
        d["cross_t"] = np.array(gate["time_s"])
    return s, cfg, d, gate


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    proto = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    paths = {k: BASE / "a_formal_v1" / k for k in ["A0", "A1"]}
    paths.update({k: BASE / "y_formal_v1" / k for k in ["Y0", "Y1", "Y2"]})
    paths.update(
        {
            k + "_half": BASE / "y_half_v1" / k
            for k in ["Y0", "Y1", "Y2"]
            if (BASE / "y_half_v1" / k / "SUMMARY.json").exists()
        }
    )
    loaded = {k: load(p) for k, p in paths.items()}
    physics = F2Physics()
    response = ResponseSection(ROOT / "results/fast_pzt/f5/inputs/f4/response.h5", physics)
    arrays = {}
    rows = {}
    provenance = {}
    for name, (s, cfg, d, gate) in loaded.items():
        for key, v in d.items():
            arrays[name + "_" + key] = v
        opts = s.get("feedback_optimization", [])
        rows[name] = {
            k: s[k]
            for k in [
                "status",
                "reason",
                "numerical_quality",
                "navigation",
                "wall_s",
                "process_body_wall_s",
                "quality",
            ]
        }
        rows[name].update(
            task_type=cfg["task_type"],
            production_mobility=cfg["mobility"],
            prediction_mobility=cfg["prediction_mobility"],
            gate=gate,
            min_largest_fraction=d["fraction"].min(),
            final_largest_fraction=d["fraction"][-1],
            all_IDs_active=bool(d["active"].all()),
            control_updates=len(opts),
            optimizer_total_s=sum(e["wall_s"] for e in opts),
            optimizer_calls=sum(e["evaluations"] for e in opts),
            optimizer_max_iterations=max([e["optimizer_iterations"] for e in opts] or [0]),
        )
        if name.startswith("Y"):
            with h5py.File(paths[name] / "trajectory.h5") as h:
                pressure = h["control/pressure_Pa"][:]
                gradient = h["control/gradient_Pa_m"][:]
            volume = response.query(pressure.ravel())[:, 0].reshape(pressure.shape)
            primary = -0.5 * np.real(volume[..., None] * gradient.conj())
            primary_velocity = primary.mean(1) / (6 * np.pi * physics.viscosity_Pa_s * physics.radius_m)
            meanflow = d["flow"].mean(1)
            arrays[name + "_sampled_primary_mean_m_s"] = primary_velocity
            arrays[name + "_sampled_CFD_mean_m_s"] = meanflow
            span = np.diff(d["t"])
            flow_displacement = (meanflow[:-1] * span[:, None]).sum(0)
            primary_displacement = (primary_velocity[:-1] * span[:, None]).sum(0)
            arrays[name + "_sampled_wall_residual_m"] = (
                d["c"]
                - d["c"][0]
                - np.vstack(
                    [
                        np.zeros((1, 3)),
                        np.cumsum((meanflow[:-1] + primary_velocity[:-1]) * span[:, None], axis=0),
                    ]
                )
            )
            rows[name]["sampled_drift_accounting"] = dict(
                CFD_left_sum_m=flow_displacement,
                primary_Stokes_left_sum_m=primary_displacement,
                observed_displacement_m=d["c"][-1] - d["c"][0],
                residual_m=d["c"][-1] - d["c"][0] - flow_displacement - primary_displacement,
                total_wall_correction_norm_per_ID_m=s["quality"]["wall"]["cumulative_correction_m"] / 32,
                scope="0.5ms sampled left sums, not per-step force integrals. M0 component balance and pairwise central forces preserve the mean free velocity; geometric wall corrections and quadrature remain. No measured-physics decomposition is claimed.",
            )
        assert rows[name]["optimizer_max_iterations"] <= 20
        assert d["amplitudes"].max() <= proto["amplitude_cap_Pa"] + 1e-10
        provenance[name] = dict(
            path=str(paths[name].relative_to(ROOT)),
            execution_commit=js(paths[name] / "IDENTITY.json")["execution_commit"],
            files={p.name: sha(p) for p in paths[name].iterdir() if p.is_file()},
        )
    # The same F7 RPY controller is reproduced despite enabling the common repair policy.
    with h5py.File(ROOT / "results/fast_pzt/f7/controller_v2/run/trajectory.h5") as h:
        exact_F7 = bool(np.array_equal(h["positions_m"][:, 0], loaded["A1"][2]["x"]))
    assert exact_F7
    p0 = loaded["A0"][0]["feedback_policy"].copy()
    p1 = loaded["A1"][0]["feedback_policy"].copy()
    assert p0.pop("prediction_mobility") == "diagonal_stokes"
    assert p1.pop("prediction_mobility") == "rpy_free_space" and p0 == p1
    for key in ["A1", "Y0", "Y1", "Y2"]:
        reference = "A0" if key == "A1" else "Y0"
        np.testing.assert_array_equal(loaded[key][2]["x"][0], loaded[reference][2]["x"][0])
    convergence = {}
    for name in ["Y0", "Y1", "Y2"]:
        if name + "_half" not in loaded:
            continue
        s, cfg, d, g = loaded[name]
        sh, _, dh, gh = loaded[name + "_half"]
        np.testing.assert_allclose(d["t"], dh["t"], rtol=0, atol=1e-13)
        delta = np.linalg.norm(d["c"] - dh["c"], axis=1)
        error = np.abs(d["error"] - dh["error"])
        rg = np.max(abs(d["rg"] - dh["rg"]) / np.maximum(d["rg"], 1e-30))
        frac = np.max(abs(d["fraction"] - dh["fraction"]))
        cross_error = (
            np.linalg.norm(d["cross_c"][1:] - dh["cross_c"][1:])
            if g["status"] == gh["status"] == "CROSSED"
            else None
        )
        convergence[name] = dict(
            max_centroid_m=delta.max(),
            max_transverse_error_m=error.max(),
            max_Rg_relative=rg,
            max_fraction_difference=frac,
            gate_transverse_difference_m=cross_error,
            gate_time_difference_s=abs(g["time_s"] - gh["time_s"]) if cross_error is not None else None,
            passed=bool(
                s["numerical_quality"] == sh["numerical_quality"] == "PASS"
                and delta.max() <= proto["timestep"]["max_centroid_m"]
                and rg <= proto["timestep"]["max_Rg_relative"]
                and frac <= proto["timestep"]["max_fraction_difference"]
            ),
        )
    design = proto["Y_design"]
    pair = {}
    d2 = loaded["Y2"][2]
    g2 = loaded["Y2"][3]
    for control in ["Y0", "Y1"]:
        dc = loaded[control][2]
        gc = loaded[control][3]
        np.testing.assert_allclose(dc["t"], d2["t"], rtol=0, atol=1e-13)
        inc = d2["c"] - dc["c"]
        arrays["Y2_minus_" + control] = inc
        gate_delta = d2["cross_c"] - dc["cross_c"] if g2["status"] == gc["status"] == "CROSSED" else None
        uncertainty = None
        if all(k in convergence for k in ["Y2", control]):
            uncertainty = (
                convergence["Y2"]["gate_transverse_difference_m"]
                + convergence[control]["gate_transverse_difference_m"]
            )
        effect_ok = bool(
            gate_delta is not None
            and uncertainty is not None
            and gate_delta[2]
            > max(design["minimum_target_increment_m"], design["numerical_error_multiplier"] * uncertainty)
        )
        pair[control] = dict(
            common_time_final_delta_m=inc[-1],
            common_time_final_error_reduction_m=dc["error"][-1] - d2["error"][-1],
            gate_delta_m=gate_delta,
            gate_target_direction_increment_m=None if gate_delta is None else gate_delta[2],
            gate_error_reduction_m=None
            if gate_delta is None
            else gc["transverse_error_m"] - g2["transverse_error_m"],
            gate_pair_discretization_bound_m=uncertainty,
            resolved_target_increment=effect_ok,
        )
    numerical = all(loaded[k][0]["numerical_quality"] == "PASS" for k in ["Y0", "Y1", "Y2"])
    nav = bool(
        numerical
        and g2["status"] == "CROSSED"
        and g2["transverse_error_m"] <= design["target_tolerance_m"]
        and g2["largest_fraction"] >= 0.9
        and all(q["resolved_target_increment"] and q["gate_error_reduction_m"] > 0 for q in pair.values())
        and all(convergence.get(k, {}).get("passed", False) for k in ["Y0", "Y1", "Y2"])
    )
    cohesion = bool(
        g2["status"] == "CROSSED"
        and g2["largest_fraction"] >= 0.9
        and loaded["Y0"][3].get("largest_fraction", 1) < 0.9
    )
    # COHESION_ONLY requires no resolved target improvement; otherwise distinguish partial steering from reaching.
    resolved = all(v["resolved_target_increment"] for v in pair.values())
    label = (
        "TARGET_NAVIGATION_PASS"
        if nav
        else "COHESION_ONLY"
        if cohesion and not resolved
        else "NAVIGATION_FAIL"
    )
    contact = js(BASE / "contact_v1/CHECKS.json") if (BASE / "contact_v1/CHECKS.json").exists() else None
    # A campaign conclusion must not label the sound-off control as acoustically cohesive.
    for name, row in rows.items():
        if name.startswith("A"):
            classification = (
                "TARGET_NAVIGATION_PASS" if row["navigation"]["arrived_and_held"] else "NAVIGATION_FAIL"
            )
            role = "single_factor_feedback"
        elif name == "Y0":
            classification, role = "NAVIGATION_FAIL", "CFD_only_control"
        elif name.startswith("Y1"):
            classification, role = "COHESION_ONLY", "fixed_PZT_control"
        else:
            classification, role = label, "target_feedback"
        row["final_evaluation"] = dict(
            classification=classification,
            role=role,
            scope="Post-run evaluation; frozen per-run SUMMARY remains the original execution record.",
        )
    final = dict(
        software_checks="PASS",
        A_only_factor_check="PASS",
        A1_bitwise_reproduces_F7=exact_F7,
        A0_navigation_pass=loaded["A0"][0]["navigation"]["arrived_and_held"],
        A1_navigation_pass=loaded["A1"][0]["navigation"]["arrived_and_held"],
        Y_numerical="NUMERICAL_PASS" if numerical else "NUMERICAL_FAIL",
        Y_navigation=label,
        Y_resolved_additional_target_motion=resolved,
        Y_connectivity_preserved=cohesion,
        Y_key_active_vs_fixed_timestep_checks=all(
            convergence.get(k, {}).get("passed", False) for k in ["Y1", "Y2"]
        ),
        Y0_half_step_status="NOT_RUN; descriptive sound-on/off difference lacks a matched Y0 discretization bound; no success inferred",
        Y_task_scale_resolution="NOT_RESOLVED: fixed-PZT gate-position step sensitivity exceeds target tolerance; broad group-level timestep PASS is insufficient for this small-effect task",
        contact="PENDING" if contact is None else contact["classification"],
        physics="PHYSICS_UNCALIBRATED",
        evidence_level="EXPLORATORY",
        figures="PENDING",
        archive="PENDING",
        progress_document="PENDING",
    )
    write(
        out / "RESULTS.json",
        dict(
            cases=rows,
            paired_Y=pair,
            timestep=convergence,
            provenance=provenance,
            partial_control_scope="Resolved +Z increments and target-region success are separate. No branching trajectory or fixed3D hold is claimed.",
        ),
    )
    np.savez_compressed(out / "plot_data.npz", **arrays)
    write(BASE / "FINAL_GATE.json", final)
    write(
        out / "AI_DATA_MANIFEST.json",
        dict(
            schema="F8_prescribed_downstream_gate",
            cases={
                k: dict(
                    **v,
                    task_type=rows[k]["task_type"],
                    numerical_label=rows[k]["numerical_quality"],
                    navigation_label=rows[k]["final_evaluation"]["classification"],
                    evaluation_role=rows[k]["final_evaluation"]["role"],
                    paired_family="original_Y_F5_300ms" if k.startswith("Y") else "F2_mature_no_flow",
                    training_ready=False,
                )
                for k, v in provenance.items()
            },
            caveats="Only two inherited initial-state families; no independent validation family. Never split frames or paired controls across train/test. All trajectories exploratory; failure/censored attempts excluded. End-run quality and full future controls are audit/labels, never causal features.",
        ),
    )
    print(json.dumps(final, indent=2))


if __name__ == "__main__":
    main()
