"""Small registered paired navigation campaign; no trajectory-driven retargeting."""

import json
import time
import h5py
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f9_ledger import (
    AuditSimulator as GateSimulator,
    AuditFeedbackSimulator as GateFeedbackSimulator,
)
from mevo_mb.fast_sim.f9_ledger import TERMS, EVENT_COLUMNS
from mevo_mb.fast_sim.metrics import frame_metrics
from mevo_mb.fast_sim.flow_gpu import AffineFlow
from mevo_mb.fast_sim.storage import save_run
from mevo_mb.fast_sim.morphodynamics import morphology_series
from mevo_mb.fast_sim.consistent_field import query_numpy
from mevo_mb.fast_sim.ai_samples import read_group_sample
from scripts.fast_pzt_f2.common import quality_pass
from .common import clean, write, identity


def arrival(times, errors, largest, proto):
    hit = (errors <= proto["target_radius_m"]) & (largest >= proto["min_largest_fraction"])
    first = None
    longest = 0.0
    start = None
    for i, value in enumerate(hit):
        if value:
            if first is None:
                first = float(times[i])
            if start is None:
                start = i
            longest = max(longest, float(times[i] - times[start]))
        else:
            start = None
    return dict(
        first_arrival_s=first,
        longest_sampled_dwell_s=longest,
        arrived_and_held=longest >= proto["dwell_s"] - 1e-12,
    )


def run_case(
    out,
    env,
    event,
    proto,
    objects_,
    *,
    dt=1e-5,
    mobility="f5_local",
    gradient="consistent_trilinear",
    duration=None,
    adaptive=False,
    feedback=False,
    audit=True,
):
    out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    p, basis, response, geometry, realflow = objects_
    duration = env["duration_s"] if duration is None else duration
    flow = (
        realflow
        if env["name"] == "Y_CFD"
        else AffineFlow(env["flow_velocity_m_s"], env["flow_gradient_s_inv"], env["initial_centroid_m"])
        if env["flow_enabled"]
        else None
    )
    if flow is not None and env["name"] != "Y_CFD":
        flow.metadata["scope"] = "SIMPLIFIED_SLOW_AFFINE_FLOW; not Y CFD"
    events = event if isinstance(event, list) else [event]
    factory = GateSimulator
    extra = dict(x_gate=env.get("x_gate_m"), audit_enabled=audit)
    if feedback:
        from mevo_mb.fast_sim.f8_navigation import MatchedFeedback

        policy = dict(
            proto["feedback_policy"],
            flow_enabled=env["flow_enabled"],
            prediction_mobility="diagonal_stokes" if feedback == "stokes" else "rpy_free_space",
            task_type=env.get("task_type", "fixed_3D_target"),
        )
        center = np.array(env["initial_centroid_m"])
        u0 = np.array(env["flow_velocity_m_s"])
        tensor = np.array(env["flow_gradient_s_inv"])
        flowfn = (
            (lambda x: realflow.sample(x)["velocity_m_s"])
            if env["name"] == "Y_CFD"
            else (lambda x: u0 + (x - center) @ tensor.T)
            if env["flow_enabled"]
            else (lambda x: np.zeros_like(x))
        )
        BoundedFeedback = MatchedFeedback
        FeedbackSimulator = GateFeedbackSimulator
        extra["controller"] = BoundedFeedback(basis, response, p, env["target_m"], flowfn, policy)
        factory = FeedbackSimulator
    sim = factory(
        np.array(env["positions_m"]),
        p,
        response,
        basis,
        **extra,
        flow=flow,
        flow_enabled=env["flow_enabled"],
        acoustic_enabled=any(events[0]["enabled"]),
        vessel=geometry,
        vessel_enabled=env["wall_enabled"],
        gradient_mode=gradient,
        mobility_mode=mobility,
    )
    config = dict(
        audit_enabled=audit,
        environment=env,
        events=events,
        dt_s=dt,
        duration_s=duration,
        adaptive=adaptive,
        mobility=mobility,
        gradient=gradient,
        feedback=feedback,
        task_type=env.get("task_type", "fixed_3D_target"),
        prediction_mobility=("diagonal_stokes" if feedback == "stokes" else "rpy_free_space")
        if feedback
        else None,
        contact_model="central_PBD",
    )
    ident = dict(identity(), basis=basis.metadata, response=response.metadata, geometry=geometry.metadata)
    write(out / "CONFIG.json", config)
    write(out / "IDENTITY.json", ident)
    gpu_start = wp.Event(enable_timing=True)
    gpu_end = wp.Event(enable_timing=True)
    wp.record_event(gpu_start)
    try:
        result = sim.run(
            duration, dt, events, proto["sample_interval_s"], wall_budget_s=360, adaptive=adaptive
        )
    except Exception as exc:
        write(
            out / "UNHANDLED_FAILURE.json",
            dict(
                reason=repr(exc),
                positions_m=sim.x.numpy(),
                quality=sim.quality(),
                feedback_events=getattr(sim, "feedback_events", []),
                failed_update=getattr(getattr(sim, "controller", None), "failed_update", None),
            ),
        )
        raise
    wp.record_event(gpu_end)
    wp.synchronize_event(gpu_end)
    gpu_elapsed = wp.get_event_elapsed_time(gpu_start, gpu_end) / 1000
    output_started = time.perf_counter()
    metrics, _ = save_run(out / "trajectory.h5", result, p, ident, config)
    if "failure_proposal" in result:
        write(out / "FAILURE_PROPOSAL.json", result["failure_proposal"])
    s = result["summary"]
    q = s["quality"]
    all_active = bool(result.get("active", np.ones(1)).all())
    passed = s["status"] == "COMPLETE" and quality_pass(q, p) and all_active
    if env["wall_enabled"]:
        passed = (
            passed
            and q["wall"]["minimum_queried_gap_m"] >= -0.01 * p.radius_m
            and q["wall"]["max_pass_correction_m"] <= 0.1 * p.radius_m
        )
    centroid = metrics["global_centroid_m"][:, 0]
    errors = np.linalg.norm(centroid - env["target_m"], axis=-1)
    fraction = metrics["largest_size"][:, 0] / sim.n
    nav = arrival(result["time_s"], errors, fraction, proto)
    nav.update(
        initial_error_m=float(errors[0]),
        final_error_m=float(errors[-1]),
        minimum_error_m=float(errors.min()),
        centroid_displacement_m=centroid[-1] - centroid[0],
        target_progress_m=float(
            (centroid[-1] - centroid[0]) @ (np.array(env["target_m"]) - centroid[0]) / errors[0]
        ),
        final_largest_fraction=float(fraction[-1]),
        final_straggler_count=int(sim.n - metrics["largest_size"][-1, 0]),
        classification="ARRIVAL_HOLD_CANDIDATE"
        if passed and nav["arrived_and_held"]
        else "NUMERICAL_FAILURE"
        if not passed
        else "NOT_REACHED_OR_NOT_HELD",
    )
    gate = sim.gate_record()
    if gate is not None:
        if gate["status"] == "CROSSED":
            gx = np.asarray(gate["derived_positions_m"])
            gm = frame_metrics(gx, p.radius_m)
            target_yz = np.asarray(env["target_m"])[1:]
            gate.update(
                transverse_error_m=float(np.linalg.norm(gx.mean(0)[1:] - target_yz)),
                bubble_fraction_in_target=float(
                    np.mean(np.linalg.norm(gx[:, 1:] - target_yz, axis=1) <= env["target_tolerance_m"])
                ),
                largest_fraction=gm["largest_size"] / sim.n,
                endpoint_largest_fractions=[
                    frame_metrics(np.asarray(gate[k]), p.radius_m)["largest_size"] / sim.n
                    for k in ["accepted_before_m", "accepted_after_m"]
                ],
            )
        write(out / "GATE_CROSSING.json", gate)
        s["gate_crossing"] = {
            k: v
            for k, v in gate.items()
            if k not in ["accepted_before_m", "accepted_after_m", "derived_positions_m"]
        }
        errors = np.linalg.norm(centroid[:, 1:] - np.asarray(env["target_m"])[1:], axis=1)
        nav = dict(
            task_type="flow_guided_transverse",
            final_transverse_error_m=float(errors[-1]),
            centroid_displacement_m=centroid[-1] - centroid[0],
            classification="PENDING_PAIRED_AND_TIMESTEP_EVALUATION" if passed else "NUMERICAL_FAILURE",
            fixed_32_ID_metrics_valid=all_active,
        )
    s.update(
        numerical_quality="PASS" if passed else "FAIL",
        navigation=nav,
        gpu_stream_elapsed_s=gpu_elapsed,
        gpu_time_scope="CUDA events around run including graph capture, synchronization and host-idle gaps; not exclusive active-kernel time",
        evidence_level="EXPLORATORY",
        fixed_ID_metrics_valid=all_active,
    )
    x = result["positions_m"][:, 0]
    morph = morphology_series(x, np.arange(sim.n), metrics["labels"][:, 0], 0.1) if all_active else {}
    # Saved-frame pressures come from original coordinates and actual controls.
    pressure = np.zeros(x.shape[:2], complex)
    grad = np.zeros(x.shape, complex)
    for ci, control in enumerate(s["controls"]):
        take = result["control_index"] == ci
        if not control["acoustic_enabled"]:
            continue
        field, _, _ = basis.combine(control["amplitude_Pa"], control["phase_rad"], control["enabled"])
        if gradient == "consistent_trilinear":
            pp, gg = query_numpy(field.pressure_Pa, basis.origin, basis.spacing, x[take].reshape(-1, 3))
        else:
            pp = field.interpolate_complex_pressure(x[take].reshape(-1, 3))
            gg = field.interpolate_complex_gradient(x[take].reshape(-1, 3))
        pressure[take] = pp.reshape(-1, sim.n)
        grad[take] = gg.reshape(-1, sim.n, 3)
    with h5py.File(out / "trajectory.h5", "r+") as h:
        h.attrs.update(
            schema="fast_pzt_f10a_v1",
            evidence_level="EXPLORATORY",
            split_group=env["split_group"],
            partition=env["partition"],
            numerical_quality=s["numerical_quality"],
        )
        h.attrs["task_type"] = config["task_type"]
        h.attrs["production_mobility"] = mobility
        h.attrs["prediction_mobility"] = config["prediction_mobility"] or "none"
        h.attrs["contact_model"] = "central_PBD"
        h.require_group("source").attrs["identity_json"] = json.dumps(clean(ident))
        h.require_group("medium").attrs["json"] = json.dumps(
            dict(
                density_kg_m3=1000, viscosity_Pa_s=0.001, radial_effective_viscosity_Pa_s=0.002, model="water"
            )
        )
        h.require_group("control").attrs["events_json"] = json.dumps(clean(s["controls"]))
        h["control/pressure_Pa"] = pressure
        h["control/gradient_Pa_m"] = grad
        for key in ["particle_id", "radii_m"]:
            h["bubbles/" + key] = h[key]
        for key in ["time_s", "positions_m", "control_index"]:
            h["trajectory/" + key] = h[key]
        h["trajectory/physical_dt_s"] = result["physical_dt_s"]
        if "adaptive_attempts" in result:
            ds = h.create_dataset("quality/adaptive_attempts", data=result["adaptive_attempts"])
            ds.attrs["columns"] = (
                "start_s,trial_dt_s,error_m,budget_m,accepted,coarse_error_bits,fine_error_bits"
            )
        for key in ["active", "wall_gap_m", "wall_stats"]:
            values = result.get(key)
            if values is None:
                values = (
                    np.ones(x.shape[:2], int)[:, None]
                    if key == "active"
                    else np.full(x.shape[:2], np.nan)[:, None]
                    if key == "wall_gap_m"
                    else np.zeros((*x.shape[:2], 12))[:, None]
                )
            h.create_dataset(key, data=values, compression="gzip")
        for key, values in result["flow"].items():
            h.create_dataset("flow/" + key, data=values, compression="gzip")
        for key, values in morph.items():
            h.create_dataset("morphology/" + key, data=values, compression="gzip")
        h["morphology/metrics"] = h["metrics"]
        h["quality/raw_F1_diagonal"] = result["quality_samples"]
        h["trajectory/boundary_events_json"] = json.dumps(result.get("boundary_events", []))
        h["trajectory/target_error_m"] = errors
        if gate is not None:
            h["trajectory"].attrs["gate_crossing_json"] = json.dumps(clean(gate))
        if "motion_ledger" in result:
            group = h.require_group("motion_ledger")
            group.attrs["terms"] = json.dumps(TERMS)
            group.attrs["step_meta_columns"] = (
                "start_s,dt_s,accepted,max_particle_residual_m,mean_pair_path_m,mean_wall_path_m,max_wall_metadata_residual_m"
            )
            group.attrs["event_columns"] = json.dumps(EVENT_COLUMNS)
            group.attrs["units"] = "SI; displacement vectors and positions m, time s; normals dimensionless"
            group.attrs["wall_face_index"] = (
                "index into original prepared_geometry faces[labels==1]; no remeshing"
            )
            for key, value in result["motion_ledger"].items():
                if key == "rejected_intervals":
                    for j, record in enumerate(value):
                        for subkey, arr in record.items():
                            group.create_dataset(f"rejected/{j}/{subkey}", data=arr, compression="gzip")
                else:
                    group.create_dataset(key, data=value, compression="gzip")
            from .validation import validate_ledger

            s["ledger_validation"] = validate_ledger(
                dict(result["motion_ledger"], rejected_intervals=[]), result["positions_m"], result["time_s"]
            )
        s["rejected_audit_intervals"] = len(result.get("motion_ledger", {}).get("rejected_intervals", []))
        s["output_and_postprocessing_s"] = time.perf_counter() - output_started
        s["process_body_wall_s"] = time.perf_counter() - started
        h.require_group("quality").attrs["summary_json"] = json.dumps(clean(s))
    s["process_body_wall_s"] = time.perf_counter() - started
    write(out / "SUMMARY.json", s)
    sample = read_group_sample(out / "trajectory.h5", len(x) - 1)
    assert sample["trajectory"]["positions_m"].shape[-2:] == (sim.n, 3)
    if s["status"] == "COMPLETE":
        np.testing.assert_allclose(
            result["physical_dt_s"].sum(), s["completed_physical_s"], rtol=1e-10, atol=1e-14
        )
    print(
        json.dumps(
            clean(
                dict(
                    run=out.name,
                    status=s["status"],
                    quality=s["numerical_quality"],
                    navigation=nav,
                    wall_s=s["wall_s"],
                )
            )
        ),
        flush=True,
    )
    return s
