"""Matched full-particle recording for original F6 and F10C execution.

All comparisons save the same 201 full32 states and diagnostics. Dynamic arrays
are never duplicated in place of an integration. Timing charges morphology/HDF5.
"""

from time import perf_counter
import json
import h5py
import numpy as np
import warp as wp
from mevo_mb.fast_sim.metrics import trajectory_metrics
from mevo_mb.fast_sim.morphodynamics import morphology_series, topology
from mevo_mb.fast_sim.f6_engine import F6Simulator
from mevo_mb.fast_sim.f10c_engine import FullParticleBatch
from .common import identity, write, clean


def read_prefix(array, count):
    shape = (count, *array.shape[1:])
    host = wp.empty(shape, dtype=array.dtype, device="cpu", pinned=True)
    if count:
        wp.copy(host, array, count=count * int(np.prod(array.shape[1:])))
        wp.synchronize()
    return host.numpy().copy()


def snapshot(s):
    shape = (s.batch, s.n)
    arrays = {
        "positions_m": s.x,
        "force_N": s.force,
        "free_velocity_m_s": s.free,
        "pressure_complex_Pa": s.pressure,
        "diag": s.diag,
        "active": s.active,
        "wall_gap_lower_bound_m": s.wall_gap,
        "wall_stats": s.wall_stats,
        "flow_velocity_m_s": s.flow_query.velocity,
        "flow_gradient_s_inv": s.flow_query.gradient,
        "flow_valid": s.flow_query.valid,
        "flow_tetra_id": s.flow_query.cache,
        "near_diag": s.near_diag,
        "free_path_m": s.free_path,
    }
    result = {}
    for name, array in arrays.items():
        value = array.numpy()
        result[name] = value.reshape(*shape, *value.shape[1:])
    return result


def run(
    out,
    x,
    events,
    obj,
    *,
    original=False,
    decomposition=False,
    duration=0.01,
    graphs=True,
    postprocess=None,
    sample_ids=None,
    dataset_metadata=None,
    identity_provider=None,
    dataset_options=None,
    physical_dt_s=1e-5,
    simulator_class=None,
    export_observer=False,
    substep_protocol=False,
    **extra,
):
    out.mkdir(parents=True, exist_ok=False)
    begin = perf_counter()
    p, b, r, g, f = obj
    xx = np.asarray(x)
    xx = xx[None] if xx.ndim == 2 else xx
    ids = list(sample_ids) if sample_ids is not None else [f"{out.name}_{i:04d}" for i in range(len(xx))]
    if len(ids) != len(xx) or len(set(ids)) != len(ids):
        raise ValueError("Unique sample identity required")
    source_identity = identity_provider or identity
    storage = (
        dict(compression="gzip", compression_opts=1) if dataset_options is None else dict(dataset_options)
    )
    dt = float(physical_dt_s)
    sample = 5e-5
    count = round(sample / dt)
    registered = (
        getattr(simulator_class, "registered_physical_steps", ()) if substep_protocol else (1e-5, 5e-6)
    )
    if dt not in registered or abs(count * dt - sample) > 1e-15:
        raise ValueError("Only registered 10/5 us steps with unchanged 50 us outputs")
    cls = F6Simulator if original else FullParticleBatch
    if "field_stride" in extra or "flow_stride" in extra:
        from mevo_mb.fast_sim.f10c_multirate import MultirateBatch

        cls = MultirateBatch
    if simulator_class is not None:
        if original or "field_stride" in extra or "flow_stride" in extra:
            raise ValueError("Explicit diagnostic class cannot override another execution mode")
        cls = simulator_class
    options = {} if original else dict(sample_ids=ids, decomposition=decomposition)
    options.update(extra)
    s = cls(xx, p, r, b, flow=f, vessel=g, graph_chunk=count, **options)
    if substep_protocol:
        s.configure_schedule(events, dt)
    records = [s.control(events[0])]
    prep = perf_counter() - begin
    t = perf_counter()
    s.graph(dt, 1 if substep_protocol else count)
    wp.synchronize()
    capture = perf_counter() - t
    np.testing.assert_array_equal(s.state(), xx)
    # Initial force fields are explicitly uncomputed; first physical-step diagnostics at t=50 us.
    t = perf_counter()
    frames = [snapshot(s)]
    readback = perf_counter() - t
    times = [0.0]
    control_indices = [0]
    integration = 0.0
    status = "COMPLETE"
    reason = None
    seg = 0
    steps = round(duration / sample)
    start_event = wp.Event(enable_timing=True)
    end_event = wp.Event(enable_timing=True)
    wp.record_event(start_event)
    for frame in range(steps):
        current = frame * sample
        if (
            not substep_protocol
            and seg + 1 < len(events)
            and abs(current - events[seg + 1]["time_s"]) < 1e-13
        ):
            seg += 1
            records.append(s.control(events[seg]))
            control_indices[-1] = seg
        t = perf_counter()
        try:
            s.advance(dt, count, graphs=graphs)
            if int(s.error.numpy()[0]):
                raise RuntimeError("GPU_ERROR_" + str(int(s.error.numpy()[0])))
        except RuntimeError as exc:
            integration += perf_counter() - t
            status = "FAILED"
            reason = str(exc)
            break
        integration += perf_counter() - t
        t = perf_counter()
        frames.append(snapshot(s))
        readback += perf_counter() - t
        times.append((frame + 1) * sample)
        control_indices.append(seg)
    if substep_protocol:
        records = s.executed_controls
        control_indices = [
            max(0, int(np.searchsorted([e["time_s"] for e in records], t + 1e-14, side="right") - 1))
            for t in times
        ]
    wp.record_event(end_event)
    wp.synchronize()
    span = wp.get_event_elapsed_time(start_event, end_event) / 1000
    failed_state = snapshot(s) if status == "FAILED" else None
    t = perf_counter()
    data = {key: np.asarray([frame[key] for frame in frames]) for key in frames[0]}
    num = int(s.event_count.numpy()[0])
    assert num <= s.event_time.shape[0]
    edge_ids = read_prefix(s.event_ids, num)
    edge_times = read_prefix(s.event_time, num)
    order = np.lexsort((edge_ids[:, 2], edge_ids[:, 1], edge_ids[:, 0], edge_times))
    edge_ids = edge_ids[order]
    edge_times = edge_times[order]
    qualities = [s.quality(i) for i in range(s.batch)]
    failures = s.failure_members() if status == "FAILED" and hasattr(s, "failure_members") else []
    postprocess_timing = None
    if postprocess is not None:
        metrics, lineage, morph, top, postprocess_timing = postprocess(
            data["positions_m"], p.radius_m, np.array(times), edge_ids, edge_times
        )
    else:
        metrics, lineage = trajectory_metrics(data["positions_m"], p.radius_m, np.array(times))
        morph = {}
        top = []
        for i in range(s.batch):
            m = morphology_series(data["positions_m"][:, i], np.arange(s.n), metrics["labels"][:, i])
            for key, value in m.items():
                morph.setdefault(key, []).append(value)
            ei = edge_ids[edge_ids[:, 0] == i].copy()
            et = edge_times[edge_ids[:, 0] == i]
            # Failed proposal events remain separate from accepted prefix topology.
            good = et <= times[-1] + 1e-14
            ei = ei[good]
            et = et[good]
            ei[:, 0] = 0
            tp, _, _ = topology(xx[i], np.array(times), metrics["labels"][:, i], et, ei, times[-1])
            top.append(tp)
        morph = {key: np.stack(value, axis=1) for key, value in morph.items()}
    topology_s = perf_counter() - t
    valid = (
        s.batch
        if status == "COMPLETE"
        and all(
            q["error_bits"] == 0
            and q["max_accepted_penetration_m"] <= 0.01 * p.radius_m
            and q["max_contact_pass_correction_m"] <= 0.1 * p.radius_m
            for q in qualities
        )
        else 0
    )
    meta = dict(
        batch_id=out.name,
        sample_id=ids,
        batch_size=s.batch,
        execution_mode="stable_reduction",
        graph_enabled=graphs,
        graph_chunk=1 if substep_protocol else count,
        field_update_stride=getattr(s, "field_stride", 1),
        flow_update_stride=getattr(s, "flow_stride", 1),
        effective_control_events=records,
        physical_dt_s=dt,
        per_sample_quality=[
            dict(sample_id=i, status="PASS" if valid else "ABORTED_BATCH", quality=q)
            for i, q in zip(ids, qualities)
        ],
        numerical_mode="FULL_PARTICLE_BASELINE" if original else "FULL_PARTICLE_BATCH",
        source=source_identity(),
        valid_sample_count=valid,
        training_ready=False,
        physics_scope="EXPLORATORY",
        labels=["BENCHMARK_ONLY", "PHYSICS_UNCALIBRATED"],
        split_group="F10A_translated_original_Y32",
        independent_scientific_initial_states=0,
        control_scope="Common synchronous control family",
        status=status,
        reason=reason,
        failures=failures,
        topology=top,
        force_timestamp="At saved t>0, force/pressure/free velocity refer to the last physical substep input at t-dt; positions are accepted t. Initial force fields are uncomputed zeros.",
        wall_gap_scope="Conservative Lipschitz lower bound from unchanged per-step final_geometry_check, not a new exact wall query.",
    )
    if dataset_metadata is not None:
        allowed = {
            "control_family_id",
            "initial_state_family",
            "split_group",
            "split_groups",
            "labels",
            "independent_scientific_initial_states",
            "postprocessing_mode",
            "original_physics_mode",
            "initial_state_sources",
        }
        if not set(dataset_metadata) <= allowed:
            raise ValueError("Dataset metadata cannot override physical state or quality")
        meta.update(dataset_metadata)
        meta["execution_status"] = status
        meta["numerical_quality"] = "NUMERICAL_PASS" if valid else "ABORTED_BATCH"
    if hasattr(s, "cache_stats"):
        meta["cache_statistics_columns"] = [
            "acoustic_reused",
            "acoustic_forced_refresh",
            "flow_reused",
            "guard_rejected",
            "flow_forced_refresh",
        ]
        meta["cache_statistics"] = s.cache_stats.numpy().reshape(s.batch, s.n, 5).sum(1).tolist()
        meta["CFD_final_proposal_query_stride"] = 1
        meta["graph_count"] = len(s.graphs)
        meta["total_graph_capture_s"] = s.capture_s
    if decomposition:
        meta["decomposition_max_residual"] = s.decomp_residual.numpy().max(0).tolist()
    t = perf_counter()
    with h5py.File(out / "trajectory.h5", "x") as h:
        h.attrs.update(schema="fast_pzt_f10c_v1", metadata_json=json.dumps(clean(meta)), training_ready=False)
        h["time_s"] = times
        h["control_index"] = control_indices
        h["particle_id"] = np.arange(s.n)
        h["radii_m"] = np.full(s.n, p.radius_m)
        h["physical_dt_s"] = np.full(round(times[-1] / dt), dt)
        for key, value in data.items():
            h.create_dataset("trajectory/" + key, data=value, **storage)
        h["positions_m"] = h["trajectory/positions_m"]
        if failed_state is not None:
            for key, value in failed_state.items():
                h.create_dataset("rejected_interval/" + key, data=value)
            proposal = s.final_proposal if substep_protocol else s.tmp1
            h["rejected_interval/proposal_m"] = proposal.numpy().reshape(s.batch, s.n, 3)
            h["rejected_interval"].attrs["scope"] = (
                "Rejected interval end/proposal, NOT accepted scientific trajectory"
            )
            h["rejected_interval"].attrs["start_s"] = times[-1]
            h["rejected_interval"].attrs["attempted_end_s"] = times[-1] + sample
        if substep_protocol:
            s.export_step_evidence(h)
        for key, value in metrics.items():
            h.create_dataset("metrics/" + key, data=value, **storage)
        for key, value in morph.items():
            h.create_dataset("morphology/" + key, data=value, **storage)
        for key, value in [
            ("source", meta["source"]),
            ("medium", r.metadata),
            ("flow", f.metadata),
            ("bubbles", p.metadata()),
            ("control", records),
            ("quality", meta["per_sample_quality"]),
        ]:
            h.require_group(key).attrs["metadata_json"] = json.dumps(clean(value))
        h["bubbles/particle_id"] = h["particle_id"]
        h["bubbles/radii_m"] = h["radii_m"]
        h["edge_events/batch_i_j_linked"] = edge_ids
        h["edge_events/time_s"] = edge_times
        h["boundary/label"] = s.exit_label.numpy().reshape(s.batch, s.n)
        h["boundary/time_s"] = s.exit_time.numpy().reshape(s.batch, s.n)
        if dataset_metadata is not None:
            h.create_dataset(
                "sampled_lineage_json", data=json.dumps(clean(lineage)), dtype=h5py.string_dtype()
            )
            for key in ["sample_id", "initial_state_family", "split_groups"]:
                h.create_dataset("source/" + key, data=meta[key], dtype=h5py.string_dtype())
    output = perf_counter() - t
    observer_output = 0.0
    if export_observer:
        t = perf_counter()
        s.export_observer(out)
        observer_output = perf_counter() - t
    timing = dict(
        model_prepare_s=prep,
        graph_capture_s=capture,
        integration_quality_wall_s=integration,
        sync_readback_wall_s=readback,
        cuda_event_trajectory_span_s=span,
        morphology_quality_output_prepare_s=topology_s,
        hdf5_output_s=output,
        resident_body_s=perf_counter() - begin,
        static_prepare_s=0.0,
        module_load_compile_observed_s=s.compile_s,
        valid_sample_count=valid,
    )
    if export_observer:
        timing["observer_output_s"] = observer_output
        meta["numerical_mode"] = s.numerical_mode
        meta["contact_iterations"] = p.contact_iterations
        meta["first_rejection_step"] = s.first_rejection
    timing["valid_samples_per_gpu_hour"] = 3600 * valid / timing["resident_body_s"]
    if postprocess_timing is not None:
        timing["postprocess_breakdown"] = postprocess_timing
    meta["effective_samples_per_gpu_hour"] = timing["valid_samples_per_gpu_hour"]
    with h5py.File(out / "trajectory.h5", "a") as h:
        h.attrs["metadata_json"] = json.dumps(clean(meta))
        h.attrs["effective_samples_per_gpu_hour"] = meta["effective_samples_per_gpu_hour"]
    write(out / "TIMING.json", timing)
    write(out / "SUMMARY.json", meta)
    write(
        out / "CONFIG.json",
        dict(
            initial_positions_m=xx,
            events=events,
            duration_s=duration,
            dt_s=dt,
            options=options,
            protocol_sha256=source_identity()["protocol_sha256"],
        ),
    )
    return timing
