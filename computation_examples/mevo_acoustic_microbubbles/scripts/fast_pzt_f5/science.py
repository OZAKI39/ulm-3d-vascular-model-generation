"""One registered continuous trajectory and its exact-state quiescent branch."""

import argparse
import copy
from time import perf_counter
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f5_engine import F5FlowSimulator
from mevo_mb.fast_sim.flow_gpu import GPUFlow
from mevo_mb.fast_sim.f4_engine import GPUVessel
from mevo_mb.fast_sim.models import ResponseSection, Basis
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim import kernels as k
from .common import BASE, INPUT, js, write, identity, sha
from .experiment import run_case


def measure(sim, event, fn, count=100):
    sim.reset()
    sim.control(event)
    sim.fields()
    sim.pairs()
    sim.tmp1.assign(sim.initial.reshape(-1, 3))
    sim._current_dt = 1e-5
    fn()
    wp.synchronize()
    with wp.ScopedCapture(device=sim.device) as captured:
        for _ in range(count):
            fn()
    start = wp.Event(device=sim.device, enable_timing=True)
    end = wp.Event(device=sim.device, enable_timing=True)
    wp.record_event(start)
    wp.capture_launch(captured.graph)
    wp.record_event(end)
    elapsed = wp.get_event_elapsed_time(start, end)
    assert sim.error.numpy()[0] == 0
    return dict(repetitions=count, GPU_mean_us=elapsed * 1000 / count)


def main(label, checks):
    assert js(BASE / checks / "RESULT.json")["status"] == "PASS"
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    start = perf_counter()
    wp.init()
    p = F2Physics()
    response = ResponseSection(INPUT / "f4/response.h5", p)
    basis = Basis(INPUT / "f4/basis.h5", p.frequency_Hz)
    g = GPUVessel(INPUT / "f4/prepared_geometry.npz", INPUT / "f4/labelled_lumen_m.vtp")
    flow = GPUFlow(INPUT / "result_121.vtu")
    cfg = js(BASE / "ACOUSTIC_RELEASE_PROTOCOL.json")
    _, result, metrics = run_case(out / "main", cfg, flow, g, basis, response)
    assert result["summary"]["numerical_quality"] == "PASS", "Do not branch from an invalid main run"
    assert result["active"].all(), "This registered experiment requires complete fixed-ID histories"
    switch = cfg["events"][1]["time_s"]
    index = int(np.flatnonzero(np.isclose(result["time_s"], switch, rtol=0, atol=1e-14))[0])
    x = result["positions_m"][index, 0].copy()
    switch_arrays = dict(
        time_s=np.array(switch),
        positions_m=x,
        particle_id=np.arange(len(x)),
        radii_m=np.full(len(x), p.radius_m),
        active=result["active"][index, 0],
        wall_gap_m=result["wall_gap_m"][index, 0],
        cluster_labels=metrics["labels"][index, 0],
        raw_quality=result["quality_samples"][index],
        wall_stats=result["wall_stats"][index],
        **{key: value[index] for key, value in result["flow"].items()},
    )
    np.savez_compressed(out / "SWITCH_STATE.npz", **switch_arrays)
    write(
        out / "SWITCH_STATE.json",
        dict(
            time_s=switch,
            frame_index=index,
            state_sha256=sha(out / "SWITCH_STATE.npz"),
            state_variables="Original positions, immutable particle IDs/radii, active mask, controls; overdamped motion has no inertial velocity state",
            reproducible_work_arrays="Mobility/forces solved afresh each step; static wall/flow caches are rebuilt without changing positions. Adjacency rebuilt from the exact switch positions.",
            actual_controls_before=cfg["events"][0],
            actual_controls_after=cfg["events"][1],
        ),
    )
    # A is the uninterrupted computed main suffix: no simulated restart or
    # synthetic animation segment is substituted for its physical history.
    suffix = dict(time_s=result["time_s"][index:] - switch, absolute_time_s=result["time_s"][index:])
    for key in ["positions_m", "active", "wall_gap_m", "wall_stats", "quality_samples"]:
        suffix[key] = result[key][index:]
    for key, val in result["flow"].items():
        suffix["flow_" + key] = val[index:]
    edge_mask = result["edge_event_time_s"] > switch + 1e-14
    suffix["edge_event_time_s"] = result["edge_event_time_s"][edge_mask] - switch
    suffix["edge_event_ids"] = result["edge_event_ids"][edge_mask]
    np.savez_compressed(out / "branch_A_actual_suffix.npz", **suffix)
    write(
        out / "BRANCH_A.json",
        dict(
            mode="ULTRASOUND_OFF_FLOW_ON",
            main_sha256=sha(out / "main/trajectory.h5"),
            main_frame_start=index,
            continuous_suffix=True,
            independent_recomputation=False,
            physical_steps=int(round((cfg["duration_s"] - switch) / cfg["dt_s"])),
        ),
    )
    b = copy.deepcopy(cfg)
    b.update(positions_m=x.tolist(), duration_s=cfg["duration_s"] - switch)
    b["events"] = [dict(cfg["events"][1], time_s=0.0, flow_enabled=False, mode="ULTRASOUND_OFF_FLOW_OFF")]
    _, rb, _ = run_case(out / "branch_B", b, flow, g, basis, response)
    assert rb["summary"]["numerical_quality"] == "PASS"
    np.testing.assert_array_equal(rb["positions_m"][0, 0], x)
    np.testing.assert_array_equal(rb["active"][0, 0], result["active"][index, 0])

    # One bounded timing session, independent of the scientific trajectories.
    sim = F5FlowSimulator(x, p, response, basis, flow=flow, vessel=g)
    on = cfg["events"][0]
    off = dict(cfg["events"][1], time_s=0.0)

    def query():
        flow.launch(sim.x, sim.active, sim.pending, sim.flow_query, sim.error)

    timings = {"flow_query_cached_32": measure(sim, on, query)}

    def uncached():
        sim.flow_query.cache.fill_(-1)
        query()

    timings["flow_query_BVH_including_cache_reset_32"] = measure(sim, on, uncached)
    timings["primary_field_response_32"] = measure(sim, on, sim.fields)
    timings["secondary_pair_force_32"] = measure(sim, on, sim.pairs)

    def contact():
        sim.launch(
            k.contact,
            sim.total,
            [sim.tmp1, sim.tmp2, sim.n, p.radius_m, p.contact_relaxation, sim.diag, sim.error],
        )

    timings["one_pair_contact_pass_and_inherited_wall_pass_32"] = measure(sim, on, contact)
    timings["full_ON_step_32"] = measure(sim, on, lambda: sim.step(cfg["dt_s"]))
    timings["full_OFF_step_32"] = measure(sim, off, lambda: sim.step(cfg["dt_s"]))
    performance = dict(
        status="PASS",
        identity=identity(),
        GPU=str(wp.get_device("cuda:0")),
        GPU_name=wp.get_device("cuda:0").name,
        warp=wp.__version__,
        dt_s=cfg["dt_s"],
        static_flow_build_s=flow.metadata["static_build_s"],
        main=result["summary"],
        branch_B=rb["summary"],
        CUDA_event_timings=timings,
        scope="Resident warmed CUDA events; component timings are not additive. Full-run wall time includes saved-frame transfers, checks and graph capture; not pure GPU time.",
        branch_A="Already included in main physical steps and runtime; never double-counted",
        actual_main_macro_steps=60000,
        actual_main_physical_substeps=60000,
        actual_B_macro_steps=30000,
        actual_B_physical_substeps=30000,
        scientific_task_body_wall_s=perf_counter() - start,
    )
    write(out / "PERFORMANCE.json", performance)
    write(BASE / "PERFORMANCE.json", performance)
    write(
        out / "RESULT.json",
        dict(
            status="COMPLETE",
            numerical_quality="PASS",
            identity=identity(),
            main="main",
            branch_A="branch_A_actual_suffix.npz",
            branch_B="branch_B",
            switch_index=index,
            switch_s=switch,
            optional_C="Not needed for the registered flow/no-flow causal comparison",
        ),
    )
    print("F5 registered experiment COMPLETE", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--checks", required=True)
    args = parser.parse_args()
    main(args.label, args.checks)
