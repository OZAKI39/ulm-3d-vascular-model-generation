"""Complete timing after the first job's post-trajectory instrumentation error.

Both original successful trajectories remain byte-for-byte unchanged.
"""

import argparse
import h5py
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f5_engine import F5FlowSimulator
from mevo_mb.fast_sim.flow_gpu import GPUFlow
from mevo_mb.fast_sim.f4_engine import GPUVessel
from mevo_mb.fast_sim.models import ResponseSection, Basis
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim import kernels as k
from .science import measure
from .common import BASE, INPUT, js, write, identity, sha


def main(label, science):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    source = BASE / science
    files = [source / name / "trajectory.h5" for name in ["main", "branch_B"]]
    before = {str(f): sha(f) for f in files}
    summaries = {name: js(source / name / "SUMMARY.json") for name in ["main", "branch_B"]}
    assert all(s["numerical_quality"] == "PASS" and s["status"] == "COMPLETE" for s in summaries.values())
    switch = js(source / "SWITCH_STATE.json")
    x = np.load(source / "SWITCH_STATE.npz")["positions_m"]
    with h5py.File(files[0]) as a, h5py.File(files[1]) as b:
        np.testing.assert_array_equal(a["positions_m"][switch["frame_index"], 0], x)
        np.testing.assert_array_equal(b["positions_m"][0, 0], x)
        assert a["time_s"][-1] == 0.6 and b["time_s"][-1] == 0.3
    cfg = js(source / "main/CONFIG.json")
    wp.init()
    p = F2Physics()
    response = ResponseSection(INPUT / "f4/response.h5", p)
    basis = Basis(INPUT / "f4/basis.h5", p.frequency_Hz)
    geometry = GPUVessel(INPUT / "f4/prepared_geometry.npz", INPUT / "f4/labelled_lumen_m.vtp")
    flow = GPUFlow(INPUT / "result_121.vtu")
    sim = F5FlowSimulator(x, p, response, basis, flow=flow, vessel=geometry)
    on, off = cfg["events"][0], dict(cfg["events"][1], time_s=0.0)

    def query():
        flow.launch(sim.x, sim.active, sim.pending, sim.flow_query, sim.error)

    def uncached():
        sim.flow_query.cache.fill_(-1)
        query()

    def contact():
        sim._pair_contact_pass = 2
        sim.launch(
            k.contact,
            sim.total,
            [sim.tmp1, sim.tmp2, sim.n, p.radius_m, p.contact_relaxation, sim.diag, sim.error],
        )

    timings = {}
    for name, event, fn in [
        ("flow_query_cached_32", on, query),
        ("flow_query_BVH_including_cache_reset_32", on, uncached),
        ("primary_field_response_32", on, sim.fields),
        ("secondary_pair_force_32", on, sim.pairs),
        ("one_pair_contact_pass_and_wall_pass_32", on, contact),
        ("full_ON_step_32", on, lambda: sim.step(cfg["dt_s"])),
        ("full_OFF_step_32", off, lambda: sim.step(cfg["dt_s"])),
    ]:
        timings[name] = measure(sim, event, fn)
    assert before == {str(f): sha(f) for f in files}
    result = dict(
        status="PASS",
        identity=identity(),
        science_source=before,
        trajectories_unchanged=True,
        GPU_name=wp.get_device("cuda:0").name,
        warp=wp.__version__,
        dt_s=cfg["dt_s"],
        static_flow_build_s=flow.metadata["static_build_s"],
        CUDA_event_timings=timings,
        **summaries,
        scope="CUDA event component timings are not additive. Full run wall_s includes saved-frame transfers, checks, graph captures and host work; it is not pure GPU compute time.",
        branch_A="Continuous suffix of main; no duplicate physical steps or runtime",
        main_macro_steps=summaries["main"]["macro_step_count"],
        main_physical_substeps=summaries["main"]["actual_substep_count"],
        B_macro_steps=summaries["branch_B"]["macro_step_count"],
        B_physical_substeps=summaries["branch_B"]["actual_substep_count"],
        recovery="Only performance instrumentation was rerun; original science_v1 main and B completed before its timing exception",
    )
    write(out / "RESULT.json", result)
    write(BASE / "PERFORMANCE.json", result)
    write(
        source / "RESULT.json",
        dict(
            status="COMPLETE",
            numerical_quality="PASS",
            identity=js(source / "main/IDENTITY.json"),
            main="main",
            branch_A="branch_A_actual_suffix.npz",
            branch_B="branch_B",
            switch_index=switch["frame_index"],
            switch_s=0.3,
            performance_recovery=label,
            original_process="Exited after completed trajectories due to post-run timing initialization; original logs retained",
        ),
    )
    print("F5 timing PASS; completed science trajectories unchanged", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--science", required=True)
    a = ap.parse_args()
    main(a.label, a.science)
