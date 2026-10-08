"""One matched N200/5ms regression and bounded resident GPU wall microtimings."""

import argparse
import h5py
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f4_engine import GPUVessel, F4VesselSimulator
from mevo_mb.fast_sim import vessel_kernels as vk
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from scripts.fast_pzt.performance import event_measure
from .common import BASE, ROOT, js, write, identity
from .experiment import run_case


def main(label):
    assert js(BASE / "aggregation_v1/RESULT.json")["status"] == "COMPLETE"
    assert js(BASE / "controls_v1/RESULT.json")["status"] == "COMPLETE"
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    with h5py.File(ROOT / "results/fast_pzt/f3/performance_v1/stable_reduction/trajectory.h5") as h:
        initial = h["positions_m"][0]
    cfg = dict(
        positions_m=initial.tolist(),
        duration_s=0.005,
        dt_s=1e-5,
        sample_interval_s=0.0005,
        basis=False,
        events=[dict(time_s=0, pressure_Pa=1000)],
        wall_budget_s=60,
        execution_mode="stable_reduction",
        flow_enabled=False,
        scope="PERFORMANCE_ONLY; original F3 N200 initial state and uniform 1000 Pa",
    )
    g = GPUVessel(BASE / "inputs/prepared_geometry.npz", BASE / "inputs/labelled_lumen_m.vtp")
    rows = {}
    for enabled in [False, True]:
        name = "wall_on" if enabled else "wall_off"
        rows[name] = run_case(out / name, dict(cfg, wall_enabled=enabled), g)
        assert rows[name]["numerical_quality"] == "PASS"
    with h5py.File(out / "wall_on/trajectory.h5") as a, h5py.File(out / "wall_off/trajectory.h5") as b:
        np.testing.assert_array_equal(a["positions_m"][:], b["positions_m"][:])
    p = F2Physics()
    response = ResponseSection(BASE / "inputs/response.h5", p)
    x = np.load(BASE / "inputs/parent_initial.npz")["positions_m"]
    sim = F4VesselSimulator(x, p, response, vessel=g)

    def query():
        sim.launch(
            vk.query_points,
            sim.total,
            [*g.args(), sim.x, sim.query_d, sim.query_n, sim.query_face, sim.query_inside],
        )

    timings = {"parent_32_exact_distance_query": event_measure(sim, query)}
    sim.tmp1.assign(sim.initial.reshape(-1, 3))
    timings["parent_32_far_clearance_wall_pass"] = event_measure(sim, lambda: sim.wall_pass(sim.tmp1))
    # Recreate the same 1 nm overlap on a verified real finite face every call;
    # includes reset-device-copy, query, refinement and projection. No subtraction.
    region = js(BASE / "GEOMETRY_FIELD_OVERLAP.json")["initialization"]["parent"]
    anchor = np.array(region["wall_anchor_m"])
    normal = np.array(region["outward_normal"])
    one = F4VesselSimulator((anchor - normal * (p.radius_m + 1e-9))[None], p, response, vessel=g)
    candidate = wp.array((anchor - normal * (p.radius_m - 1e-9))[None], dtype=wp.vec3d, device=one.device)

    def correction():
        wp.copy(one.tmp1, candidate)
        one.wall_pass(one.tmp1)

    timings["one_contact_query_and_projection_including_candidate_copy"] = event_measure(one, correction)
    assert one.quality()["wall"]["correction_passes"] == 100
    assert one.error.numpy()[0] == 0
    write(
        out / "RESULT.json",
        dict(
            status="PASS",
            identity=identity(),
            rows=rows,
            GPU_timings=timings,
            wall_on_off_bitwise_equal=True,
            N200_wall_ratio=rows["wall_on"]["wall_s"] / rows["wall_off"]["wall_s"],
            physical_steps_each=500,
            additional_substeps=0,
            scope="One N200/5ms pair; no new physical validation. Timings include device-copy where stated and are not additive. Distance query and contact projection share a kernel; no invented separate projection cost.",
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    main(ap.parse_args().label)
