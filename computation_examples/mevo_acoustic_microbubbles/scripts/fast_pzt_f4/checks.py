"""Bounded actual-mesh geometry tests; imposed velocities are never acoustic data."""

import argparse
import numpy as np
from mevo_mb.fast_sim.f4_engine import GPUVessel, F4VesselSimulator
from mevo_mb.fast_sim.f3_engine import F3Simulator
from mevo_mb.fast_sim.models import ResponseSection, Basis
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import ROOT, BASE, js, write, identity


def main(label):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    p = F2Physics()
    r = p.radius_m
    response = ResponseSection(BASE / "inputs/response.h5", p)
    basis = Basis(BASE / "inputs/basis.h5", p.frequency_Hz)
    vessel = GPUVessel(BASE / "inputs/prepared_geometry.npz", BASE / "inputs/labelled_lumen_m.vtp")
    reference = np.load(BASE / "inputs/geometry_validation_reference.npz")
    gpu = vessel.query(reference["positions_m"])
    discrepancy = abs(gpu["distance_m"] - reference["distance_m"])
    np.savez_compressed(
        out / "query_accuracy.npz",
        **gpu,
        positions_m=reference["positions_m"],
        reference_distance_m=reference["distance_m"],
    )
    assert discrepancy.max() < 0.01 * r, discrepancy.max()
    reliable = abs(reference["distance_m"]) > 1e-8
    assert np.array_equal(gpu["inside"][reliable], reference["inside"][reliable])
    rows = {}
    region = js(BASE / "GEOMETRY_FIELD_OVERLAP.json")["initialization"]["parent"]
    anchor = np.array(region["wall_anchor_m"])
    normal = np.array(region["outward_normal"])
    tangent = np.array(region["tangent_1"])

    def run(name, x, velocity, duration=0.001, field=None, wall=True, expected="COMPLETE"):
        sim = F4VesselSimulator(
            x,
            p,
            response,
            field,
            vessel=vessel,
            wall_contact_enabled=wall,
            geometry_unit_test_velocity=velocity,
        )
        event = (
            dict(time_s=0, pressure_Pa=0)
            if field is None
            else dict(time_s=0, amplitude_Pa=[0] * 4, phase_rad=[0] * 4, enabled=[False] * 4)
        )
        data = sim.run(duration, 1e-5, [event], 1e-5, wall_budget_s=60)
        np.savez_compressed(
            out / (name + ".npz"),
            time_s=data["time_s"],
            positions_m=data["positions_m"],
            active=data["active"],
            wall_gap_m=data["wall_gap_m"],
            wall_stats=data["wall_stats"],
        )
        rows[name] = dict(
            summary=data["summary"],
            boundary_events=data["boundary_events"],
            initial_m=x,
            final_m=data["positions_m"][-1, 0],
            test_velocity_m_s=velocity,
            scope="GEOMETRY_UNIT_TEST",
        )
        write(out / "PARTIAL.json", dict(identity=identity(), rows=rows))
        assert data["summary"]["status"] == expected, (name, data["summary"]["reason"])
        if expected == "COMPLETE":
            assert data["summary"]["quality"]["max_accepted_penetration_m"] <= 0.01 * r
            assert data["summary"]["quality"]["max_contact_pass_correction_m"] <= 0.1 * r
            assert data["summary"]["quality"]["wall"]["max_pass_correction_m"] <= 0.1 * r
            if wall:
                assert data["wall_gap_m"][data["active"].astype(bool)].min() >= -0.01 * r
        return data

    x = (anchor - normal * (r + 1e-9))[None]
    hit = run("wall_approach", x, normal * 1e-4)
    assert hit["summary"]["quality"]["wall"]["correction_passes"] > 0
    slide = run("wall_slide", x, tangent * 1e-4)
    moved = slide["positions_m"][-1, 0, 0] - x[0]
    assert np.dot(moved, tangent) > 0.99 * 1e-4 * 0.001
    pair = np.array(
        [
            anchor - normal * (r + 2e-9) - tangent * (r + 0.15e-6),
            anchor - normal * (r + 2e-9) + tangent * (r + 0.15e-6),
        ]
    )
    velocities = np.array([normal * 5e-5 + tangent * 3e-4, normal * 5e-5 - tangent * 3e-4])
    run("pair_wall_coupling", pair, velocities)
    angle = np.pi / 6
    cloud = np.array([[0.004 * np.cos(angle), -0.004 * np.sin(angle), 0]])
    run("reject_cross_bifurcation_wall", cloud, [0, 0.004 / 1e-5, 0], duration=1e-5, expected="FAILED")
    planes = js(BASE / "inputs/PLANES.json")
    for idx, name in enumerate(["INLET", "OUTLET_A", "OUTLET_B"], start=2):
        o = np.array(planes[name]["origin_m"])
        n = np.array(planes[name]["normal"])
        d = run("exit_" + name, (o - n * 20e-9)[None], n * 1e-3, duration=5e-5)
        assert len(d["boundary_events"]) == 1 and d["boundary_events"][0]["label"] == idx
        assert abs(d["boundary_events"][0]["time_s"] - 2e-5) < 2e-12
        assert not d["active"][-1].any()
        np.testing.assert_allclose(d["positions_m"][-1, 0, 0], o, atol=1e-14, rtol=0)
    bounds = np.array(js(BASE / "GEOMETRY_FIELD_OVERLAP.json")["usable_center_bounds_m"])
    d = run(
        "field_domain_exit",
        np.array([[bounds[0, 0] + 20e-9, 0, 0]]),
        [-1e-3, 0, 0],
        duration=5e-5,
        field=basis,
    )
    assert len(d["boundary_events"]) == 1 and d["boundary_events"][0]["event"] == "FIELD_DOMAIN_EXIT"
    assert not d["active"][-1].any()

    old = ROOT / "results/fast_pzt/f2"
    mature = np.load(old / "mature_v1/MATURE_STATE.npz")["positions_m"]
    control = js(old / "mature_v1/FIELD_SELECTION.json")["baseline"]
    f3 = F3Simulator(mature, p, response, basis).run(0.0005, 1e-5, [control], 0.0001)
    f4 = F4VesselSimulator(mature, p, response, basis, vessel_enabled=False).run(
        0.0005, 1e-5, [control], 0.0001
    )
    np.testing.assert_array_equal(f3["positions_m"], f4["positions_m"])
    np.testing.assert_array_equal(f3["edge_event_ids"], f4["edge_event_ids"])
    np.savez_compressed(
        out / "disabled_equivalence.npz",
        f3_positions=f3["positions_m"],
        f4_positions=f4["positions_m"],
        time_s=f3["time_s"],
    )
    # Noninteracting inactive ID must not change the remaining bubble's motion.
    active_test = F4VesselSimulator(
        np.array([anchor - normal * (r + 5e-6), anchor - normal * (r + 10e-6)]), p, response, vessel=vessel
    )
    active_test.active.assign(np.array([1, 0], np.int32))
    active_test.control(dict(time_s=0, pressure_Pa=1000))
    active_test.step(1e-5)
    np.testing.assert_array_equal(active_test.state()[0], active_test.initial[0])
    write(
        out / "RESULT.json",
        dict(
            status="PASS",
            identity=identity(),
            query_max_abs_error_m=float(discrepancy.max()),
            query_points=len(discrepancy),
            query_vs_VTK_RMS_m=float(np.sqrt(np.mean(discrepancy**2))),
            signs_match_away_from_surface=True,
            disabled_vessel_F3_bitwise=True,
            inactive_particle_excluded=True,
            rows=rows,
            scope="GEOMETRY_UNIT_TEST; prescribed test velocities are not acoustic trajectories",
        ),
    )
    print("F4 geometry checks PASS", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    main(a.label)
