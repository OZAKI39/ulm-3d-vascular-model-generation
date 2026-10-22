"""Necessary geometry, response, zero-flow and GPU wall negative controls."""

import copy
import numpy as np
from mevo_mb.coupling.straight_tube import StraightTube, apparatus, geometry_audit, ZeroFlow
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import BASE, ROOT, configs, identity, write


def must_fail(fn):
    try:
        fn()
    except (ValueError, RuntimeError):
        return True
    raise AssertionError("Invalid input unexpectedly admitted")


def run():
    c = configs()
    g = c["demo"]["geometry"]
    f = c["acoustic"]["acoustics"]["frequency_Hz"]
    tube, elements = apparatus(g, f)
    results = {}
    assert tube.contains(tube.world([[0, 0, 0], [0, tube.radius_m - 1.78e-6, 0]]), 1.78e-6).all()
    assert not tube.contains(tube.world([[0, tube.radius_m, 0]]), 1.78e-6).any()
    assert not tube.contains(tube.world([[tube.length_m / 2 + 1e-9, 0, 0]]), 1.78e-6).any()
    results["finite_radius_radial_and_axial_guards"] = True
    rot = copy.deepcopy(g)
    rot.update(
        rotation_local_to_world=[[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        center_m=[0.002, -0.001, 0.003],
        end_sign=-1,
    )
    other, ee = apparatus(rot, f)
    assert geometry_audit(other, ee, g["gap_m"])["status"] == "PASS"
    points = np.array([[0, 0, 0], [0.005, 0.0001, 0]])
    np.testing.assert_allclose(other.local(other.world(points)), points, atol=1e-17)
    results["rotation_translation_and_end_flip"] = True
    bad = copy.deepcopy(g)
    bad["pzt_overrides"] = [dict(position_m=e["position_m"], normal=e["normal"]) for e in elements]
    bad["pzt_overrides"][0]["normal"] = [0, 1, 0]
    results["reversed_source_rejected"] = must_fail(lambda: geometry_audit(*apparatus(bad, f), g["gap_m"]))
    badz = copy.deepcopy(g)
    badz["pzt_overrides"] = copy.deepcopy(bad["pzt_overrides"])
    badz["pzt_overrides"][0] = dict(position_m=[0, 0, 0.0039], normal=[0, 0, -1])
    results["upper_source_rejected"] = must_fail(lambda: geometry_audit(*apparatus(badz, f), g["gap_m"]))
    badwall = copy.deepcopy(g)
    badwall["wall_thickness_m"] = 0
    results["unknown_wall_not_silently_zero"] = (
        c["lab_parameters"]["parameters"]["wall_thickness_m"]["value"] is None
    )
    results["zero_demo_wall_rejected"] = must_fail(lambda: StraightTube.from_config(badwall))
    rng = np.random.default_rng(111)
    q = rng.uniform(-1, 1, (4096, 3)) * [tube.length_m / 2, tube.radius_m / 2, tube.radius_m / 2]
    velocity, gradient = ZeroFlow().query(tube.world(q))
    assert np.count_nonzero(velocity) == np.count_nonzero(gradient) == 0
    results["zero_flow_4096_positions"] = True
    response = ResponseSection(ROOT / c["bubbles"]["response_path"], F2Physics())
    for name, value in [("above", 40000.01), ("below", 0.01)]:
        results["response_" + name + "_rejected"] = must_fail(
            lambda value=value: response.query(np.array([value + 0j]))
        )
    results["frequency_mismatch_rejected"] = must_fail(
        lambda: ResponseSection(ROOT / c["bubbles"]["response_path"], F2Physics(frequency_Hz=2e6))
    )
    assert np.count_nonzero(response.query(np.array([0j]))) == 0
    results["zero_pressure_exact_quiescent"] = True
    try:
        import warp as wp
    except ImportError:
        write(BASE / "analysis/CPU_CHECKS.json", dict(status="PASS", checks=results, identity=identity()))
        return
    from mevo_mb.fast_sim.straight_tube_engine import cylinder_project, GPUZeroFlow
    from mevo_mb.fast_sim.flow_gpu import QueryBuffers

    wp.init()
    dev = "cuda:0"
    rad = 1.78e-6

    def project(points, tube):
        x = wp.array(tube.world(points), dtype=wp.vec3d, device=dev)
        stats = wp.zeros((len(points), 5), dtype=wp.float64, device=dev)
        error = wp.zeros(1, dtype=int, device=dev)
        wp.launch(
            cylinder_project,
            dim=len(points),
            inputs=[
                x,
                wp.vec3d(*tube.center),
                wp.mat33d(*tube.rotation.ravel()),
                tube.radius_m,
                tube.length_m / 2,
                rad,
                int(tube.ends[0] == "closed"),
                int(tube.ends[1] == "closed"),
                stats,
                error,
            ],
            device=dev,
        )
        return tube.local(x.numpy()), stats.numpy(), int(error.numpy()[0])

    x, stats, err = project([[0, tube.radius_m - rad + 0.05 * rad, 0]], tube)
    assert err == 0 and abs(x[0, 1] - (tube.radius_m - rad)) < 1e-16
    assert np.isclose(stats[0, 2], 0.05 * rad, rtol=1e-8)
    results["wall_small_projection_and_record"] = True
    _, _, err = project([[0, tube.radius_m - rad + 0.11 * rad, 0]], tube)
    assert err & 1024
    results["wall_original_0p1R_budget_rejects"] = True
    _, _, err = project([[tube.length_m / 2 + 0.01 * rad, 0, 0]], tube)
    assert err & 32768
    results["open_end_stops_without_reflection"] = True
    closed = copy.deepcopy(tube)
    closed.ends = ("closed", "closed")
    xx, _, err = project([[tube.length_m / 2 - rad + 0.05 * rad, 0, 0]], closed)
    assert err == 0 and np.isclose(xx[0, 0], tube.length_m / 2 - rad, rtol=0, atol=1e-16)
    results["closed_end_finite_radius_constraint"] = True
    xr, _, err = project([[0, other.radius_m - rad + 0.05 * rad, 0]], other)
    np.testing.assert_allclose(xr, x, atol=1e-17)
    assert err == 0
    results["GPU_wall_rotation_equivariance"] = True
    gpu_x = wp.array(tube.world(q), dtype=wp.vec3d, device=dev)
    active = wp.ones(len(q), dtype=int, device=dev)
    pending = wp.zeros_like(active)
    error = wp.zeros(1, dtype=int, device=dev)
    buffers = QueryBuffers(len(q), dev)
    GPUZeroFlow().launch(gpu_x, active, pending, buffers, error)
    assert np.count_nonzero(buffers.velocity.numpy()) == np.count_nonzero(buffers.gradient.numpy()) == 0
    assert buffers.valid.numpy().all()
    results["GPU_zero_flow_4096_exact"] = True
    from mevo_mb.acoustics.kwave_backend import extract_complex

    times = np.arange(333, 333 + 96) / (24 * f)
    phasor = np.array([3 + 4j, -2 + 1j])
    raw = np.real(phasor[:, None] * np.exp(2j * np.pi * f * times))
    np.testing.assert_allclose(extract_complex(raw, times, f), phasor, rtol=1e-6, atol=1e-6)
    results["absolute_time_positive_phasor_convention"] = True
    from mevo_mb.fast_sim.models import Basis
    from mevo_mb.fast_sim.f6_engine import F6Simulator
    from mevo_mb.fast_sim.straight_tube_engine import StraightTubeSimulator

    basis = Basis(BASE / "acoustic/basis.h5", f)
    initial = np.load(BASE / "initial_state.npz")["positions_m"]
    event = dict(c["acoustic"]["demo_control"])
    old = F6Simulator(
        initial, F2Physics(), response, basis, vessel_enabled=False, flow=None, flow_enabled=False
    )
    new = StraightTubeSimulator(initial, F2Physics(), response, basis, tube=tube)
    for sim in (old, new):
        sim.control(event)
        sim.advance(1e-6, 10, graphs=False)
    np.testing.assert_array_equal(old.state(), new.state())
    results["far_wall_M0_ten_steps_bitwise_legacy_equivalent"] = True
    results["four_source_control_rejected"] = must_fail(
        lambda: new.control(dict(event, amplitude_Pa=[1] * 4))
    )
    results["nonzero_background_flow_rejected"] = must_fail(
        lambda: new.control(dict(event, flow_enabled=True))
    )
    write(BASE / "analysis/ENGINE_CHECKS.json", dict(status="PASS", checks=results, identity=identity()))
    print(results, flush=True)


if __name__ == "__main__":
    run()
