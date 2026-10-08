"""Small preregistered interpolation, coupling, and domain-semantics checks."""

import argparse
import numpy as np
from scipy.linalg import expm
import warp as wp
from mevo_mb.fast_sim.flow_gpu import GPUFlow, AffineFlow
from mevo_mb.fast_sim.f5_engine import F5FlowSimulator
from mevo_mb.fast_sim.f4_engine import F4VesselSimulator, GPUVessel
from mevo_mb.fast_sim.models import ResponseSection, Basis
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.flow.fields import SimVascularFlowField
from scripts.fast_pzt_f3.snapshot import evaluate
from .common import BASE, INPUT, js, write, identity


def main(label):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    wp.init()
    p = F2Physics()
    response = ResponseSection(INPUT / "f4/response.h5", p)
    basis = Basis(INPUT / "f4/basis.h5", p.frequency_Hz)
    geometry = GPUVessel(INPUT / "f4/prepared_geometry.npz", INPUT / "f4/labelled_lumen_m.vtp")
    flow = GPUFlow(INPUT / "result_121.vtu")
    cpu = SimVascularFlowField(INPUT / "result_121.vtu", spatial_method="tetra_p1")
    groups = np.load(INPUT / "flow_validation_points.npz")
    interpolation = {}
    for name in groups.files:
        x = groups[name]
        gpu = flow.sample(x)
        ref = cpu.velocity(x, 0)
        assert gpu["valid"].all() and gpu["error"] == 0
        error = float(np.linalg.norm(gpu["velocity_m_s"] - ref, axis=1).max())
        assert error < 2e-12
        grad = (
            np.swapaxes(
                cpu.sampler.values[gpu["tetra_id"], 1:] - cpu.sampler.values[gpu["tetra_id"], :1], 1, 2
            )
            @ cpu.sampler.inverse[gpu["tetra_id"]]
        )
        grad_error = float(np.abs(grad - gpu["gradient_s_inv"]).max())
        assert grad_error < 1e-10
        fd = np.empty_like(grad)
        for a in range(3):
            dx = np.eye(3)[a] * 1e-9
            fd[:, :, a] = (cpu.velocity(x + dx, 0) - cpu.velocity(x - dx, 0)) / (2e-9)
        fd_error = float(np.abs(fd - gpu["gradient_s_inv"]).max())
        assert fd_error < 1e-6
        interpolation[name] = dict(
            points=len(x),
            max_velocity_abs_error_m_s=error,
            max_gradient_abs_error_s_inv=grad_error,
            CPU_finite_difference_gradient_error_s_inv=fd_error,
        )
        np.savez_compressed(out / (name + "_interpolation.npz"), positions_m=x, CPU_velocity_m_s=ref, **gpu)
    outside = flow.sample(np.array([[0.0, 0.0, 0.003]]))
    assert not outside["valid"].any() and outside["error"] & 16384
    validation = dict(
        status="PASS",
        identity=identity(),
        flow=flow.metadata,
        groups=interpolation,
        outside_rejected=True,
        scope="Interpolation of stored CFD values, not CFD mesh/physical validation",
    )
    write(out / "FLOW_GPU_VALIDATION.json", validation)
    write(BASE / "FLOW_GPU_VALIDATION.json", validation)

    rows = {}
    zero = dict(
        time_s=0.0, pressure_Pa=0.0, flow_enabled=True, secondary_enabled=False, mode="ULTRASOUND_OFF_FLOW_ON"
    )

    def execute(name, x, provider, events=None, duration=0.002, real_basis=False, vessel=False):
        e = events or [zero]
        sim = F5FlowSimulator(
            x,
            p,
            response,
            basis if real_basis else None,
            flow=provider,
            acoustic_enabled=bool(real_basis and any(e[0].get("enabled", []))),
            flow_enabled=e[0]["flow_enabled"],
            vessel=geometry,
            vessel_enabled=vessel,
        )
        result = sim.run(duration, 1e-5, e, 0.0001, 90)
        assert result["summary"]["status"] == "COMPLETE", result["summary"]
        np.savez_compressed(out / (name + ".npz"), time_s=result["time_s"], positions_m=result["positions_m"])
        rows[name] = dict(summary=result["summary"])
        return sim, result

    uniform = AffineFlow([0.002, -0.001, 0.0005])
    for n in [1, 4]:
        x = np.zeros((n, 3))
        x[:, 0] = -0.002
        x[:, 1] = np.arange(n) * 2.02 * p.radius_m
        sim, r = execute("uniform_" + str(n), x, uniform)
        expected = x[None] + r["time_s"][:, None, None] * np.array([0.002, -0.001, 0.0005])
        err = float(np.abs(r["positions_m"][:, 0] - expected).max())
        rel = r["positions_m"][:, 0] - r["positions_m"][:, 0].mean(1)[:, None]
        rg = np.sqrt(np.mean(np.sum(rel * rel, axis=-1), axis=1))
        rgerr = float(np.abs(rg - rg[0]).max())
        assert err < 5e-14 and rgerr < 5e-14
        rows["uniform_" + str(n)].update(max_position_error_m=err, max_rg_change_m=rgerr)
    shear = np.zeros((3, 3))
    shear[0, 1] = 8.0
    x = np.array([[-0.002, -1e-5, 0], [-0.002, 1e-5, 0]])
    _, r = execute("linear_shear", x, AffineFlow([0.001, 0, 0], shear))
    expected = x[None] + r["time_s"][:, None, None] * np.array([[0.001 - 8e-5, 0, 0], [0.001 + 8e-5, 0, 0]])
    err = float(np.abs(r["positions_m"][:, 0] - expected).max())
    assert err < 5e-14
    rows["linear_shear"].update(
        max_position_error_m=err,
        relative_displacement_m=(
            r["positions_m"][-1, 0, 1] - r["positions_m"][-1, 0, 0] - (x[1] - x[0])
        ).tolist(),
    )

    # Euler integration error against the exact local affine P1 solution, only
    # when the short path remains in its initial tetrahedron.
    x = groups["parent"][6:7]
    sample = flow.sample(x)
    cell = int(sample["tetra_id"][0])
    G = sample["gradient_s_inv"][0]
    u = sample["velocity_m_s"][0]
    _, r = execute("CFD_single_particle", x, flow, duration=0.0005, vessel=True)
    tensor = np.zeros((4, 4))
    tensor[:3, :3] = G
    tensor[:3, 3] = u - G @ x[0]
    exact = (expm(tensor * 0.0005) @ np.r_[x[0], 1.0])[:3]
    xx = r["positions_m"][:, 0, 0]
    weights = cpu.sampler.weights(xx, np.full(len(xx), cell))
    assert weights.min() > 0
    local = x[0].copy()
    for _ in range(50):
        local = local + 1e-5 * cpu.velocity(local[None], 0)[0]
    agreement = float(np.linalg.norm(xx[-1] - local))
    integration_error = float(np.linalg.norm(xx[-1] - exact))
    assert agreement < 5e-14 and integration_error < 5e-8
    rows["CFD_single_particle"].update(
        CPU_same_dt_agreement_m=agreement,
        Euler_integration_error_vs_exact_local_P1_m=integration_error,
        exact_reference_stays_in_same_tetra=True,
    )

    protocol = js(BASE / "ACOUSTIC_RELEASE_PROTOCOL.json")
    on = protocol["events"][0]
    off = dict(protocol["events"][1], time_s=0.0)
    parent = np.asarray(protocol["positions_m"])
    semantic = F5FlowSimulator(parent, p, response, basis, flow=flow, vessel=geometry)
    secondary_off = dict(on, secondary_enabled=False, mode="SECONDARY_OFF_PZT_ON")
    ss = evaluate(semantic, secondary_off)
    assert np.linalg.norm(ss["primary_N"]) > 0 and np.array_equal(
        ss["secondary_N"], np.zeros_like(ss["secondary_N"])
    )
    zz = evaluate(semantic, off)
    assert (
        np.array_equal(zz["total_N"], np.zeros_like(zz["total_N"]))
        and np.linalg.norm(zz["free_velocity_m_s"]) > 0
    )
    rows["OFF_semantics"] = dict(
        SECONDARY_OFF_PZT_ON_primary_nonzero=True,
        ULTRASOUND_OFF_FLOW_ON_acoustics_zero=True,
        flow_continues=True,
    )

    x = np.array([[0.006, 0.006 * np.tan(np.pi / 6), 0.0]])
    outside_field, r = execute(
        "OFF_beyond_acoustic_grid", x, flow, [off], duration=0.001, real_basis=True, vessel=True
    )
    assert r["active"].all() and r["summary"]["field_exit_count"] == 0
    assert np.linalg.norm(r["positions_m"][-1] - r["positions_m"][0]) > 1e-6
    rows["OFF_beyond_acoustic_grid"]["displacement_m"] = float(
        np.linalg.norm(r["positions_m"][-1] - r["positions_m"][0])
    )
    y = float(np.asarray(outside_field.field_upper)[1]) - 5e-8
    x = np.array([[y / np.tan(np.pi / 6), y, 0.0]])
    _, r = execute("ON_acoustic_domain_exit", x, flow, [on], duration=0.001, real_basis=True, vessel=True)
    assert r["summary"]["field_exit_count"] == 1 and r["boundary_events"][0]["label"] == 100
    rows["ON_acoustic_domain_exit"]["events"] = r["boundary_events"]
    planes = js(INPUT / "f4/PLANES.json")
    cap = planes["OUTLET_A"]
    normal = np.array(cap["normal"])
    center = np.array(cap["origin_m"])
    _, r = execute(
        "OFF_physical_outlet_exit",
        (center - 1e-7 * normal)[None],
        flow,
        [off],
        duration=0.001,
        real_basis=True,
        vessel=True,
    )
    assert r["summary"]["terminal_exit_count"] == 1 and r["boundary_events"][0]["label"] == 3
    rows["OFF_physical_outlet_exit"]["events"] = r["boundary_events"]

    base_event = dict(on, flow_enabled=False, mode="ACOUSTIC_ON_FLOW_OFF")
    f4 = F4VesselSimulator(parent, p, response, basis, vessel=geometry)
    r4 = f4.run(0.001, 1e-5, [base_event], 0.0001, 30)
    f5 = F5FlowSimulator(parent, p, response, basis, flow=flow, flow_enabled=False, vessel=geometry)
    r5 = f5.run(0.001, 1e-5, [base_event], 0.0001, 30)
    np.testing.assert_array_equal(r4["positions_m"], r5["positions_m"])
    rows["flow_disabled_F4_compatibility"] = dict(bitwise_positions=True, steps=100)
    write(out / "RESULT.json", dict(status="PASS", rows=rows, identity=identity(), warp=wp.__version__))
    print("F5 basic checks PASS", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    main(ap.parse_args().label)
