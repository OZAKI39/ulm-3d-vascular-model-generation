"""Independent NumPy one/pair oracle, controls, graphs, real budget rejection."""

import argparse
from dataclasses import replace
import numpy as np
from mevo_mb.fast_sim.f10f_smallsteps import SmallStepsBatch, MODES
from mevo_mb.fast_sim.flow_gpu import AffineFlow
from mevo_mb.fast_sim.consistent_field import query_numpy
from mevo_mb.fast_sim.nearfield import pair_reference
from .common import BASE, HISTORY, ROOT, objects, identity, js, write


def pbd(x, radius, iterations, relaxation=0.5):
    x = x.copy()
    rounds = []
    for _ in range(iterations):
        dx = np.zeros_like(x)
        for i in range(len(x)):
            for j in range(len(x)):
                if i == j:
                    continue
                v = x[i] - x[j]
                d = np.linalg.norm(v)
                if d < 2 * radius:
                    dx[i] += 0.5 * (2 * radius - d) * v / d
        x += relaxation * dx
        rounds.append(x.copy())
    return x, rounds


def cpu_step(x, dt, p, r, b, ev, velocity):
    if np.any(np.asarray(ev["amplitude_Pa"])):
        field, _, _ = b.combine(ev["amplitude_Pa"], ev["phase_rad"])
        pressure, grad = query_numpy(field.pressure_Pa, b.origin, b.spacing, x)
        v = r.query(pressure)
        force = -0.5 * np.real(v[:, 0, None] * grad.conj())
        if ev.get("secondary_enabled", True):
            force += pair_reference(
                x, v, pressure, p, r.query(np.array([r.knots[0] + 0j]))[0, 0] / r.knots[0]
            )
    else:
        force = np.zeros_like(x)
    free = force / (6 * np.pi * p.viscosity_Pa_s * p.radius_m) + velocity
    old = free.copy()

    def phi(e):
        return 0.5 / e - 0.45 * np.log(e) - (3 / 56) * e * np.log(e)

    for _ in range(2):
        new = old.copy()
        for i in range(len(x)):
            mat = np.eye(3)
            rhs = free[i].copy()
            for j in range(len(x)):
                if i == j:
                    continue
                d = np.linalg.norm(x[j] - x[i])
                gap = (d - 2 * p.radius_m) / p.radius_m
                if gap < 0.1:
                    n = (x[j] - x[i]) / d
                    c = 0.5 * (phi(max(gap, 0.001)) - phi(0.1))
                    mat += c * np.outer(n, n)
                    rhs += c * np.dot(n, old[j]) * n
            new[i] = 0.5 * (old[i] + np.linalg.solve(mat, rhs))
        old = new
    if len(x) == 2 and np.linalg.norm(x[1] - x[0]) < 2.1 * p.radius_m:
        old += (free - old).mean(0)
    y, _ = pbd(x + dt * old, p.radius_m, p.contact_iterations)
    return y, force


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    p, b, r, g, f = objects()
    selection = js(HISTORY / "SELECTION.json")["branch"]
    cfg = js(ROOT / selection["relative_path"] / "CONFIG.json")
    center = np.array(cfg["initial_positions_m"])[0].mean(0)
    initial_event = cfg["events"][0]
    velocity = np.array([2e-4, 0, 0])
    flow = AffineFlow(velocity)
    rows = []
    for mode, (dt, iterations) in MODES.items():
        pp = replace(p, contact_iterations=iterations)
        for n in [1, 2]:
            x = np.repeat(center[None], n, axis=0)
            if n == 2:
                x[:, 0] += np.array([-1, 1]) * (p.radius_m - 5e-9)
            events = [
                initial_event | {"time_s": 0, "flow_enabled": True},
                initial_event | {"time_s": 2 * dt, "amplitude_Pa": [0] * 4, "flow_enabled": True},
                initial_event | {"time_s": 3 * dt, "flow_enabled": True},
            ]
            gpu = []
            for graphs in [False, True]:
                s = SmallStepsBatch(x, pp, r, b, vessel=g, flow=flow, numerical_mode=mode, diagnostic=False)
                s.configure_schedule(events, dt)
                s.control(events[0])
                s.graph(dt, 1)
                np.testing.assert_array_equal(s.x.numpy(), x)
                s.advance(dt, 8, graphs=graphs)
                gpu.append(s.x.numpy().copy())
                assert s.accepted_steps == 8 and s.physical_steps == 8
                assert [e["time_s"] for e in s.executed_controls] == [0, 2 * dt, 3 * dt]
                count = int(s.event_count.numpy()[0])
                et = s.event_time.numpy()[:count]
                assert np.all(np.abs(et / dt - np.round(et / dt)) < 1e-10)
            np.testing.assert_array_equal(*gpu)
            cpu = x.copy()
            for step in range(8):
                ev = events[0 if step < 2 else 1 if step < 3 else 2]
                cpu, force = cpu_step(cpu, dt, pp, r, b, ev, velocity)
            error = float(np.max(np.abs(cpu - gpu[-1])))
            force_error = float(np.max(np.abs(force - s.force.numpy())))
            assert error <= 1e-12 and force_error <= 1e-22, (error, force_error)
            rows.append(
                dict(
                    mode=mode,
                    particle_count=n,
                    CPU_GPU_position_max_m=error,
                    CPU_GPU_force_max_N=force_error,
                    graph_raw_positions_equal=True,
                    steps=8,
                    actual_end_s=8 * dt,
                    control_times_s=[0, 2 * dt, 3 * dt],
                    quality=s.quality(),
                )
            )
    # Fail the original single-pass budget even when final PBD residual converges.
    x = np.repeat(center[None], 2, axis=0)
    x[:, 0] += np.array([-1, 1]) * (p.radius_m - 0.4e-6)
    s = SmallStepsBatch(x, p, r, b, vessel=g, flow=flow, numerical_mode="B0", diagnostic=False)
    off = initial_event | {"time_s": 0, "amplitude_Pa": [0] * 4, "flow_enabled": True}
    s.configure_schedule([off], 1e-5)
    s.control(off)
    try:
        s.advance(1e-5, 5, graphs=True)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Real correction budget must reject")
    assert s.physical_steps == 1 and s.accepted_steps == 0
    assert s.quality()["max_contact_pass_correction_m"] > 178e-9
    np.testing.assert_array_equal(s.x.numpy(), x)
    first = s.first_rejection.copy()
    try:
        s.advance(1e-5, 5, graphs=True)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Rejected object cannot advance")
    assert s.physical_steps == 1 and s.first_rejection == first
    # A user can't pass an unaligned switch to this fixed-substep mode.
    try:
        s.configure_schedule([off, off | {"time_s": 1.1e-6}], 1e-5)
    except ValueError:
        pass
    else:
        raise AssertionError("Unaligned control must be explicitly rejected")
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            cases=rows,
            real_single_pass_rejection=first,
            post_failure_no_motion=True,
            nonaligned_controls_rejected=True,
            scope="SOFTWARE_TEST_ONLY; affine flow for independent CPU oracle; real CFD covered by deterministic scientific replay",
            training_ready=False,
        ),
    )
    print("PASS engineering checks", flush=True)


if __name__ == "__main__":
    main()
