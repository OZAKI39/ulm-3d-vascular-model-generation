"""Bounded fixed-state GPU microbenchmarks, not trajectory throughput."""

import argparse
from time import perf_counter
import numpy as np
import warp as wp
from mevo_mb.fast_sim import kernels as k
from mevo_mb.fast_sim import vessel_active_kernels as masked
from mevo_mb.fast_sim import f3_kernels as ordered
from mevo_mb.fast_sim.f10g_contact import resistance, jacobi, balance, pbd
from .common import BASE, CASES, state, identity, write


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    wp.init()
    device = "cuda:0"
    wp.load_module(module=k, device=device)
    wp.load_module(module=masked, device=device)
    wp.load_module(module=ordered, device=device)
    rows = []
    for case in CASES:
        _, t, _ = state(case)
        x0 = np.array(t["initial_positions_m"])
        u = np.array(t["free_velocity_m_s"])
        r = 1.78e-6
        op = resistance(x0, r)
        x = wp.array(x0, dtype=wp.vec3d, device=device)
        free = wp.array(u, dtype=wp.vec3d, device=device)
        active = wp.ones(32, dtype=int, device=device)
        v1 = wp.zeros_like(x)
        v2 = wp.zeros_like(x)
        p1 = wp.zeros_like(x)
        p2 = wp.zeros_like(x)
        parent = wp.zeros(32, dtype=int, device=device)
        count = wp.zeros(32, dtype=int, device=device)
        shift = wp.zeros_like(x)
        diag = wp.zeros((32, 11), dtype=wp.float64, device=device)
        err = wp.zeros(1, dtype=int, device=device)

        def local(n):
            velocity = free
            for it in range(n):
                nxt = v1 if it % 2 == 0 else v2
                wp.launch(
                    masked.local_mobility,
                    32,
                    [x, free, velocity, nxt, 32, r, 0.1, 0.001, active],
                    device=device,
                )
                velocity = nxt
            return velocity

        for n in [2, 4, 8]:
            for scope in ["lubrication_only", "lubrication_balance_prediction_PBD12"]:
                start = perf_counter()
                with wp.ScopedCapture(device=device) as capture:
                    vel = local(n)
                    if scope != "lubrication_only":
                        wp.launch(k.component_init, 32, [parent, shift, count], device=device)
                        wp.launch(
                            masked.component_links,
                            (1, 32, 32),
                            [x, parent, 32, 2.1 * r, active],
                            device=device,
                        )
                        wp.launch(
                            ordered.component_balance_ordered,
                            32,
                            [parent, free, vel, shift, count, 32],
                            device=device,
                        )
                        wp.launch(
                            k.predict, 32, [x, vel, p1, 1e-5, diag, parent, shift, count], device=device
                        )
                        for it in range(12):
                            src, dst = (p1, p2) if it % 2 == 0 else (p2, p1)
                            wp.launch(
                                masked.contact, 32, [src, dst, 32, r, 0.5, diag, err, active], device=device
                            )
                capture_s = perf_counter() - start
                wp.capture_launch(capture.graph)
                wp.synchronize_device(device)
                wall = []
                gpu = []
                for _ in range(25):
                    aev = wp.Event(enable_timing=True)
                    bev = wp.Event(enable_timing=True)
                    start = perf_counter()
                    wp.record_event(aev)
                    wp.capture_launch(capture.graph)
                    wp.record_event(bev)
                    wp.synchronize_event(bev)
                    wall.append(perf_counter() - start)
                    gpu.append(wp.get_event_elapsed_time(aev, bev) / 1000)
                start = perf_counter()
                actual = vel.numpy()
                read = perf_counter() - start
                error = float(np.max(np.abs(actual - jacobi(op, u, n))))
                assert error < 1e-12
                position_error = None
                if scope != "lubrication_only":
                    predicted = x0 + 1e-5 * balance(op, u, jacobi(op, u, n))
                    expected, _ = pbd(predicted, r)
                    position_error = float(np.max(np.abs(p1.numpy() - expected)))
                    assert position_error < 1e-15
                rows.append(
                    dict(
                        case=case,
                        iterations=n,
                        scope=scope,
                        capture_s=capture_s,
                        gpu_median_s=float(np.median(gpu)),
                        host_launch_sync_median_s=float(np.median(wall)),
                        gpu_samples_s=gpu,
                        host_samples_s=wall,
                        readback_s=read,
                        velocity_error_m_s=error,
                        position_error_m=position_error,
                        kernel_family="F4 masked active + F3 ordered component balance, identical current production arithmetic",
                        replays=25,
                        physical_steps_advanced=0,
                    )
                )
    write(
        out / "RESULTS.json",
        dict(
            rows=rows,
            interpretation="Original kernels on fixed real inputs; excludes force assembly, geometry, acceptance, morphology and output. No production throughput claim.",
            training_ready=False,
        ),
    )
    print("Fixed-state GPU kernel benchmark PASS", len(rows))


if __name__ == "__main__":
    main()
