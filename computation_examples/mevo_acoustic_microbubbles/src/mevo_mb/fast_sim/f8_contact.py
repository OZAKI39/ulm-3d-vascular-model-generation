"""Bounded frozen-mobility normal contact prototype; never a production default.

Dimensionless mobility H=M/mu0, B=H J^T. Projected coordinate updates solve
0 <= lambda perpendicular to gap + J B lambda >= 0, with lambda in metres
(physical contact impulse=lambda/mu0). No full constraint matrix or QP solver.
The small-state CPU assembly cost is explicitly included, not hidden as GPU work.
"""

import time
import numpy as np
import warp as wp
from .rpy import matrix_numpy


@wp.kernel
def projected_contacts(
    ij: wp.array2d(dtype=int),
    normals: wp.array(dtype=wp.vec3d),
    gap: wp.array(dtype=wp.float64),
    columns: wp.array2d(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.float64),
    multipliers: wp.array(dtype=wp.float64),
    displacement: wp.array(dtype=wp.float64),
    trace: wp.array2d(dtype=wp.float64),
    passes: int,
    dimensions: int,
    count: int,
):
    for c in range(count):
        multipliers[c] = wp.float64(0.0)
    for k in range(dimensions):
        displacement[k] = wp.float64(0.0)
    for it in range(passes):
        largest_update = wp.float64(0.0)
        for c in range(count):
            i = ij[c, 0]
            j = ij[c, 1]
            value = gap[c]
            for a in range(3):
                value += normals[c][a] * (displacement[3 * i + a] - displacement[3 * j + a])
            old = multipliers[c]
            new = wp.max(wp.float64(0.0), old - value / diagonal[c])
            change = new - old
            multipliers[c] = new
            for k in range(dimensions):
                delta = columns[c, k] * change
                displacement[k] += delta
                largest_update = wp.max(largest_update, wp.abs(delta))
        residual = wp.float64(0.0)
        for c in range(count):
            value = gap[c]
            for a in range(3):
                value += normals[c][a] * (displacement[3 * ij[c, 0] + a] - displacement[3 * ij[c, 1] + a])
            # Projected-coordinate fixed point residual has units of metres.
            residual = wp.max(
                residual,
                wp.abs(multipliers[c] - wp.max(wp.float64(0.0), multipliers[c] - value / diagonal[c])),
            )
        trace[it, 0] = residual
        trace[it, 1] = largest_update


def prepare(reference, predicted, radius, viscosity, device="cuda:0"):
    begin = time.perf_counter()
    x = np.asarray(reference)
    free = np.asarray(predicted)
    n = len(x)
    mu0 = 1 / (6 * np.pi * viscosity * radius)
    H = matrix_numpy(x, radius, viscosity) / mu0
    i, j = np.triu_indices(n, 1)
    ij = np.column_stack([i, j]).astype(np.int32)
    delta = free[i] - free[j]
    distance = np.linalg.norm(delta, axis=1)
    if np.any(distance <= 0):
        raise ValueError("Coincident centers cannot define a normal")
    normal = delta / distance[:, None]
    J = np.zeros((len(i), 3 * n))
    for c in range(len(i)):
        J[c, 3 * i[c] : 3 * i[c] + 3] = normal[c]
        J[c, 3 * j[c] : 3 * j[c] + 3] = -normal[c]
    B = H @ J.T
    diag = np.einsum("ck,kc->c", J, B)
    if np.any(diag <= 0):
        raise ValueError("Nonpositive contact self response")
    wp.load_module(module=projected_contacts.module, device=device)
    args = [
        wp.array(ij, dtype=int, device=device),
        wp.array(normal, dtype=wp.vec3d, device=device),
        wp.array(distance - 2 * radius, dtype=wp.float64, device=device),
        wp.array(B.T.copy(), dtype=wp.float64, device=device),
        wp.array(diag, dtype=wp.float64, device=device),
        wp.zeros(len(i), dtype=wp.float64, device=device),
        wp.zeros(3 * n, dtype=wp.float64, device=device),
    ]
    wp.synchronize_device(device)
    return dict(
        H=H,
        J=J,
        B=B,
        gap=distance - 2 * radius,
        mu0=mu0,
        args=args,
        n=n,
        ij=ij,
        setup_wall_s=time.perf_counter() - begin,
        device=device,
    )


def solve(prepared, passes=16):
    d = prepared
    trace = wp.zeros((passes, 2), dtype=wp.float64, device=d["device"])
    start, end = wp.Event(enable_timing=True), wp.Event(enable_timing=True)
    begin = time.perf_counter()
    wp.record_event(start)
    wp.launch(
        projected_contacts,
        dim=1,
        inputs=[*d["args"], trace, passes, 3 * d["n"], len(d["ij"])],
        device=d["device"],
    )
    wp.record_event(end)
    wp.synchronize_event(end)
    return dict(
        displacement_m=d["args"][6].numpy().reshape(-1, 3),
        lambda_m=d["args"][5].numpy(),
        trace=trace.numpy(),
        gpu_launch_span_s=wp.get_event_elapsed_time(start, end) / 1000,
        host_wall_s=time.perf_counter() - begin,
    )
