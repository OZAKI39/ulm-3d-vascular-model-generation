"""Derivative of the same complex trilinear pressure; no field extrapolation.

Cell faces use the cell on the positive side, except the outer maximum uses the
last cell. Derivatives are one-sided there, and generally jump across cell faces.
The stored k-Wave samples have NOT gained spatial resolution or convergence.
"""

import numpy as np
import warp as wp


def query_numpy(pressure, origin, spacing, positions):
    p = np.asarray(pressure)
    x = np.atleast_2d(positions)
    q = (x - origin) / spacing
    shape = np.array(p.shape)
    if not np.isfinite(q).all() or np.any(x < origin) or np.any(x > origin + spacing * (shape - 1)):
        raise ValueError("Outside stored field; no extrapolation")
    q = np.clip(q, 0, shape - 1)  # only quotient roundoff AFTER strict physical-domain check
    cell = np.minimum(np.floor(q).astype(int), shape - 2)
    w = q - cell
    value = np.zeros(len(x), complex)
    grad = np.zeros((len(x), 3), complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                corner = np.array([a, b, c])
                weights = np.where(corner, w, 1 - w)
                v = p[tuple((cell + corner).T)]
                value += v * np.prod(weights, axis=1)
                for axis in range(3):
                    other = [d for d in range(3) if d != axis]
                    grad[:, axis] += (
                        v * (2 * corner[axis] - 1) / spacing[axis] * np.prod(weights[:, other], axis=1)
                    )
    return value, grad


@wp.kernel
def consistent_primary(
    x: wp.array(dtype=wp.vec3d),
    grid: wp.array2d(dtype=wp.vec2d),
    origin: wp.vec3d,
    spacing: wp.vec3d,
    upper: wp.vec3d,
    nx: int,
    ny: int,
    nz: int,
    active: wp.array(dtype=int),
    harmonics: wp.array2d(dtype=wp.vec2d),
    force: wp.array(dtype=wp.vec3d),
    pressure: wp.array(dtype=wp.vec2d),
    gradients: wp.array2d(dtype=wp.vec2d),
    error: wp.array(dtype=int),
):
    i = wp.tid()
    force[i] = wp.vec3d(wp.float64(0.0))
    if active[i] == 0:
        return
    q = wp.cw_div(x[i] - origin, spacing)
    if (
        not wp.isfinite(q[0])
        or not wp.isfinite(q[1])
        or not wp.isfinite(q[2])
        or x[i][0] < origin[0]
        or x[i][1] < origin[1]
        or x[i][2] < origin[2]
        or x[i][0] > upper[0]
        or x[i][1] > upper[1]
        or x[i][2] > upper[2]
    ):
        wp.atomic_or(error, 0, 1)
        return
    q = wp.vec3d(
        wp.clamp(q[0], wp.float64(0.0), wp.float64(nx - 1)),
        wp.clamp(q[1], wp.float64(0.0), wp.float64(ny - 1)),
        wp.clamp(q[2], wp.float64(0.0), wp.float64(nz - 1)),
    )
    ix = wp.min(int(wp.floor(q[0])), nx - 2)
    iy = wp.min(int(wp.floor(q[1])), ny - 2)
    iz = wp.min(int(wp.floor(q[2])), nz - 2)
    w = q - wp.vec3d(wp.float64(ix), wp.float64(iy), wp.float64(iz))
    p = wp.vec2d(wp.float64(0.0))
    gx = p
    gy = p
    gz = p
    for a in range(2):
        for b in range(2):
            for c in range(2):
                wa = wp.float64(1.0) - w[0]
                wb = wp.float64(1.0) - w[1]
                wc = wp.float64(1.0) - w[2]
                if a == 1:
                    wa = w[0]
                if b == 1:
                    wb = w[1]
                if c == 1:
                    wc = w[2]
                idx = ((ix + a) * ny + iy + b) * nz + iz + c
                v = grid[idx, 0]
                p = p + wa * wb * wc * v
                gx = gx + wp.float64(2 * a - 1) * wb * wc / spacing[0] * v
                gy = gy + wp.float64(2 * b - 1) * wa * wc / spacing[1] * v
                gz = gz + wp.float64(2 * c - 1) * wa * wb / spacing[2] * v
    pressure[i] = p
    gradients[i, 0] = gx
    gradients[i, 1] = gy
    gradients[i, 2] = gz
    v1 = harmonics[i, 0]
    force[i] = -wp.float64(0.5) * wp.vec3d(wp.dot(v1, gx), wp.dot(v1, gy), wp.dot(v1, gz))
