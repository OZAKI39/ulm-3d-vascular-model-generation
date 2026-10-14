"""Cubic pressure and its exact derivative; inherited strict APECSS interval lookup."""

import warp as wp
from .kernels import cmul


@wp.func
def weight(t: wp.float64, j: int):
    v = wp.float64(0.0)
    if j == 0:
        v = -wp.float64(0.5) * t + t * t - wp.float64(0.5) * t * t * t
    elif j == 1:
        v = wp.float64(1.0) - wp.float64(2.5) * t * t + wp.float64(1.5) * t * t * t
    elif j == 2:
        v = wp.float64(0.5) * t + wp.float64(2.0) * t * t - wp.float64(1.5) * t * t * t
    else:
        v = -wp.float64(0.5) * t * t + wp.float64(0.5) * t * t * t
    return v


@wp.func
def derivative(t: wp.float64, j: int):
    v = wp.float64(0.0)
    if j == 0:
        v = -wp.float64(0.5) + wp.float64(2.0) * t - wp.float64(1.5) * t * t
    elif j == 1:
        v = -wp.float64(5.0) * t + wp.float64(4.5) * t * t
    elif j == 2:
        v = wp.float64(0.5) + wp.float64(4.0) * t - wp.float64(4.5) * t * t
    else:
        v = -t + wp.float64(1.5) * t * t
    return v


@wp.kernel
def cubic_response(
    x: wp.array(dtype=wp.vec3d),
    grid: wp.array2d(dtype=wp.vec2d),
    origin: wp.vec3d,
    spacing: wp.vec3d,
    nx: int,
    ny: int,
    nz: int,
    knots: wp.array(dtype=wp.float64),
    values: wp.array3d(dtype=wp.vec2d),
    valid: wp.array(dtype=int),
    kv: wp.array2d(dtype=wp.vec2d),
    kval: wp.array(dtype=int),
    harmonics: wp.array2d(dtype=wp.vec2d),
    force: wp.array(dtype=wp.vec3d),
    pressure: wp.array(dtype=wp.vec2d),
    gradients: wp.array2d(dtype=wp.vec2d),
    diag: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
):
    i = wp.tid()
    if error[0] != 0:
        return
    q = wp.cw_div(x[i] - origin, spacing)
    if (
        not wp.isfinite(q[0])
        or not wp.isfinite(q[1])
        or not wp.isfinite(q[2])
        or q[0] < wp.float64(1.0)
        or q[1] < wp.float64(1.0)
        or q[2] < wp.float64(1.0)
        or q[0] >= wp.float64(nx - 2)
        or q[1] >= wp.float64(ny - 2)
        or q[2] >= wp.float64(nz - 2)
    ):
        wp.atomic_or(error, 0, 1)
        return
    ix = int(wp.floor(q[0]))
    iy = int(wp.floor(q[1]))
    iz = int(wp.floor(q[2]))
    t = q - wp.vec3d(wp.float64(ix), wp.float64(iy), wp.float64(iz))
    p = wp.vec2d(wp.float64(0.0))
    gx = p
    gy = p
    gz = p
    for a in range(4):
        for b in range(4):
            for c in range(4):
                wa = weight(t[0], a)
                wb = weight(t[1], b)
                wc = weight(t[2], c)
                idx = ((ix + a - 1) * ny + iy + b - 1) * nz + iz + c - 1
                v = grid[idx, 0]
                p = p + wa * wb * wc * v
                gx = gx + derivative(t[0], a) * wb * wc / spacing[0] * v
                gy = gy + wa * derivative(t[1], b) * wc / spacing[1] * v
                gz = gz + wa * wb * derivative(t[2], c) / spacing[2] * v
    gradients[i, 0] = gx
    gradients[i, 1] = gy
    gradients[i, 2] = gz
    pressure[i] = p
    amp = wp.length(p)
    diag[i, 0] = wp.min(diag[i, 0], amp)
    diag[i, 1] = wp.max(diag[i, 1], amp)
    cell = int(-1)
    knot = int(-1)
    for k in range(knots.shape[0]):
        if amp == knots[k]:
            knot = k
        if k < knots.shape[0] - 1:
            if amp > knots[k] and amp < knots[k + 1]:
                cell = k
    ok = int(1)
    if amp > wp.float64(0.0):
        if knot >= 0:
            ok = kval[knot]
        elif cell >= 0:
            ok = valid[cell]
        else:
            ok = 0
    if ok == 0:
        wp.atomic_or(error, 0, 2)
    theta = wp.atan2(p[1], p[0])
    for h in range(8):
        v = wp.vec2d(wp.float64(0.0))
        if amp > wp.float64(0.0) and ok == 1:
            if knot >= 0:
                v = kv[knot, h]
            else:
                fraction = (amp - knots[cell]) / (knots[cell + 1] - knots[cell])
                v = (wp.float64(1.0) - fraction) * values[cell, 0, h] + fraction * values[cell, 1, h]
            angle = wp.float64(h + 1) * theta
            v = cmul(v, wp.vec2d(wp.cos(angle), wp.sin(angle)))
        harmonics[i, h] = v
    v = harmonics[i, 0]
    force[i] = -wp.float64(0.5) * wp.vec3d(wp.dot(v, gx), wp.dot(v, gy), wp.dot(v, gz))
