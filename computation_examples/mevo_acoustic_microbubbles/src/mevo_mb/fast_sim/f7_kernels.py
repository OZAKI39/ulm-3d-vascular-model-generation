"""Cached linear basis samples at actual positions, production response and force kernels.

Only the four complex controls vary inside a candidate evaluation. Spatial samples
are rebuilt at every observed state. No candidate moves particles or changes radii.
"""

import warp as wp
from .kernels import cmul


@wp.kernel
def cached_response(
    samples: wp.array2d(dtype=wp.vec2d),
    gradients: wp.array3d(dtype=wp.vec2d),
    controls: wp.array(dtype=wp.vec2d),
    knots: wp.array(dtype=wp.float64),
    values: wp.array3d(dtype=wp.vec2d),
    valid: wp.array(dtype=int),
    kv: wp.array2d(dtype=wp.vec2d),
    kval: wp.array(dtype=int),
    harmonics: wp.array2d(dtype=wp.vec2d),
    force: wp.array(dtype=wp.vec3d),
    pressure: wp.array(dtype=wp.vec2d),
    diag: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
):
    i = wp.tid()
    p = wp.vec2d(wp.float64(0.0))
    gx = p
    gy = p
    gz = p
    for j in range(4):
        p = p + cmul(samples[i, j], controls[j])
        gx = gx + cmul(gradients[i, 0, j], controls[j])
        gy = gy + cmul(gradients[i, 1, j], controls[j])
        gz = gz + cmul(gradients[i, 2, j], controls[j])
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
