"""Same F2 pair expression, stored once; fixed-ID reductions change only summation order."""

import warp as wp
from .kernels import cmul, root_of
from .f2_kernels import cdiv


@wp.kernel
def pair_buffer(
    x: wp.array(dtype=wp.vec3d),
    v: wp.array2d(dtype=wp.vec2d),
    pressure: wp.array(dtype=wp.vec2d),
    pairs: wp.array2d(dtype=wp.vec3d),
    n: int,
    radius: wp.float64,
    omega: wp.float64,
    rho: wp.float64,
    transition: wp.float64,
    soft: wp.float64,
    end: wp.float64,
    chi_zero: wp.vec2d,
    pmin: wp.float64,
    pmax: wp.float64,
    near_diag: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
):
    b, i, j = wp.tid()
    if j > i:
        ii = b * n + i
        jj = b * n + j
        delta = x[jj] - x[ii]
        d = wp.length(delta)
        if d <= wp.float64(0.0):
            wp.atomic_or(error, 0, 4)
        else:
            corr = wp.float64(0.0)
            for h in range(8):
                corr = corr + wp.float64((h + 1) * (h + 1)) * wp.dot(v[ii, h], v[jj, h])
            s = wp.clamp(
                (d / (wp.float64(2.0) * radius) - wp.float64(1.0)) / (end - wp.float64(1.0)),
                wp.float64(0.0),
                wp.float64(1.0),
            )
            w = wp.float64(1.0) - s * s * (wp.float64(3.0) - wp.float64(2.0) * s)
            if w > wp.float64(0.0):
                ci = chi_zero
                cj = chi_zero
                if wp.dot(pressure[ii], pressure[ii]) > wp.float64(0.0):
                    ci = cdiv(v[ii, 0], pressure[ii])
                if wp.dot(pressure[jj], pressure[jj]) > wp.float64(0.0):
                    cj = cdiv(v[jj, 0], pressure[jj])
                k = rho * omega * omega / (wp.float64(12.566370614359172) * d)
                a = k * ci
                c = k * cj
                det = wp.vec2d(wp.float64(1.0), wp.float64(0.0)) - cmul(a, c)
                determinant = wp.length(det)
                wp.atomic_min(near_diag, ii, 2, determinant)
                wp.atomic_min(near_diag, jj, 2, determinant)
                if determinant < wp.float64(1e-6):
                    wp.atomic_or(error, 0, 64)
                else:
                    vi = cdiv(v[ii, 0] - cmul(a, v[jj, 0]), det)
                    vj = cdiv(v[jj, 0] - cmul(c, v[ii, 0]), det)
                    pi = wp.length(pressure[ii] - k * vj)
                    pj = wp.length(pressure[jj] - k * vi)
                    wp.atomic_min(near_diag, ii, 0, pi)
                    wp.atomic_min(near_diag, jj, 0, pj)
                    wp.atomic_max(near_diag, ii, 1, pi)
                    wp.atomic_max(near_diag, jj, 1, pj)
                    if (
                        pi > pmax
                        or pj > pmax
                        or (pi > wp.float64(0.0) and pi < pmin)
                        or (pj > wp.float64(0.0) and pj < pmin)
                    ):
                        wp.atomic_or(error, 0, 128)
                    corr = corr + w * (wp.dot(vi, vj) - wp.dot(v[ii, 0], v[jj, 0]))
            q = wp.min(d / (transition * wp.float64(2.0) * radius), wp.float64(1.0))
            core = soft * wp.float64(2.0) * radius * (wp.float64(1.0) - q) * (wp.float64(1.0) - q)
            d2 = d * d + core * core
            f = rho * omega * omega / (wp.float64(25.132741228718345) * d2 * wp.sqrt(d2)) * corr * delta
            pairs[ii, j] = f
            pairs[jj, i] = -f


@wp.kernel
def pair_sum(pairs: wp.array2d(dtype=wp.vec3d), force: wp.array(dtype=wp.vec3d), n: int):
    i = wp.tid()
    value = force[i]
    for j in range(n):
        if j != i % n:
            value = value + pairs[i, j]
    force[i] = value


@wp.kernel
def component_balance_ordered(
    parent: wp.array(dtype=int),
    free: wp.array(dtype=wp.vec3d),
    velocity: wp.array(dtype=wp.vec3d),
    delta: wp.array(dtype=wp.vec3d),
    count: wp.array(dtype=int),
    n: int,
):
    i = wp.tid()
    if root_of(parent, i) == i:
        start = (i // n) * n
        total = wp.vec3d(wp.float64(0.0))
        size = int(0)
        for j in range(start, start + n):
            if root_of(parent, j) == i:
                total = total + (free[j] - velocity[j])
                size = size + 1
        delta[i] = total
        count[i] = size
