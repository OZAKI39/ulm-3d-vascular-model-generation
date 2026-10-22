"""F2 extensions only; frozen F1 kernels are imported, never modified."""

import warp as wp
from .kernels import cmul, lub_phi


@wp.func
def cdiv(a: wp.vec2d, b: wp.vec2d):
    return wp.vec2d(a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1]) / wp.dot(b, b)


@wp.kernel
def pair_force(
    x: wp.array(dtype=wp.vec3d),
    v: wp.array2d(dtype=wp.vec2d),
    pressure: wp.array(dtype=wp.vec2d),
    force: wp.array(dtype=wp.vec3d),
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
            wp.atomic_add(force, ii, f)
            wp.atomic_add(force, jj, -f)


@wp.kernel
def mobility(
    x: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.mat33d),
    n: int,
    radius: wp.float64,
    activation: wp.float64,
    floor: wp.float64,
):
    i = wp.tid()
    start = (i // n) * n
    mat = wp.identity(n=3, dtype=wp.float64)
    for j in range(start, start + n):
        if i != j:
            delta = x[j] - x[i]
            d = wp.length(delta)
            e = (d - wp.float64(2.0) * radius) / radius
            if e < activation and d > wp.float64(0.0):
                normal = delta / d
                c = wp.float64(0.5) * (lub_phi(wp.max(e, floor)) - lub_phi(activation))
                mat = mat + c * wp.outer(normal, normal)
    out[i] = wp.inverse(mat)


@wp.kernel
def contact(
    x: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
    n: int,
    radius: wp.float64,
    relax: wp.float64,
    diag: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
    m: wp.array(dtype=wp.mat33d),
):
    i = wp.tid()
    start = (i // n) * n
    dx = wp.vec3d(wp.float64(0.0))
    penetration = wp.float64(0.0)
    for j in range(start, start + n):
        if j != i:
            delta = x[i] - x[j]
            d = wp.length(delta)
            overlap = wp.float64(2.0) * radius - d
            if overlap > wp.float64(0.0):
                penetration = wp.max(penetration, overlap)
                if d > wp.float64(0.0):
                    normal = delta / d
                    h = wp.float64(0.5) * (m[i] + m[j])
                    response = h * normal
                    metric = wp.dot(normal, response)
                    if metric > wp.float64(0.0):
                        dx = dx + wp.float64(0.5) * overlap * response / metric
                    else:
                        wp.atomic_or(error, 0, 256)
                else:
                    wp.atomic_or(error, 0, 4)
    dx = relax * dx
    out[i] = x[i] + dx
    diag[i, 2] = wp.max(diag[i, 2], penetration)
    diag[i, 4] = diag[i, 4] + wp.length(dx)
    diag[i, 6] = wp.max(diag[i, 6], wp.length(dx))


@wp.kernel
def free_path(v: wp.array(dtype=wp.vec3d), path: wp.array(dtype=wp.float64), dt: wp.float64):
    i = wp.tid()
    path[i] = path[i] + dt * wp.length(v[i])
