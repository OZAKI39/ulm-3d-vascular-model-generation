"""F4 active-ID masks around unchanged F1/F2/F3 arithmetic.
Generated once from archived functions; disabled vessel uses original kernels.
Inactive IDs have zero response and no pair/contact interactions.
"""

import warp as wp
from .kernels import cmul, lub_phi, root_of
from .f2_kernels import cdiv


@wp.kernel
def field_response(
    x: wp.array(dtype=wp.vec3d),
    grid: wp.array2d(dtype=wp.vec2d),
    origin: wp.vec3d,
    spacing: wp.vec3d,
    nx: int,
    ny: int,
    nz: int,
    uniform: int,
    uniform_p: wp.array(dtype=wp.vec2d),
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
    active: wp.array(dtype=int),
):
    i = wp.tid()
    if active[i] == 0:
        force[i] = wp.vec3d(wp.float64(0.0))
        pressure[i] = wp.vec2d(wp.float64(0.0))
        for h in range(8):
            harmonics[i, h] = wp.vec2d(wp.float64(0.0))
        return
    p = wp.vec2d(wp.float64(0.0))
    gx = p
    gy = p
    gz = p
    if uniform == 1:
        p = uniform_p[0]
    else:
        q = wp.cw_div(x[i] - origin, spacing)
        if (
            q[0] < wp.float64(0.0)
            or q[1] < wp.float64(0.0)
            or q[2] < wp.float64(0.0)
            or (q[0] > wp.float64(nx - 1))
            or (q[1] > wp.float64(ny - 1))
            or (q[2] > wp.float64(nz - 1))
        ):
            wp.atomic_or(error, 0, 1)
        else:
            ix = wp.min(int(wp.floor(q[0])), nx - 2)
            iy = wp.min(int(wp.floor(q[1])), ny - 2)
            iz = wp.min(int(wp.floor(q[2])), nz - 2)
            w = q - wp.vec3d(wp.float64(ix), wp.float64(iy), wp.float64(iz))
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
                        weight = wa * wb * wc
                        idx = ((ix + a) * ny + iy + b) * nz + iz + c
                        p = p + weight * grid[idx, 0]
                        gx = gx + weight * grid[idx, 1]
                        gy = gy + weight * grid[idx, 2]
                        gz = gz + weight * grid[idx, 3]
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


@wp.kernel
def local_mobility(
    x: wp.array(dtype=wp.vec3d),
    free: wp.array(dtype=wp.vec3d),
    old: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
    n: int,
    radius: wp.float64,
    activation: wp.float64,
    floor: wp.float64,
    active: wp.array(dtype=int),
):
    i = wp.tid()
    if active[i] == 0:
        out[i] = wp.vec3d(wp.float64(0.0))
        return
    start = i // n * n
    mat = wp.identity(n=3, dtype=wp.float64)
    rhs = free[i]
    for j in range(start, start + n):
        if j != i and active[j] == 1:
            delta = x[j] - x[i]
            d = wp.length(delta)
            e = (d - wp.float64(2.0) * radius) / radius
            if e < activation and d > wp.float64(0.0):
                normal = delta / d
                c = wp.float64(0.5) * (lub_phi(wp.max(e, floor)) - lub_phi(activation))
                mat = mat + c * wp.outer(normal, normal)
                rhs = rhs + c * wp.dot(normal, old[j]) * normal
    out[i] = wp.float64(0.5) * (old[i] + wp.inverse(mat) * rhs)


@wp.kernel
def component_links(
    x: wp.array(dtype=wp.vec3d),
    parent: wp.array(dtype=int),
    n: int,
    reach: wp.float64,
    active: wp.array(dtype=int),
):
    b, i, j = wp.tid()
    if j > i:
        ii = b * n + i
        jj = b * n + j
        if active[ii] == 0 or active[jj] == 0:
            return
        if wp.length(x[ii] - x[jj]) < reach:
            ri = root_of(parent, ii)
            rj = root_of(parent, jj)
            while ri != rj:
                high = wp.max(ri, rj)
                low = wp.min(ri, rj)
                previous = wp.atomic_cas(parent, high, high, low)
                if previous == high:
                    break
                ri = root_of(parent, ii)
                rj = root_of(parent, jj)


@wp.kernel
def contact(
    x: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
    n: int,
    radius: wp.float64,
    relaxation: wp.float64,
    diag: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
    active: wp.array(dtype=int),
):
    i = wp.tid()
    if active[i] == 0:
        out[i] = x[i]
        return
    start = i // n * n
    dx = wp.vec3d(wp.float64(0.0))
    penetration = wp.float64(0.0)
    for j in range(start, start + n):
        if i != j and active[j] == 1:
            delta = x[i] - x[j]
            d = wp.length(delta)
            overlap = wp.float64(2.0) * radius - d
            if overlap > wp.float64(0.0):
                penetration = wp.max(penetration, overlap)
                if d > wp.float64(0.0):
                    dx = dx + wp.float64(0.5) * overlap / d * delta
                else:
                    wp.atomic_or(error, 0, 4)
    dx = relaxation * dx
    out[i] = x[i] + dx
    diag[i, 2] = wp.max(diag[i, 2], penetration)
    diag[i, 4] = diag[i, 4] + wp.length(dx)
    diag[i, 6] = wp.max(diag[i, 6], wp.length(dx))


@wp.kernel
def accept(
    x: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
    n: int,
    radius: wp.float64,
    diag: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
    adjacency: wp.array2d(dtype=int),
    event_count: wp.array(dtype=int),
    event_ids: wp.array2d(dtype=int),
    event_time: wp.array(dtype=wp.float64),
    step_count: wp.array(dtype=int),
    time_origin: wp.array(dtype=wp.float64),
    dt: wp.float64,
    active: wp.array(dtype=int),
):
    if error[0] != 0:
        return
    i = wp.tid()
    start = i // n * n
    penetration = wp.float64(0.0)
    for j in range(start, start + n):
        if i != j:
            gap = (
                wp.length(x[i] - x[j]) - wp.float64(2.0) * radius
                if active[i] == 1 and active[j] == 1
                else wp.float64(1.0)
            )
            penetration = wp.max(penetration, -gap)
            if j > i:
                linked = int(0)
                if gap <= wp.float64(0.1) * radius:
                    linked = 1
                if linked != adjacency[i, j - start]:
                    idx = wp.atomic_add(event_count, 0, 1)
                    if idx < event_time.shape[0]:
                        event_ids[idx, 0] = i // n
                        event_ids[idx, 1] = i - start
                        event_ids[idx, 2] = j - start
                        event_ids[idx, 3] = linked
                        event_time[idx] = time_origin[0] + wp.float64(step_count[0] + 1) * dt
                    else:
                        wp.atomic_or(error, 0, 16)
                    adjacency[i, j - start] = linked
    if not wp.isfinite(x[i][0]) or not wp.isfinite(x[i][1]) or (not wp.isfinite(x[i][2])):
        wp.atomic_or(error, 0, 8)
    diag[i, 3] = wp.max(diag[i, 3], penetration)
    out[i] = x[i]


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
    active: wp.array(dtype=int),
):
    b, i, j = wp.tid()
    if j > i:
        ii = b * n + i
        jj = b * n + j
        if active[ii] == 0 or active[jj] == 0:
            pairs[ii, j] = wp.vec3d(wp.float64(0.0))
            pairs[jj, i] = wp.vec3d(wp.float64(0.0))
            return
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
                if determinant < wp.float64(1e-06):
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
    active: wp.array(dtype=int),
):
    b, i, j = wp.tid()
    if j > i:
        ii = b * n + i
        jj = b * n + j
        if active[ii] == 0 or active[jj] == 0:
            return
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
                if determinant < wp.float64(1e-06):
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
