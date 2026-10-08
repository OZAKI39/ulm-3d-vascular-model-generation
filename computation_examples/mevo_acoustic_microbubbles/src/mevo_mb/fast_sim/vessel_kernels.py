"""GPU static-BVH broad phase and local binary64 finite-triangle wall constraints.

No signed-distance assumption is made on the open WALL subset. Sign comes from
the closed labelled lumen. Exact-input rational checks live in the validation,
not in the physical step. Broad-phase padding is not a contact tolerance.
"""

import warp as wp


@wp.func
def f32(x: wp.vec3d):
    return wp.vec3(float(x[0]), float(x[1]), float(x[2]))


@wp.func
def f64(x: wp.vec3):
    return wp.vec3d(wp.float64(x[0]), wp.float64(x[1]), wp.float64(x[2]))


@wp.func
def closest_segment(p: wp.vec3d, a: wp.vec3d, b: wp.vec3d):
    e = b - a
    t = wp.clamp(
        wp.dot(p - a, e) / wp.max(wp.dot(e, e), wp.float64(1e-100)), wp.float64(0.0), wp.float64(1.0)
    )
    return a + t * e


@wp.func
def closest_triangle(p: wp.vec3d, a: wp.vec3d, b: wp.vec3d, c: wp.vec3d):
    e = b - a
    f = c - a
    v = p - a
    aa = wp.dot(e, e)
    ab = wp.dot(e, f)
    bb = wp.dot(f, f)
    den = aa * bb - ab * ab
    if den > wp.float64(0.0):
        s = (wp.dot(v, e) * bb - wp.dot(v, f) * ab) / den
        t = (wp.dot(v, f) * aa - wp.dot(v, e) * ab) / den
        if s >= wp.float64(0.0) and t >= wp.float64(0.0) and s + t <= wp.float64(1.0):
            return a + s * e + t * f
    q = closest_segment(p, a, b)
    r = closest_segment(p, b, c)
    if wp.length_sq(p - r) < wp.length_sq(p - q):
        q = r
    r = closest_segment(p, c, a)
    if wp.length_sq(p - r) < wp.length_sq(p - q):
        q = r
    return q


@wp.func
def segment_segment_distance2(p: wp.vec3d, q: wp.vec3d, a: wp.vec3d, b: wp.vec3d):
    best = wp.min(wp.length_sq(p - closest_segment(p, a, b)), wp.length_sq(q - closest_segment(q, a, b)))
    best = wp.min(best, wp.length_sq(a - closest_segment(a, p, q)))
    best = wp.min(best, wp.length_sq(b - closest_segment(b, p, q)))
    u = q - p
    v = b - a
    w = p - a
    aa = wp.dot(u, u)
    bb = wp.dot(u, v)
    cc = wp.dot(v, v)
    dd = wp.dot(u, w)
    ee = wp.dot(v, w)
    den = aa * cc - bb * bb
    if den > wp.float64(0.0):
        s = (bb * ee - cc * dd) / den
        t = (aa * ee - bb * dd) / den
        if s >= wp.float64(0.0) and s <= wp.float64(1.0) and t >= wp.float64(0.0) and t <= wp.float64(1.0):
            best = wp.min(best, wp.length_sq(p + s * u - a - t * v))
    return best


@wp.func
def swept_distance2(p: wp.vec3d, q: wp.vec3d, a: wp.vec3d, b: wp.vec3d, c: wp.vec3d):
    best = wp.min(
        wp.length_sq(p - closest_triangle(p, a, b, c)), wp.length_sq(q - closest_triangle(q, a, b, c))
    )
    normal = wp.cross(b - a, c - a)
    den = wp.dot(normal, q - p)
    if wp.abs(den) > wp.float64(1e-100):
        t = wp.dot(normal, a - p) / den
        if t >= wp.float64(0.0) and t <= wp.float64(1.0):
            hit = p + t * (q - p)
            # Same finite-triangle barycentric predicate as the closest-point function.
            closest = closest_triangle(hit, a, b, c)
            best = wp.min(best, wp.length_sq(hit - closest))
    best = wp.min(best, segment_segment_distance2(p, q, a, b))
    best = wp.min(best, segment_segment_distance2(p, q, b, c))
    best = wp.min(best, segment_segment_distance2(p, q, c, a))
    return best


@wp.struct
class WallHit:
    distance: wp.float64
    normal: wp.vec3d
    closest: wp.vec3d
    face: int
    valid: int
    closed_inside: int
    triangle_tests: int


@wp.func
def wall_query(
    closed: wp.uint64,
    wall: wp.uint64,
    vertices: wp.array(dtype=wp.vec3d),
    faces: wp.array2d(dtype=int),
    tags: wp.array(dtype=int),
    wall_faces: wp.array2d(dtype=int),
    origin: wp.vec3d,
    padding: wp.float64,
    p: wp.vec3d,
):
    out = WallHit()
    out.valid = 0
    out.triangle_tests = 0
    point = f32(p - origin)
    coarse = wp.mesh_query_point_no_sign(wall, point, 1.0)
    side = wp.mesh_query_point(closed, point, 1.0)
    if coarse.result and side.result:
        c32 = wp.mesh_eval_position(wall, coarse.face, coarse.u, coarse.v)
        bound = wp.float64(wp.length(point - c32)) + padding
        query = wp.mesh_query_aabb(wall, f32(p - origin - wp.vec3d(bound)), f32(p - origin + wp.vec3d(bound)))
        face = int(0)
        best = wp.float64(1.0)
        q = wp.vec3d(wp.float64(0.0))
        winner = int(-1)
        while wp.mesh_query_aabb_next(query, face):
            a = vertices[wall_faces[face, 0]]
            b = vertices[wall_faces[face, 1]]
            c = vertices[wall_faces[face, 2]]
            near = closest_triangle(p, a, b, c)
            dd = wp.length_sq(p - near)
            if dd < best:
                best = dd
                q = near
                winner = face
            out.triangle_tests = out.triangle_tests + 1
        if winner >= 0:
            d = wp.sqrt(best)
            sign = wp.float64(1.0)
            out.closed_inside = int(side.sign < 0.0)
            # Outside an open cap is not penetration through a physical wall.
            # Rim crossings still require swept WALL clearance before acceptance.
            if side.sign > 0.0 and tags[side.face] == 1:
                sign = wp.float64(-1.0)
            n = sign * (q - p) / wp.max(d, wp.float64(1e-100))
            if d < wp.float64(1e-15):
                a = vertices[wall_faces[winner, 0]]
                b = vertices[wall_faces[winner, 1]]
                c = vertices[wall_faces[winner, 2]]
                n = wp.normalize(wp.cross(b - a, c - a))
            out.distance = sign * d
            out.normal = n
            out.closest = q
            out.face = winner
            out.valid = 1
    return out


@wp.kernel
def query_points(
    closed: wp.uint64,
    wall: wp.uint64,
    vertices: wp.array(dtype=wp.vec3d),
    faces: wp.array2d(dtype=int),
    tags: wp.array(dtype=int),
    wall_faces: wp.array2d(dtype=int),
    origin: wp.vec3d,
    padding: wp.float64,
    x: wp.array(dtype=wp.vec3d),
    distance: wp.array(dtype=wp.float64),
    normal: wp.array(dtype=wp.vec3d),
    face: wp.array(dtype=int),
    inside: wp.array(dtype=int),
):
    i = wp.tid()
    h = wall_query(closed, wall, vertices, faces, tags, wall_faces, origin, padding, x[i])
    distance[i] = h.distance
    normal[i] = h.normal
    face[i] = h.face
    inside[i] = h.closed_inside
    if h.valid == 0:
        face[i] = -1


@wp.kernel
def wall_project(
    closed: wp.uint64,
    wall: wp.uint64,
    vertices: wp.array(dtype=wp.vec3d),
    faces: wp.array2d(dtype=int),
    tags: wp.array(dtype=int),
    wall_faces: wp.array2d(dtype=int),
    origin: wp.vec3d,
    padding: wp.float64,
    x: wp.array(dtype=wp.vec3d),
    radius: wp.float64,
    active: wp.array(dtype=int),
    cache_x: wp.array(dtype=wp.vec3d),
    cache_d: wp.array(dtype=wp.float64),
    stats: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
):
    i = wp.tid()
    if active[i] == 0 or error[0] != 0:
        return
    lower = cache_d[i] - wp.length(x[i] - cache_x[i])
    if lower >= radius + padding:
        stats[i, 9] = stats[i, 9] + wp.float64(1.0)
        return
    h = wall_query(closed, wall, vertices, faces, tags, wall_faces, origin, padding, x[i])
    stats[i, 4] = stats[i, 4] + wp.float64(1.0)
    stats[i, 5] = stats[i, 5] + wp.float64(h.triangle_tests)
    if h.valid == 0:
        wp.atomic_or(error, 0, 512)
        return
    cache_x[i] = x[i]
    cache_d[i] = h.distance
    correction = wp.max(radius - h.distance, wp.float64(0.0))
    stats[i, 1] = wp.max(stats[i, 1], correction)
    if correction > wp.float64(0.1) * radius:
        wp.atomic_or(error, 0, 1024)
        return
    if correction > wp.float64(0.0):
        x[i] = x[i] - correction * h.normal
        stats[i, 2] = wp.max(stats[i, 2], correction)
        stats[i, 3] = stats[i, 3] + correction
        stats[i, 7] = stats[i, 7] + wp.float64(1.0)


@wp.kernel
def clip_boundaries(
    closed: wp.uint64,
    tags: wp.array(dtype=int),
    origin: wp.vec3d,
    padding: wp.float64,
    old: wp.array(dtype=wp.vec3d),
    x: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    error: wp.array(dtype=int),
    cap_origins: wp.array(dtype=wp.vec3d),
    cap_normals: wp.array(dtype=wp.vec3d),
    field_lower: wp.vec3d,
    field_upper: wp.vec3d,
    field_enabled: int,
    pending: wp.array(dtype=int),
    fractions: wp.array(dtype=wp.float64),
):
    i = wp.tid()
    pending[i] = 0
    fractions[i] = wp.float64(1.0)
    if active[i] == 0 or error[0] != 0:
        return
    p = old[i]
    q = x[i]
    dx = q - p
    alpha = wp.float64(1.0)
    label = int(0)
    for cap in range(3):
        a = wp.dot(p - cap_origins[cap], cap_normals[cap])
        b = wp.dot(q - cap_origins[cap], cap_normals[cap])
        if a <= wp.float64(0.0) and b > wp.float64(0.0):
            t = -a / (b - a)
            hit = p + t * dx
            near = wp.mesh_query_point_no_sign(closed, f32(hit - origin), float(padding * wp.float64(4.0)))
            if near.result and tags[near.face] == cap + 2 and t <= alpha:
                alpha = t
                label = cap + 2
    if field_enabled == 1:
        for axis in range(3):
            if q[axis] < field_lower[axis] and dx[axis] < wp.float64(0.0):
                t = (field_lower[axis] - p[axis]) / dx[axis]
                if t >= wp.float64(0.0) and t < alpha:
                    alpha = t
                    label = 100
            if q[axis] > field_upper[axis] and dx[axis] > wp.float64(0.0):
                t = (field_upper[axis] - p[axis]) / dx[axis]
                if t >= wp.float64(0.0) and t < alpha:
                    alpha = t
                    label = 100
    if label != 0:
        x[i] = p + alpha * dx
    pending[i] = label
    fractions[i] = alpha


@wp.kernel
def final_geometry_check(
    closed: wp.uint64,
    wall: wp.uint64,
    vertices: wp.array(dtype=wp.vec3d),
    faces: wp.array2d(dtype=int),
    tags: wp.array(dtype=int),
    wall_faces: wp.array2d(dtype=int),
    origin: wp.vec3d,
    padding: wp.float64,
    old: wp.array(dtype=wp.vec3d),
    x: wp.array(dtype=wp.vec3d),
    n: int,
    radius: wp.float64,
    active: wp.array(dtype=int),
    cache_x: wp.array(dtype=wp.vec3d),
    cache_d: wp.array(dtype=wp.float64),
    stats: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
    cap_origins: wp.array(dtype=wp.vec3d),
    cap_normals: wp.array(dtype=wp.vec3d),
    field_lower: wp.vec3d,
    field_upper: wp.vec3d,
    field_enabled: int,
    wall_enabled: int,
    pending: wp.array(dtype=int),
    fractions: wp.array(dtype=wp.float64),
    gap: wp.array(dtype=wp.float64),
):
    i = wp.tid()
    if active[i] == 0 or error[0] != 0:
        return
    p = old[i]
    q = x[i]
    if not wp.isfinite(q[0]) or not wp.isfinite(q[1]) or not wp.isfinite(q[2]):
        wp.atomic_or(error, 0, 8)
        return
    start = (i // n) * n
    for j in range(start, start + n):
        if i != j and active[j] == 1:
            if wp.length(q - x[j]) < wp.float64(1.99) * radius:
                wp.atomic_or(error, 0, 2048)
    if wall_enabled == 1:
        lower = cache_d[i] - wp.length(p - cache_x[i])
        if lower - radius > wp.length(q - p) + padding:
            stats[i, 9] = stats[i, 9] + wp.float64(1.0)
            gap[i] = lower - radius - wp.length(q - p)
        else:
            h = wall_query(closed, wall, vertices, faces, tags, wall_faces, origin, padding, q)
            stats[i, 4] = stats[i, 4] + wp.float64(1.0)
            stats[i, 5] = stats[i, 5] + wp.float64(h.triangle_tests)
            if h.valid == 0 or h.distance < wp.float64(0.99) * radius:
                wp.atomic_or(error, 0, 4096)
            gap[i] = h.distance - radius
            stats[i, 0] = wp.min(stats[i, 0], gap[i])
            cache_x[i] = q
            cache_d[i] = h.distance
            if gap[i] <= wp.float64(0.01) * radius:
                stats[i, 8] = stats[i, 8] + wp.float64(1.0)
            reach = radius + padding
            query = wp.mesh_query_aabb(
                wall,
                f32(wp.min(p, q) - origin - wp.vec3d(reach)),
                f32(wp.max(p, q) - origin + wp.vec3d(reach)),
            )
            face = int(0)
            threshold = wp.float64(0.99) * radius
            while wp.mesh_query_aabb_next(query, face):
                a = vertices[wall_faces[face, 0]]
                b = vertices[wall_faces[face, 1]]
                c = vertices[wall_faces[face, 2]]
                dd = swept_distance2(p, q, a, b, c)
                stats[i, 6] = stats[i, 6] + wp.float64(1.0)
                if dd < threshold * threshold:
                    wp.atomic_or(error, 0, 8192)


@wp.kernel
def commit_events(
    active: wp.array(dtype=int),
    pending: wp.array(dtype=int),
    fractions: wp.array(dtype=wp.float64),
    exit_label: wp.array(dtype=int),
    exit_time: wp.array(dtype=wp.float64),
    stats: wp.array2d(dtype=wp.float64),
    error: wp.array(dtype=int),
    step_count: wp.array(dtype=int),
    time_origin: wp.array(dtype=wp.float64),
    dt: wp.float64,
):
    i = wp.tid()
    if error[0] == 0 and active[i] == 1 and pending[i] != 0:
        active[i] = 0
        exit_label[i] = pending[i]
        exit_time[i] = time_origin[0] + (wp.float64(step_count[0]) + fractions[i]) * dt
        if pending[i] == 100:
            stats[i, 11] = stats[i, 11] + wp.float64(1.0)
        else:
            stats[i, 10] = stats[i, 10] + wp.float64(1.0)


@wp.kernel
def prescribed_prediction(
    x: wp.array(dtype=wp.vec3d),
    velocity: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
    dt: wp.float64,
    active: wp.array(dtype=int),
):
    i = wp.tid()
    out[i] = x[i]
    if active[i] == 1:
        out[i] = x[i] + velocity[i] * dt
