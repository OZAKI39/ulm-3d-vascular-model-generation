"""Exact rational squared distances for binary64 input geometry.

No tolerance is added and no square root participates in the safety predicate.
The segment/triangle minimum is attained at a crossing, a segment endpoint, or
one of the three finite triangle edges. Parallel interior minima also attain
the same value on one of these boundaries. Degenerate triangles reduce to edges.
Used for difficult predicates and independent diagnosis, not dense sampling.
"""

from fractions import Fraction as F
from decimal import Decimal, localcontext


def vector(x):
    return tuple(v if isinstance(v, F) else F(float(v)) for v in x)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(a, s):
    return tuple(x * s for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def point_segment(p, a, b):
    v = sub(b, a)
    vv = dot(v, v)
    t = max(F(0), min(F(1), dot(sub(p, a), v) / vv)) if vv else F(0)
    q = add(a, mul(v, t))
    d = sub(p, q)
    return dot(d, d), t, q


def point_triangle(p, a, b, c):
    e, f = sub(b, a), sub(c, a)
    aa, ab, bb = dot(e, e), dot(e, f), dot(f, f)
    v = sub(p, a)
    ve, vf = dot(v, e), dot(v, f)
    den = aa * bb - ab * ab
    if den:
        s, t = (ve * bb - vf * ab) / den, (vf * aa - ve * ab) / den
        if s >= 0 and t >= 0 and s + t <= 1:
            q = add(a, add(mul(e, s), mul(f, t)))
            d = sub(p, q)
            return dot(d, d), q, "face", (1 - s - t, s, t)
    rows = []
    for i, (u, v) in enumerate([(a, b), (b, c), (c, a)]):
        d, t, q = point_segment(p, u, v)
        weights = [F(0)] * 3
        weights[i] = 1 - t
        weights[(i + 1) % 3] = t
        rows.append((d, q, "edge" if 0 < t < 1 else "vertex", tuple(weights)))
    return min(rows, key=lambda z: z[0])


def segment_segment(p, q, a, b):
    u, v, w = sub(q, p), sub(b, a), sub(p, a)
    aa, bb, cc, dd, ee = dot(u, u), dot(u, v), dot(v, v), dot(u, w), dot(v, w)
    rows = []
    for s, pt in [(F(0), p), (F(1), q)]:
        d, t, close = point_segment(pt, a, b)
        rows.append((d, s, t, pt, close))
    for t, pt in [(F(0), a), (F(1), b)]:
        d, s, close = point_segment(pt, p, q)
        rows.append((d, s, t, close, pt))
    den = aa * cc - bb * bb
    if den:
        s, t = (bb * ee - cc * dd) / den, (aa * ee - bb * dd) / den
        if 0 <= s <= 1 and 0 <= t <= 1:
            left, right = add(p, mul(u, s)), add(a, mul(v, t))
            d = sub(left, right)
            rows.append((dot(d, d), s, t, left, right))
    return min(rows, key=lambda z: z[0])


def segment_triangle(p, q, a, b, c):
    p, q, a, b, c = map(vector, [p, q, a, b, c])
    n = cross(sub(b, a), sub(c, a))
    u = sub(q, p)
    den = dot(n, u)
    if den:
        s = dot(n, sub(a, p)) / den
        if 0 <= s <= 1:
            hit = add(p, mul(u, s))
            d, close, feature, weights = point_triangle(hit, a, b, c)
            if not d:
                return dict(
                    distance2=d,
                    fraction=s,
                    point_path=hit,
                    point_triangle=close,
                    feature=feature,
                    barycentric=weights,
                )
    rows = []
    for s, pt in [(F(0), p), (F(1), q)]:
        d, close, feature, weights = point_triangle(pt, a, b, c)
        rows.append(
            dict(
                distance2=d,
                fraction=s,
                point_path=pt,
                point_triangle=close,
                feature=feature,
                barycentric=weights,
            )
        )
    for i, (left, right) in enumerate([(a, b), (b, c), (c, a)]):
        d, s, t, pt, close = segment_segment(p, q, left, right)
        weights = [F(0)] * 3
        weights[i] = 1 - t
        weights[(i + 1) % 3] = t
        rows.append(
            dict(
                distance2=d,
                fraction=s,
                point_path=pt,
                point_triangle=close,
                feature="edge" if 0 < t < 1 else "vertex",
                barycentric=tuple(weights),
            )
        )
    return min(rows, key=lambda z: z["distance2"])


def gap_decimal(distance2, radius):
    with localcontext() as ctx:
        ctx.prec = 80
        d = (Decimal(distance2.numerator) / Decimal(distance2.denominator)).sqrt()
        return str(d - Decimal.from_float(float(radius)))


def record(result, radius, tolerance):
    d = result["distance2"]
    threshold = F(float(radius)) - F(float(tolerance))
    gap = gap_decimal(d, radius)
    return dict(
        classification="SAFE" if threshold <= 0 or d >= threshold * threshold else "UNSAFE",
        predicate="exact rational min squared distance >= exact (radius - tolerance)^2",
        gap_m=float(gap),
        gap_decimal_m=gap,
        minimum_fraction=float(result["fraction"]),
        feature=result["feature"],
        barycentric=[float(x) for x in result["barycentric"]],
        closest_path_point_m=[float(x) for x in result["point_path"]],
        closest_triangle_point_m=[float(x) for x in result["point_triangle"]],
        squared_margin_numerator=str((d - threshold * threshold).numerator),
        squared_margin_denominator=str((d - threshold * threshold).denominator),
    )
