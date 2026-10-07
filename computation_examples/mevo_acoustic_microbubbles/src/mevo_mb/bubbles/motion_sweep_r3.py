"""R3 continuous WALL certificate: separating bounds, then exact predicates.

For ANY nonzero direction n, min(segment.n)-max(triangle.n), divided by |n|,
is a lower distance bound when positive. All arithmetic in that bound is rounded
outward with nextafter. An approximate closest-point direction can only weaken
the certificate, not create a false PASS. Unresolved candidates use exact rational
segment/finite-triangle squared distance. No proof query alters the physical path.
"""

from fractions import Fraction as F
import numpy as np
from numba import njit
from .motion_sweep_r2 import WallSweep
from .geometry_exact_r3 import segment_triangle, record, vector, sub, point_segment


@njit(cache=True)
def closest_point(p, tri):
    a, b, c = tri
    e = b - a
    f = c - a
    v = p - a
    aa = np.dot(e, e)
    ab = np.dot(e, f)
    bb = np.dot(f, f)
    den = aa * bb - ab * ab
    if den > 0:
        s = (np.dot(v, e) * bb - np.dot(v, f) * ab) / den
        t = (np.dot(v, f) * aa - np.dot(v, e) * ab) / den
        if s >= 0 and t >= 0 and s + t <= 1:
            return a + s * e + t * f
    best = np.inf
    answer = a.copy()
    for i in range(3):
        left = tri[i]
        right = tri[(i + 1) % 3]
        d = right - left
        dd = np.dot(d, d)
        t = max(0.0, min(1.0, np.dot(p - left, d) / dd)) if dd > 0 else 0.0
        q = left + t * d
        dist = np.dot(p - q, p - q)
        if dist < best:
            best = dist
            answer = q
    return answer


@njit(cache=True)
def projected_interval(x, origin, n):
    lo = 0.0
    hi = 0.0
    for k in range(3):
        d = x[k] - origin[k]
        dl = np.nextafter(d, -np.inf)
        dh = np.nextafter(d, np.inf)
        a = dl * n[k]
        b = dh * n[k]
        pl = np.nextafter(min(a, b), -np.inf)
        ph = np.nextafter(max(a, b), np.inf)
        lo = np.nextafter(lo + pl, -np.inf)
        hi = np.nextafter(hi + ph, np.inf)
    return lo, hi


@njit(cache=True)
def separation_bound(p, q, tri, n):
    norm2 = 0.0
    for k in range(3):
        norm2 = np.nextafter(norm2 + np.nextafter(n[k] * n[k], np.inf), np.inf)
    if norm2 <= 0 or not np.isfinite(norm2):
        return -np.inf
    origin = tri[0]
    seg = min(projected_interval(p, origin, n)[0], projected_interval(q, origin, n)[0])
    face = -np.inf
    for i in range(3):
        face = max(face, projected_interval(tri[i], origin, n)[1])
    gap = np.nextafter(seg - face, -np.inf)
    if gap <= 0:
        return 0.0
    return np.nextafter(gap / np.nextafter(np.sqrt(norm2), np.inf), -np.inf)


@njit(cache=True)
def finite_triangle_lower_bound(p, q, tri):
    middle = p + (q - p) / 2
    n = middle - closest_point(middle, tri)
    normal = np.cross(tri[1] - tri[0], tri[2] - tri[0])
    return max(
        separation_bound(p, q, tri, n),
        separation_bound(p, q, tri, normal),
        separation_bound(p, q, tri, -normal),
    )


@njit(cache=True)
def pair_lower_bounds(x, y, r, ii, jj):
    result = np.empty(len(ii))
    for k in range(len(ii)):
        i, j = ii[k], jj[k]
        delta = x[j] - x[i]
        change = (y[j] - x[j]) - (y[i] - x[i])
        den = np.dot(change, change)
        at = max(0.0, min(1.0, -np.dot(delta, change) / den)) if den else 0.0
        n = delta + at * change
        norm2 = 0.0
        for a in range(3):
            norm2 = np.nextafter(norm2 + np.nextafter(n[a] * n[a], np.inf), np.inf)
        support = min(projected_interval(x[j], x[i], n)[0], projected_interval(y[j], y[i], n)[0])
        lower = (
            np.nextafter(support / np.nextafter(np.sqrt(norm2), np.inf), -np.inf)
            if support > 0 and norm2 > 0
            else 0.0
        )
        result[k] = np.nextafter(np.nextafter(lower - r[i], -np.inf) - r[j], -np.inf)
    return result


class WallSweepR3(WallSweep):
    def __init__(self, geometry, wall_path):
        super().__init__(geometry, wall_path)
        self.triangle_low = self.triangles.min(1)
        self.triangle_high = self.triangles.max(1)
        self.rational_calls = 0
        self.last_details = []
        finite_triangle_lower_bound(np.zeros(3), np.ones(3), self.triangles[0])

    def candidates(self, p, q, reach):
        low = np.nextafter(np.minimum(p, q) - reach, -np.inf)
        high = np.nextafter(np.maximum(p, q) + reach, np.inf)
        self.locator.FindCellsWithinBounds(np.column_stack([low, high]).ravel(), self.ids)
        cells = np.array([self.ids.GetId(k) for k in range(self.ids.GetNumberOfIds())], int)
        return cells[
            np.all(self.triangle_high[cells] >= low, axis=1)
            & np.all(self.triangle_low[cells] <= high, axis=1)
        ]

    def __call__(self, x, y, r, tolerance):
        self.calls += 1
        self.last_details = []
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            return dict(pass_gate=False, reason="NONFINITE_INPUT")
        initial = self.geometry.wall_gap_and_normal(x, r)[0]
        final = self.geometry.wall_gap_and_normal(y, r)[0]
        guard = 64 * np.finfo(float).eps * max(np.max(abs(x)), np.max(abs(y)), r.max())
        bound = initial - np.linalg.norm(y - x, axis=1) - guard
        # Signed initial-centre membership is required. Caps do not enter WALL
        # candidate sets. This +5 ms protocol starts inside and has no port exit.
        if np.any(initial + r < -guard):
            return dict(pass_gate=False, reason="INITIAL_SIDE_NOT_CERTIFIED")
        count = exact = 0
        wall_safe = True
        for i in range(len(x)):
            cells = self.candidates(x[i], y[i], np.nextafter(r[i] + tolerance, np.inf))
            value = float(tolerance)
            for cell in cells:
                tri = self.triangles[cell]
                lower = np.nextafter(finite_triangle_lower_bound(x[i], y[i], tri) - r[i], -np.inf)
                detail = None
                if lower < -tolerance:
                    result = segment_triangle(x[i], y[i], *tri)
                    detail = record(result, r[i], tolerance)
                    exact += 1
                    self.rational_calls += 1
                    if detail["classification"] == "SAFE":
                        # The threshold itself is now an EXACTLY proved bound,
                        # not a clipped raw gap; raw high-precision gap is retained.
                        lower = -tolerance
                    else:
                        lower = detail["gap_m"]
                        wall_safe = False
                    self.last_details.append(
                        dict(
                            particle_index=int(i),
                            wall_triangle_id=int(cell),
                            directed_float_lower_gap_m=float(
                                np.nextafter(finite_triangle_lower_bound(x[i], y[i], tri) - r[i], -np.inf)
                            ),
                            **detail,
                        )
                    )
                value = min(value, lower)
            bound[i] = value
            count += len(cells)
        ii, jj = np.triu_indices(len(x), 1)
        pair_bounds = pair_lower_bounds(x, y, r, ii, jj)
        pair = float(np.min(pair_bounds, initial=np.inf))
        if pair < -tolerance:
            ii, jj = np.triu_indices(len(x), 1)
            delta = x[jj] - x[ii]
            change = (y - x)[jj] - (y - x)[ii]
            den = np.einsum("ij,ij->i", change, change)
            at = np.clip(-np.einsum("ij,ij->i", delta, change) / np.maximum(den, np.finfo(float).tiny), 0, 1)
            gaps = np.linalg.norm(delta + at[:, None] * change, axis=1) - r[ii] - r[jj]
            pair = float(np.min(gaps - guard))
            for k in np.flatnonzero(pair_bounds < -tolerance):
                i, j = ii[k], jj[k]
                p = sub(vector(x[j]), vector(x[i]))
                q = sub(vector(y[j]), vector(y[i]))
                distance2, _, _ = point_segment((F(0), F(0), F(0)), p, q)
                target = F(float(r[i])) + F(float(r[j])) - F(float(tolerance))
                if distance2 < target * target:
                    return dict(
                        pass_gate=False,
                        reason="EXACT_PAIR_UNSAFE",
                        particle_pair_indices=[int(i), int(j)],
                        minimum_pair_swept_gap_m=float(gaps[k]),
                        minimum_wall_swept_gap_bound_m=float(bound.min()),
                        exact_wall_details=self.last_details,
                    )
            # Every ambiguous pair was proved; all remaining pair bounds pass.
            pair = -tolerance
        return dict(
            pass_gate=bool(wall_safe and min(pair, bound.min()) >= -tolerance),
            minimum_pair_swept_gap_m=float(pair) if np.isfinite(pair) else None,
            pair_gap_semantics="certified continuous lower bound, not a clipped raw gap"
            if len(x) > 1
            else "NOT_APPLICABLE: fewer than two particles",
            minimum_wall_swept_gap_bound_m=float(bound.min()),
            minimum_wall_endpoint_gap_m=float(final.min()),
            wall_certificate="signed initial membership; expanded-AABB coverage / outward-rounded separating bound / exact rational finite-triangle predicate",
            wall_triangle_candidates=count,
            exact_wall_predicates=exact,
            exact_wall_details=self.last_details,
            projection_max_m=0.0,
            roundoff_guard_m=guard,
        )
