from fractions import Fraction as F
import numpy as np
import pytest
from mevo_mb.bubbles.geometry_exact_r3 import segment_triangle, record


@pytest.mark.parametrize(
    "p,q,expected,feature",
    [
        ([0.25, 0.25, 1], [0.5, 0.25, 1], F(1), "face"),
        ([0.25, 0.25, 1], [0.25, 0.25, -1], F(0), "face"),
        ([0.5, -1, 1], [0.5, -1, -1], F(1), "edge"),
        ([-1, -1, 1], [-1, -1, -1], F(2), "vertex"),
        ([-1, 0.25, 1], [2, 0.25, 1], F(1), None),
        ([0.25, 0.25, 1], [0.25, 0.25, 1], F(1), "face"),
    ],
)
def test_exact_finite_triangle(p, q, expected, feature):
    out = segment_triangle(p, q, [0, 0, 0], [1, 0, 0], [0, 1, 0])
    assert out["distance2"] == expected
    if feature:
        assert out["feature"] == feature


def test_finite_triangle_is_not_an_infinite_wall():
    out = segment_triangle([3, 3, 1], [3, 3, -1], [0, 0, 0], [1, 0, 0], [0, 1, 0])
    assert out["distance2"] == F(25, 2)
    assert record(out, 1.0, 1e-11)["classification"] == "SAFE"


def test_exact_threshold_preserves_negative_tolerance_band():
    tolerance = 2.0**-20
    radius = 1.0
    gap = radius - tolerance
    tri = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    for z, expected in [
        (gap, "SAFE"),
        (np.nextafter(gap, -np.inf), "UNSAFE"),
        (np.nextafter(gap, np.inf), "SAFE"),
    ]:
        out = segment_triangle([0.25, 0.25, z], [0.5, 0.25, z], *tri)
        result = record(out, radius, tolerance)
        assert result["classification"] == expected
        assert result["gap_m"] < 0


def test_endpoint_clear_but_middle_crosses():
    out = segment_triangle([0.25, 0.25, 1], [0.25, 0.25, -1], [0, 0, 0], [1, 0, 0], [0, 1, 0])
    assert out["fraction"] == F(1, 2)
    assert record(out, 0.1, 1e-11)["classification"] == "UNSAFE"


def test_binary_exact_common_translation():
    tri = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    p = np.array([0.5, -1, 1])
    q = np.array([0.5, -1, -1])
    shift = np.array([2.0**20, -(2.0**18), 2.0**16])
    a = segment_triangle(p, q, *tri)
    b = segment_triangle(p + shift, q + shift, *(tri + shift))
    assert a["distance2"] == b["distance2"] and a["fraction"] == b["fraction"]


def test_degenerate_triangle_and_parallel_segments():
    out = segment_triangle([0.25, 1, 0], [0.5, 1, 0], [0, 0, 0], [1, 0, 0], [2, 0, 0])
    assert out["distance2"] == 1


def test_directed_separating_bound_never_exceeds_exact_distance():
    from mevo_mb.bubbles.motion_sweep_r3 import finite_triangle_lower_bound

    rng = np.random.default_rng(38219)
    for _ in range(60):
        tri = rng.normal(size=(3, 3)) * 1e-4 + np.array([0.002, 0.001, -0.001])
        p = tri.mean(0) + rng.normal(size=3) * 1e-6
        q = p + rng.normal(size=3) * 1e-7
        lower = finite_triangle_lower_bound(p, q, tri)
        exact = segment_triangle(p, q, *tri)["distance2"]
        assert lower >= 0 and F(float(lower)) ** 2 <= exact


def test_all_pair_directed_bounds_cover_exact_relative_motion():
    from mevo_mb.bubbles.motion_sweep_r3 import pair_lower_bounds
    from mevo_mb.bubbles.geometry_exact_r3 import vector, sub, point_segment

    rng = np.random.default_rng(902)
    x = 0.01 + rng.normal(size=(12, 3)) * 1e-5
    y = x + rng.normal(size=x.shape) * 1e-7
    r = np.full(12, 1.78e-6)
    ii, jj = np.triu_indices(len(x), 1)
    lower = pair_lower_bounds(x, y, r, ii, jj)
    for k, (i, j) in enumerate(zip(ii, jj)):
        d, _, _ = point_segment((F(0),) * 3, sub(vector(x[j]), vector(x[i])), sub(vector(y[j]), vector(y[i])))
        dist = F(float(lower[k])) + F(float(r[i])) + F(float(r[j]))
        assert dist <= 0 or dist * dist <= d
