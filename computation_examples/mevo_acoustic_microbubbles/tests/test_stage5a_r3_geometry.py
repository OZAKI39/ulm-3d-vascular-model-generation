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
