from __future__ import annotations

import numpy as np

from utils.cfd_flow.dimensionless_geometry_kernel import dimensionless_ray_triangle


def _triangle(z: float) -> np.ndarray:
    return np.asarray(((-1.0, -1.0, z), (2.0, -1.0, z), (-1.0, 2.0, z)))


def test_orthogonal_q_is_fraction_of_non_normalized_link() -> None:
    result = dimensionless_ray_triangle(
        np.zeros(3), np.asarray((0.0, 0.0, 2.0)), _triangle(0.75)
    )

    assert result.hit is True
    assert result.classification == "HIT_INTERIOR"
    assert np.isclose(result.q, 0.375)


def test_diagonal_link_q_preserves_treelm_fraction_contract() -> None:
    triangle = np.asarray(((0.5, -1.0, 0.0), (0.5, 2.0, 0.0), (0.5, -1.0, 2.0)))

    result = dimensionless_ray_triangle(
        np.zeros(3), np.asarray((1.0, 1.0, 0.0)), triangle
    )

    assert result.hit is True
    assert np.isclose(result.q, 0.5)


def test_q_below_and_above_half_are_not_truncated() -> None:
    below = dimensionless_ray_triangle(
        np.zeros(3), np.asarray((0, 0, 1)), _triangle(0.2)
    )
    above = dimensionless_ray_triangle(
        np.zeros(3), np.asarray((0, 0, 1)), _triangle(0.8)
    )

    assert below.hit and np.isclose(below.q, 0.2)
    assert above.hit and np.isclose(above.q, 0.8)
