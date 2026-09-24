from __future__ import annotations

import math

import numpy as np

from utils.cfd_flow.wall_qvalue_oracle import ray_cylinder_fraction


def test_axis_aligned_ray_cylinder_q() -> None:
    origins = np.asarray(((0.0, 0.5, 0.0), (0.0, 0.0, -0.25)))
    directions = np.asarray(((0.0, 1.0, 0.0), (0.0, 0.0, -0.5)))
    q = ray_cylinder_fraction(origins, directions, axis=(1.0, 0.0, 0.0), radius_m=1.0)
    assert np.allclose(q, (0.5, 1.5))


def test_rotated_ray_cylinder_q() -> None:
    axis = np.asarray((1.0, 1.0, 0.0)) / math.sqrt(2.0)
    radial = np.asarray((-1.0, 1.0, 0.0)) / math.sqrt(2.0)
    origins = np.asarray((0.25 * radial,))
    directions = np.asarray((0.5 * radial,))
    q = ray_cylinder_fraction(origins, directions, axis=axis, radius_m=1.0)
    assert np.allclose(q, (1.5,))


def test_parallel_ray_has_no_side_intersection() -> None:
    q = ray_cylinder_fraction(
        np.asarray(((0.0, 0.5, 0.0),)),
        np.asarray(((1.0, 0.0, 0.0),)),
        axis=(1.0, 0.0, 0.0),
        radius_m=1.0,
    )
    assert math.isnan(float(q[0]))
