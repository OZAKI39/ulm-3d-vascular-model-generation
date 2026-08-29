from __future__ import annotations

import numpy as np
import pytest

from utils.cfd_flow.inlet_rim_audit import (
    classify_inlet_normals,
    points_in_polygon_xy,
    seam_distances_xy,
)


def _square() -> np.ndarray:
    return np.asarray(((0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)))


def test_seam_distance_uses_polygon_segments() -> None:
    points = np.asarray(((1.0, 1.0), (2.5, 1.0), (2.5, 2.5), (0.0, 0.5)))
    assert seam_distances_xy(points, _square()).tolist() == pytest.approx(
        (1.0, 0.5, np.sqrt(0.5), 0.0)
    )


def test_polygon_inside_outside_and_boundary() -> None:
    points = np.asarray(((1.0, 1.0), (2.5, 1.0), (2.0, 0.5), (-0.1, -0.1)))
    assert points_in_polygon_xy(points, _square()).tolist() == [True, False, True, False]


def test_existing_287_213_74_classification_is_preserved() -> None:
    normal_indices = np.asarray([2] * 213 + [6] * 66 + [14] * 8, dtype=np.int64)
    result = classify_inlet_normals(normal_indices)
    assert result["total_count"] == 287
    assert result["cardinal_count"] == 213
    assert result["diagonal_count"] == 74
    assert result["other_count"] == 0
