from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from utils.rodent_vasculature.roi_surface import (
    TUBE_SIDES,
    _interpolate_branch,
    build_roi_display_tubes,
)


def _roi(points, radii, edges):
    points = np.asarray(points, dtype=float)
    radii = np.asarray(radii, dtype=float)
    edges = np.asarray(edges, dtype=np.int64)
    return SimpleNamespace(
        local_node_ids=np.arange(len(points)),
        local_node_positions_um=points,
        local_node_radius_um=radii,
        local_edges=edges,
        local_edge_points_um=points[edges],
        local_edge_radius_um=radii[edges],
    )


def test_interpolation_keeps_exact_samples_without_radius_overshoot():
    points = np.asarray([[0, 0, 0], [3, 0, 0], [4, 2, 0], [7, 3, 1], [10, 3, 1]], dtype=float)
    radii = np.asarray([1, 3, 0.4, 2, 1], dtype=float)
    sampled, radius, knots, arc = _interpolate_branch(points, radii)
    np.testing.assert_array_equal(sampled[knots], points)
    np.testing.assert_array_equal(radius[knots], radii)
    assert np.all(np.diff(arc) > 0)
    for i, (start, end) in enumerate(zip(knots[:-1], knots[1:])):
        assert np.all(radius[start:end + 1] >= min(radii[i:i + 2]))
        assert np.all(radius[start:end + 1] <= max(radii[i:i + 2]))
    # Interpolation preserves a changing slope through an increasing knot,
    # unlike independent per-edge easing with a forced zero slope at every node.
    straight = np.column_stack((np.arange(3, dtype=float), np.zeros((3, 2))))
    _, increasing, ix, parameter = _interpolate_branch(straight, np.asarray([1., 2., 3.]))
    slopes = np.diff(increasing) / np.diff(parameter)
    np.testing.assert_allclose(slopes[ix[1]-1:ix[1]+1], [1, 1], atol=1e-12)


def test_continuous_tube_has_exact_geometric_radii_and_no_internal_caps():
    roi = _roi([[0, 0, 0], [5, 0, 0], [10, 3, 0], [15, 3, 2]], [1, 2, 0.5, 1.5], [[0, 1], [1, 2], [2, 3]])
    before = {key: value.copy() for key, value in vars(roi).items()}
    line, tube = build_roi_display_tubes(roi)
    assert line.n_lines == 1
    assert tube.GetNumberOfPolys() == 2
    assert tube.n_open_edges == 0
    source_ids = np.asarray(line["source_local_node_id"])
    for node in range(4):
        sample, = np.flatnonzero(source_ids == node)
        np.testing.assert_array_equal(line.points[sample], roi.local_node_positions_um[node])
        assert line["radius_um"][sample] == roi.local_node_radius_um[node]
        assert np.count_nonzero(tube["display_sample_id"] == sample) == TUBE_SIDES
    sample_ids = np.asarray(tube["display_sample_id"])
    actual_radii = np.linalg.norm(tube.points - line.points[sample_ids], axis=1)
    np.testing.assert_allclose(actual_radii, line["radius_um"][sample_ids], rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(tube["diameter_um"], 2 * tube["radius_um"])
    for key, value in before.items():
        np.testing.assert_array_equal(getattr(roi, key), value)


def test_bifurcation_preserves_branches_and_has_no_junction_end_disk():
    roi = _roi([[0, 0, 0], [4, 0, 0], [8, 0, 0], [12, 4, 0], [12, -4, 0]], [2, 2, 1.5, 1, 0.5], [[0, 1], [1, 2], [2, 3], [2, 4]])
    line, tube = build_roi_display_tubes(roi)
    assert line.n_lines == 3
    assert tube.GetNumberOfPolys() == 3
    assert np.count_nonzero(line["source_local_node_id"] == 2) == 3
    faces = tube.faces.reshape((-1, TUBE_SIDES + 1))[:, 1:]
    cap_nodes = tube["source_local_node_id"][faces]
    assert set(np.unique(cap_nodes)) == {0, 3, 4}
    assert np.all(tube["radius_um"][tube["source_local_node_id"] == 2] == 1.5)


def test_closed_branch_preserves_knots_without_end_caps():
    roi = _roi([[0, 0, 0], [5, 0, 0], [5, 5, 0], [0, 5, 0]], [1, 2, 1, 0.5], [[0, 1], [1, 2], [2, 3], [3, 0]])
    line, tube = build_roi_display_tubes(roi)
    assert line.n_lines == 1
    assert tube.GetNumberOfPolys() == 0
    assert set(line["source_local_node_id"]) == {-1, 0, 1, 2, 3}
    for node in range(4):
        assert np.all(tube["radius_um"][tube["source_local_node_id"] == node] == roi.local_node_radius_um[node])


def test_conflicting_radii_are_rejected_instead_of_averaged():
    roi = _roi([[0, 0, 0], [5, 0, 0]], [1, 2], [[0, 1]])
    roi.local_edge_radius_um[0, 1] = 3
    with pytest.raises(ValueError, match="refusing to average"):
        build_roi_display_tubes(roi)
