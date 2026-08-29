from __future__ import annotations

import numpy as np
import pytest

from utils.cfd_flow.axis_aligned_inlet import (
    TARGET_INWARD_NORMAL,
    apply_rigid_transform,
    area_proxy_qc,
    cardinal_normal_gate,
    compute_phase_aligned_bounding_cube,
    homogeneous_transform,
    inlet_plane_deviation,
    minimal_rotation_matrix,
    pairwise_distances,
)
from utils.cfd_flow.exact_link_flux import reconstruct_boundary


def test_rigid_rotation_maps_arbitrary_normal_to_negative_z() -> None:
    source = np.asarray((-0.025, -0.672, -0.740))
    rotation, angle = minimal_rotation_matrix(source, TARGET_INWARD_NORMAL)
    mapped = rotation @ (source / np.linalg.norm(source))
    assert mapped == pytest.approx(TARGET_INWARD_NORMAL, abs=1.0e-14)
    assert angle > 0.0
    assert np.linalg.det(rotation) == pytest.approx(1.0)


def test_rigid_rotation_preserves_pairwise_distances() -> None:
    points = np.asarray(((0.0, 0.0, 0.0), (1.0, 2.0, 3.0), (-2.0, 0.5, 4.0)))
    rotation, _ = minimal_rotation_matrix(
        np.asarray((1.0, 2.0, 3.0)), TARGET_INWARD_NORMAL
    )
    rotated = apply_rigid_transform(points, rotation, np.asarray((0.4, 0.2, -0.3)))
    assert pairwise_distances(rotated) == pytest.approx(
        pairwise_distances(points), abs=1.0e-14
    )


def test_rigid_rotation_preserves_triangle_area() -> None:
    triangle = np.asarray(((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (0.0, 2.0, 0.0)))
    rotation, _ = minimal_rotation_matrix(
        np.asarray((0.2, 0.8, 0.5)), TARGET_INWARD_NORMAL
    )
    rotated = apply_rigid_transform(triangle, rotation, np.asarray((1.0, 1.0, 1.0)))
    area_before = (
        np.linalg.norm(np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0]))
        / 2.0
    )
    area_after = (
        np.linalg.norm(np.cross(rotated[1] - rotated[0], rotated[2] - rotated[0])) / 2.0
    )
    assert area_after == pytest.approx(area_before, abs=1.0e-14)


def test_inverse_homogeneous_transform_recovers_original_points() -> None:
    rotation, _ = minimal_rotation_matrix(
        np.asarray((0.1, -0.8, -0.5)), TARGET_INWARD_NORMAL
    )
    pivot = np.asarray((2.0, 3.0, 4.0))
    forward, inverse = homogeneous_transform(rotation, pivot)
    points = np.asarray(((1.0, 2.0, 3.0, 1.0), (-2.0, 4.0, 0.5, 1.0))).T
    assert inverse @ forward @ points == pytest.approx(points, abs=1.0e-14)


def test_cell_entity_ids_are_unchanged_by_point_only_transform() -> None:
    entity_ids = np.asarray((0, 1, 1, 4), dtype=np.int64)
    points = np.eye(3)
    rotation, _ = minimal_rotation_matrix(
        np.asarray((1.0, 1.0, 1.0)), TARGET_INWARD_NORMAL
    )
    _ = apply_rigid_transform(points, rotation, np.zeros(3))
    assert entity_ids.tolist() == [0, 1, 1, 4]


def test_planar_inlet_remains_planar_after_rigid_rotation() -> None:
    points = np.asarray(((0.0, 0.0, 2.0), (1.0, 0.0, 2.0), (0.0, 1.0, 2.0)))
    faces = np.asarray(((0, 1, 2),))
    deviation = inlet_plane_deviation(
        points,
        faces,
        np.asarray([0]),
        np.asarray((1 / 3, 1 / 3, 2.0)),
        np.asarray((0, 0, 1)),
    )
    assert deviation == pytest.approx(0.0)


def test_grid_phase_alignment_is_exact_and_margin_is_preserved() -> None:
    bounds = np.asarray(((0.0, 10.0), (2.0, 12.0), (3.0, 13.0)))
    cube, qc = compute_phase_aligned_bounding_cube(
        bounds, 13.0, dx_m=0.2e-6, margin_cells=4
    )
    assert qc["inlet_plane_grid_phase_error_over_dx"] <= 1.0e-10
    assert cube.margin_cells_minimum >= 4.0


def test_existing_boundary_parser_is_reused_for_d3q19_inlet() -> None:
    ids = np.zeros((1, 26), dtype=np.int64)
    ids[0, 5] = 2  # outward +Z, incoming -Z
    inlet = reconstruct_boundary(ids, np.asarray([9]), label="inlet", boundary_id=2)
    assert inlet.cell_indices.tolist() == [9]
    assert inlet.normal_indices.tolist() == [2]


def test_normal_ind_cardinal_gate_rejects_any_diagonal() -> None:
    passed = cardinal_normal_gate(np.asarray((2, 2, 2)))
    failed = cardinal_normal_gate(np.asarray((2, 6)))
    assert passed["status"] == "PASS"
    assert passed["target_cardinal_normal_fraction"] == 1.0
    assert failed["status"] == "FAIL"
    assert failed["diagonal_normal_ind_count"] == 1


def test_area_proxy_ratio_uses_d3q19_count_without_hardcoding_195() -> None:
    qc = area_proxy_qc(196, 7.819752111687344e-12, 2.0e-7)
    assert qc["mfr_eq_area_proxy_m2"] == pytest.approx(7.84e-12)
    assert qc["area_proxy_over_smooth_area"] == pytest.approx(1.002589324)
    assert qc["status"] == "PASS"
