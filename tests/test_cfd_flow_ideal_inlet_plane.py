from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from utils.cfd_flow.axis_aligned_inlet import (
    EXPECTED_DX_M,
    EXPECTED_SMOOTH_INLET_AREA_M2,
    area_proxy_qc,
    cardinal_normal_gate,
)
from utils.cfd_flow.geometry import BoundingCube
from utils.cfd_flow.ideal_inlet_plane import (
    PREVIOUS_FLUID_CELL_COUNT,
    assess_plane_coverage_safety,
    build_grid_snapped_rectangle,
    connected_fluid_region_count,
    file_snapshot,
    fluid_count_comparison_qc,
    generate_ideal_plane_seeder_lua,
    plane_vertices,
    write_two_triangle_plane_stl,
)


def _rectangle() -> tuple[dict[str, object], np.ndarray]:
    physical = np.asarray(((1.05e-6, 2.03e-6), (3.07e-6, 4.01e-6)))
    rectangle = build_grid_snapped_rectangle(
        physical,
        np.asarray((0.0, 0.0)),
        dx_m=EXPECTED_DX_M,
        margin_cells=2,
    )
    return rectangle, physical


def test_plane_z_normal_margin_and_outward_grid_snap(tmp_path: Path) -> None:
    rectangle, physical = _rectangle()
    expanded = np.asarray(rectangle["expanded_xy_bounds_before_snap_m"])
    snapped = np.asarray(rectangle["snapped_xy_bounds_m"])
    assert physical[:, 0] - expanded[:, 0] == pytest.approx(2 * EXPECTED_DX_M)
    assert expanded[:, 1] - physical[:, 1] == pytest.approx(2 * EXPECTED_DX_M)
    assert np.all(snapped[:, 0] <= expanded[:, 0])
    assert np.all(snapped[:, 1] >= expanded[:, 1])
    assert rectangle["maximum_grid_face_phase_error_over_dx"] == pytest.approx(0.0)

    z = 8.2e-6
    vertices = plane_vertices(snapped, z)
    source_copy = vertices.copy()
    qc = write_two_triangle_plane_stl(tmp_path / "numerical_inlet_plane.stl", vertices)
    assert np.array_equal(vertices, source_copy)
    assert np.all(vertices[:, 2] == z)
    assert qc["triangle_count"] == 2
    assert qc["plane_normal"] == [0.0, 0.0, 1.0]
    assert qc["triangle_orientation_consistent"] is True


def test_seeder_replaces_only_the_physical_inlet_reference() -> None:
    cube = BoundingCube(np.zeros(3), 102.4e-6, 9, 512, 4.0)
    text = generate_ideal_plane_seeder_lua(
        seed_m=np.asarray((1.0e-6, 2.0e-6, 3.0e-6)),
        cube=cube,
        rectangle_xy_bounds_m=np.asarray(((1.0e-6, 3.0e-6), (2.0e-6, 5.0e-6))),
        plane_z_m=8.0e-6,
    )
    assert "geometry_solver_m/inlet.stl" not in text
    assert "label = 'inlet'" in text
    assert "kind = 'canoND'" in text
    assert "vec = {" in text
    for label in ("wall", "outlet_01", "outlet_02", "outlet_03"):
        assert f"geometry_solver_m/{label}.stl" in text


def test_geometry_safety_allows_seam_and_rejects_unrelated_section() -> None:
    points = np.asarray(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, -1.0),
            (3.0, 3.0, -1.0),
            (3.0, 3.0, 1.0),
            (4.0, 3.0, 0.0),
        )
    )
    faces = np.asarray(((0, 1, 2), (0, 1, 3), (4, 5, 6)))
    safe = assess_plane_coverage_safety(
        points,
        faces,
        np.asarray((0,)),
        plane_z_m=0.0,
        rectangle_xy_bounds_m=np.asarray(((-0.2, 1.2), (-0.2, 1.2))),
    )
    unsafe = assess_plane_coverage_safety(
        points,
        faces,
        np.asarray((0,)),
        plane_z_m=0.0,
        rectangle_xy_bounds_m=np.asarray(((-0.2, 4.2), (-0.2, 3.2))),
    )
    assert safe["status"] == "PASS"
    assert unsafe["unsafe_component_count"] == 1


def test_source_snapshot_remains_unchanged(tmp_path: Path) -> None:
    source = tmp_path / "rotated_geometry.stl"
    source.write_bytes(b"frozen-rotated-geometry")
    before = file_snapshot((source,))
    _ = build_grid_snapped_rectangle(
        np.asarray(((0.0, 1.0), (0.0, 1.0))),
        np.zeros(2),
        dx_m=0.1,
        margin_cells=2,
    )
    assert file_snapshot((source,)) == before


def test_physical_area_and_area_proxy_contract() -> None:
    assert EXPECTED_SMOOTH_INLET_AREA_M2 == 7.819752111687344e-12
    qc = area_proxy_qc(196, EXPECTED_SMOOTH_INLET_AREA_M2, EXPECTED_DX_M)
    assert qc["mfr_eq_area_proxy_m2"] == pytest.approx(7.84e-12)
    assert qc["area_proxy_over_smooth_area"] == pytest.approx(1.002589324)
    assert qc["status"] == "PASS"


def test_normal_gate_requires_all_negative_z() -> None:
    assert cardinal_normal_gate(np.asarray((2, 2, 2)))["status"] == "PASS"
    failed = cardinal_normal_gate(np.asarray((2, 6)))
    assert failed["status"] == "FAIL"
    assert failed["target_cardinal_normal_fraction"] == 0.5
    assert failed["diagonal_normal_ind_count"] == 1


def test_fluid_count_gate_is_half_percent() -> None:
    passed = fluid_count_comparison_qc(PREVIOUS_FLUID_CELL_COUNT + 1_000)
    failed = fluid_count_comparison_qc(PREVIOUS_FLUID_CELL_COUNT + 1_200)
    assert passed["status"] == "PASS"
    assert failed["status"] == "FAIL"
    assert passed["maximum_relative_difference"] == 0.005


def test_uniform_fluid_connectivity_counts_regions() -> None:
    # At level 2, tree IDs are first_id(2)=9 plus Morton codes.
    connected = np.asarray((9, 10, 11), dtype=np.int64)
    disconnected = np.asarray((9, 10, 72), dtype=np.int64)
    assert connected_fluid_region_count(connected, 2) == 1
    assert connected_fluid_region_count(disconnected, 2) == 2
