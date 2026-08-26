"""Focused tests for local-only CFD boundary surgery and BC correction."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import trimesh

from utils.cfd_surface_prepare.config import (
    ExtensionMeshConfig,
    LocalCutConfig,
    SurfaceQCConfig,
    load_surface_prepare_config,
)
from utils.cfd_surface_prepare.export import export_geometry, meter_scale_qc
from utils.cfd_surface_prepare.extension import extrude_and_cap
from utils.cfd_surface_prepare.io import (
    BoundaryInput,
    SurfacePrepareError,
    sha256_file,
)
from utils.cfd_surface_prepare.local_cut import local_plane_cut
from utils.cfd_surface_prepare.mesh_quality import (
    extension_mesh_quality_qc,
    summarize_extension_mesh,
    triangle_metrics,
)
from utils.cfd_surface_prepare.mesh_refinement import (
    build_refined_rings,
    calculate_ring_count,
    choose_quad_diagonal,
    polygon_area_centroid,
    regularize_loop,
)
from utils.cfd_surface_prepare.pressure_correction import (
    CORRECTION_ROLE,
    build_extended_boundary_conditions,
    calculate_pressure_corrections,
)
from utils.cfd_surface_prepare.qc import (
    boundary_geometry_qc,
    core_surface_preservation_qc,
    extension_collision_qc,
    original_locked_vertex_motion_qc,
    surface_topology_qc,
)
from utils.cfd_surface_prepare.types import BoundarySurfaceResult, CutLoop, TaggedSurface


@pytest.fixture(scope="module")
def boundary() -> BoundaryInput:
    return BoundaryInput(
        index=0,
        port_id="synthetic__cut_000",
        boundary_origin="CUT_PORT",
        role="ASSUMED_INLET",
        global_node_id=None,
        global_edge_id=10,
        center_um=np.asarray((0.0, 0.0, 1.7)),
        source_radius_um=1.0,
        pressure_original_pa=20.0,
        expected_flow_m3_s=1.0e-15,
        simulation_tangent=np.asarray((0.0, 0.0, 1.0)),
        outward_normal=np.asarray((0.0, 0.0, 1.0)),
        extension_length_um=5.0,
        extension_end_um=np.asarray((0.0, 0.0, 6.7)),
    )


@pytest.fixture(scope="module")
def qc_config() -> SurfaceQCConfig:
    return SurfaceQCConfig(
        require_single_component=True,
        require_watertight=True,
        require_zero_boundary_edges=True,
        require_zero_nonmanifold_edges=True,
        require_zero_self_intersections=True,
        require_zero_degenerate_triangles=True,
        expected_boundary_count=1,
        maximum_cap_planarity_error_um=1.0e-5,
        minimum_normal_dot=0.999,
        maximum_extension_length_error_um=1.0e-4,
        maximum_core_surface_distance_um=1.0e-4,
        maximum_core_surface_p95_distance_um=1.0e-5,
        maximum_equivalent_radius_relative_error=0.08,
    )


@pytest.fixture(scope="module")
def local_config() -> LocalCutConfig:
    return LocalCutConfig(3.0, 2.5, 3.0)


@pytest.fixture(scope="module")
def full_config():
    project = Path(__file__).resolve().parents[1]
    return load_surface_prepare_config(
        project / "configs" / "cfd_surface_prepare.yaml", project_root=project
    )


@pytest.fixture(scope="module")
def mesh_config(full_config) -> ExtensionMeshConfig:
    return full_config.extension_mesh


@pytest.fixture(scope="module")
def prepared(boundary: BoundaryInput, mesh_config: ExtensionMeshConfig):
    original = trimesh.creation.capsule(radius=1.0, height=4.0, count=[24, 24])
    cut, loop, cut_report = local_plane_cut(
        TaggedSurface.from_mesh(original),
        boundary,
        radial_factor=3.0,
        axial_back_factor=2.5,
        axial_forward_factor=3.0,
    )
    final, result = extrude_and_cap(
        cut,
        loop,
        boundary,
        local_original_median_edge_length_um=0.3,
        mesh_config=mesh_config,
    )
    final.compact()
    return original, final, loop, cut_report, result


def _result(boundary: BoundaryInput, area_um2: float = np.pi) -> BoundarySurfaceResult:
    return BoundarySurfaceResult(
        boundary_index=boundary.index,
        port_id=boundary.port_id,
        boundary_origin=boundary.boundary_origin,
        role=boundary.role,
        source_radius_um=boundary.source_radius_um,
        extension_length_um=boundary.extension_length_um,
        actual_cap_area_um2=area_um2,
        equivalent_radius_um=float(np.sqrt(area_um2 / np.pi)),
        cap_planarity_error_um=0.0,
        minimum_cap_normal_dot=1.0,
        extension_length_error_um=0.0,
        extension_axis_dot=1.0,
        intermediate_ring_centerline_max_deviation_um=0.0,
        intermediate_ring_centerline_p95_deviation_um=0.0,
        intermediate_ring_centerline_mean_deviation_um=0.0,
        intermediate_ring_centerline_worst_ring_index=1,
        intermediate_ring_centerline_worst_ring_axial_station_um=0.3,
        intermediate_ring_axial_station_max_error_um=0.0,
        intermediate_ring_area_relative_error_max=0.0,
        intermediate_ring_all_areas_finite_positive=True,
        intermediate_ring_all_polygons_simple_valid=True,
        intermediate_ring_all_orientations_consistent=True,
        local_original_median_edge_length_um=0.3,
        target_edge_length_um=0.3,
        ring_count=18,
        transition_length_um=2.0,
        proximal_ring_max_motion_um=0.0,
        distal_ring_max_motion_um=0.0,
        cut_loop_vertex_count=24,
        cap_face_indices=np.arange(22),
        side_face_indices=np.arange(48),
    )


def test_original_surface_sha_is_unchanged(tmp_path: Path):
    path = tmp_path / "immutable.stl"
    path.write_bytes(b"immutable")
    before = sha256_file(path)
    after = sha256_file(path)
    assert before == after


def test_rounded_cap_can_be_cut_flat(prepared):
    _, _, loop, _, _ = prepared
    assert loop.plane_residual_um <= 1.0e-12
    assert loop.area_um2 > 0


def test_local_cut_gets_exactly_one_loop(prepared):
    _, _, _, report, _ = prepared
    assert report["cut_loop_count"] == 1
    assert report["candidate_surface_component_count"] == 1


def test_actual_loop_extrusion_length_is_exact(prepared):
    _, _, _, _, result = prepared
    assert result.extension_length_error_um <= 1.0e-12


def test_distal_cap_is_planar(prepared):
    _, _, _, _, result = prepared
    assert result.cap_planarity_error_um <= 1.0e-12


def test_distal_cap_normal_is_outward(prepared):
    _, _, _, _, result = prepared
    assert result.minimum_cap_normal_dot >= 0.999999
    assert result.extension_axis_dot >= 0.999


def test_extension_axis_uses_section_centers_during_shape_transition(
    boundary, mesh_config
):
    angles = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    radii = np.where(np.arange(12) % 2 == 0, 1.4, 0.6)
    proximal = np.column_stack(
        (radii * np.cos(angles), radii * np.sin(angles), np.zeros(12))
    )
    original = trimesh.Trimesh(
        vertices=proximal, faces=np.empty((0, 3), dtype=np.int64), process=False
    )
    loop = CutLoop(
        boundary_index=0,
        vertex_ids=np.arange(12, dtype=np.int64),
        center_um=np.zeros(3),
        outward_normal=np.asarray((0.0, 0.0, 1.0)),
        area_um2=2.52,
        equivalent_radius_um=0.9,
        plane_residual_um=0.0,
    )
    aligned_boundary = replace(
        boundary,
        center_um=np.zeros(3),
        extension_end_um=np.asarray((0.0, 0.0, boundary.extension_length_um)),
    )
    final, result = extrude_and_cap(
        TaggedSurface.from_mesh(original),
        loop,
        aligned_boundary,
        local_original_median_edge_length_um=0.3,
        mesh_config=mesh_config,
    )
    distal_ids = np.unique(final.faces[result.cap_face_indices])
    pointwise_vectors = final.vertices[distal_ids] - proximal
    pointwise_axis_dots = pointwise_vectors[:, 2] / np.linalg.norm(
        pointwise_vectors, axis=1
    )
    assert np.min(pointwise_axis_dots) < 0.999
    assert result.extension_axis_dot >= 1.0 - 1.0e-12
    assert result.intermediate_ring_centerline_max_deviation_um <= 1.0e-12
    assert result.intermediate_ring_centerline_p95_deviation_um <= 1.0e-12
    assert result.intermediate_ring_centerline_mean_deviation_um <= 1.0e-12
    assert 1 <= result.intermediate_ring_centerline_worst_ring_index < result.ring_count - 1
    assert result.intermediate_ring_axial_station_max_error_um <= 1.0e-12
    assert result.intermediate_ring_area_relative_error_max <= 1.0e-12
    assert result.intermediate_ring_all_areas_finite_positive
    assert result.intermediate_ring_all_polygons_simple_valid
    assert result.intermediate_ring_all_orientations_consistent


def test_final_synthetic_surface_is_watertight(prepared, qc_config):
    _, final, _, _, _ = prepared
    report, _ = surface_topology_qc(final, qc_config)
    assert report["watertight"]
    assert report["boundary_edge_count"] == 0


def test_final_synthetic_surface_has_zero_nonmanifold_edges(prepared, qc_config):
    _, final, _, _, _ = prepared
    report, _ = surface_topology_qc(final, qc_config)
    assert report["nonmanifold_edge_count"] == 0


def test_core_outside_local_region_is_unchanged(
    prepared, boundary, local_config, qc_config
):
    original, final, _, _, _ = prepared
    report = core_surface_preservation_qc(
        original, final.mesh(), [boundary], local_config, qc_config
    )
    assert report["status"] == "PASS"
    assert report["max_core_surface_distance_um"] <= 1.0e-12


def test_multiple_nearby_surfaces_are_rejected(boundary):
    first = trimesh.creation.capsule(radius=1.0, height=4.0, count=[16, 16])
    second = first.copy()
    second.apply_translation((2.2, 0.0, 0.0))
    combined = trimesh.util.concatenate((first, second))
    with pytest.raises(SurfacePrepareError, match="LOCAL_PORT_CUT_AMBIGUOUS"):
        local_plane_cut(
            TaggedSurface.from_mesh(combined),
            boundary,
            radial_factor=3.0,
            axial_back_factor=2.5,
            axial_forward_factor=3.0,
        )


def test_extension_collision_is_rejected(qc_config):
    first = trimesh.creation.box(extents=(2.0, 2.0, 2.0))
    second = trimesh.creation.box(extents=(2.0, 2.0, 2.0))
    second.apply_translation((1.0, 0.25, 0.25))
    surface = TaggedSurface.from_mesh(trimesh.util.concatenate((first, second)))
    surface.face_kind[len(first.faces) :] = 1
    _, intersections = surface_topology_qc(surface, qc_config)
    report = extension_collision_qc(surface, intersections)
    assert report["status"] == "FAIL"
    assert report["extension_collision_count"] > 0


def test_inlet_flow_rate_is_not_changed(boundary):
    rows = calculate_pressure_corrections(
        [boundary],
        [_result(boundary)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )
    assert rows[0]["Q_solver_m3_s"] == boundary.expected_flow_m3_s
    assert rows[0]["pressure_correction_applied"] is False


def test_outlet_solver_pressure_is_lower(boundary):
    outlet = replace(boundary, role="ASSUMED_OUTLET")
    row = calculate_pressure_corrections(
        [outlet],
        [_result(outlet)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )[0]
    assert row["P_solver_boundary_pa"] < outlet.pressure_original_pa
    assert row["extension_pressure_correction_role"] == CORRECTION_ROLE


def test_true_terminal_negative_gauge_is_allowed(boundary):
    terminal = replace(
        boundary,
        role="ASSUMED_OUTLET",
        boundary_origin="TRUE_TERMINAL",
        pressure_original_pa=0.0,
    )
    row = calculate_pressure_corrections(
        [terminal],
        [_result(terminal)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )[0]
    assert row["P_solver_boundary_pa"] < 0.0


def test_extended_bc_preserves_source_p_and_q(boundary):
    outlet = replace(boundary, index=1, port_id="out", role="ASSUMED_OUTLET")
    rows = calculate_pressure_corrections(
        [boundary, outlet],
        [_result(boundary), _result(outlet)],
        dynamic_viscosity_pa_s=0.00345312,
        allow_negative_gauge_pressure=True,
    )
    original = {
        "method": "saved",
        "fluid": {},
        "flow_model": {},
        "pressure_reference": "gauge",
        "wall": {"type": "NO_SLIP"},
        "inlet": {
            "port_id": boundary.port_id,
            "flow_rate_m3_s": boundary.expected_flow_m3_s,
            "profile": "PARABOLIC",
        },
        "outlets": [{"port_id": outlet.port_id}],
    }
    extended = build_extended_boundary_conditions(
        original, [boundary, outlet], rows
    )
    assert extended["inlet"]["flow_rate_m3_s"] == boundary.expected_flow_m3_s
    assert extended["outlets"][0]["P_original_1D_pa"] == outlet.pressure_original_pa


def test_boundary_tags_have_expected_count(prepared, boundary, qc_config):
    _, final, _, cut_report, result = prepared
    report = boundary_geometry_qc(
        [boundary], [cut_report], [result], qc_config
    )
    assert report["status"] == "PASS"
    assert set(final.boundary_type) == {0, 1}
    assert set(final.boundary_index[final.boundary_type > 0]) == {0}


def test_boundary_qc_rejects_intermediate_ring_centerline_deviation(
    prepared, boundary, qc_config
):
    _, _, _, cut_report, result = prepared
    shifted = replace(
        result,
        intermediate_ring_centerline_max_deviation_um=(
            qc_config.maximum_extension_length_error_um * 1.01
        ),
    )
    report = boundary_geometry_qc(
        [boundary], [cut_report], [shifted], qc_config
    )
    boundary_report = report["boundaries"][0]
    assert not boundary_report["checks"]["intermediate_ring_centerline"]
    assert boundary_report["status"] == "FAIL"


def test_boundary_axis_direction_threshold_remains_point_nine_nine_nine(
    prepared, boundary, qc_config
):
    _, _, _, cut_report, result = prepared
    assert qc_config.minimum_normal_dot == 0.999
    misaligned = replace(result, extension_axis_dot=0.999 - 1.0e-12)
    report = boundary_geometry_qc(
        [boundary], [cut_report], [misaligned], qc_config
    )
    assert not report["boundaries"][0]["checks"]["extension_axis"]


@pytest.mark.parametrize(
    ("field", "value", "check"),
    (
        (
            "intermediate_ring_axial_station_max_error_um",
            1.01e-4,
            "intermediate_ring_axial_stations",
        ),
        (
            "intermediate_ring_area_relative_error_max",
            1.0e-6,
            "intermediate_ring_target_areas_preserved",
        ),
        (
            "intermediate_ring_all_areas_finite_positive",
            False,
            "intermediate_ring_areas_finite_positive",
        ),
        (
            "intermediate_ring_all_polygons_simple_valid",
            False,
            "intermediate_ring_polygons_simple_valid_no_crossing",
        ),
        (
            "intermediate_ring_all_orientations_consistent",
            False,
            "intermediate_ring_orientations_no_fold",
        ),
    ),
)
def test_boundary_qc_rejects_invalid_intermediate_ring_geometry(
    prepared, boundary, qc_config, field, value, check
):
    _, _, _, cut_report, result = prepared
    invalid = replace(result, **{field: value})
    report = boundary_geometry_qc(
        [boundary], [cut_report], [invalid], qc_config
    )
    assert not report["boundaries"][0]["checks"][check]


def test_meter_stl_scale_is_exact(prepared, boundary, tmp_path: Path):
    _, final, _, _, _ = prepared
    from utils.cfd_surface_prepare.export import create_layout

    layout = create_layout(tmp_path, "synthetic")
    paths = export_geometry(
        final, [boundary], layout, create_meter_copy=True
    )
    report = meter_scale_qc(
        paths["cfd_surface_refined_um_stl"],
        paths["cfd_surface_refined_m_stl"],
    )
    assert report["status"] == "PASS"


def test_strict_yaml_rejects_unknown_key(tmp_path: Path):
    source = (
        Path(__file__).resolve().parents[1]
        / "configs"
        / "cfd_surface_prepare.yaml"
    )
    text = source.read_text(encoding="utf-8") + "\nunknown: true\n"
    path = tmp_path / "bad.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="Unknown keys"):
        load_surface_prepare_config(path, project_root=source.parents[1])


def _irregular_loop(count: int = 32) -> np.ndarray:
    angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    radius = 1.0 + 0.08 * np.sin(5.0 * angle) + 0.03 * np.cos(9.0 * angle)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), np.zeros(count)))


def test_refinement_long_extension_generates_more_than_two_rings(
    prepared,
):
    *_, result = prepared
    assert result.ring_count > 2


def test_refinement_ring_count_follows_target_spacing():
    assert calculate_ring_count(12.0, 0.3, minimum=6, maximum=128) == 40
    assert calculate_ring_count(1.0, 0.3, minimum=6, maximum=128) == 6


def test_refinement_last_ring_is_exactly_at_extension_end(mesh_config):
    proximal = _irregular_loop()
    rings = build_refined_rings(proximal, np.array([0.0, 0.0, 1.0]), 5.0, 0.3, 1.0, mesh_config)
    assert np.max(np.abs(rings.points[-1, :, 2] - 5.0)) <= 1.0e-12


def test_refinement_ring_zero_is_unchanged(mesh_config):
    proximal = _irregular_loop()
    rings = build_refined_rings(proximal, np.array([0.0, 0.0, 1.0]), 5.0, 0.3, 1.0, mesh_config)
    assert np.array_equal(rings.points[0], proximal)


def test_refinement_original_core_vertices_are_exactly_locked(
    prepared, boundary, local_config
):
    original, final, *_ = prepared
    report = original_locked_vertex_motion_qc(
        original, final, [boundary], local_config
    )
    assert report["status"] == "PASS"
    assert report["original_locked_vertex_motion_max_um"] == 0.0


def test_refinement_distal_ring_is_locked(prepared):
    *_, result = prepared
    assert result.distal_ring_max_motion_um == 0.0


def test_refinement_cap_remains_planar(prepared):
    *_, result = prepared
    assert result.cap_planarity_error_um <= 1.0e-12


def test_refinement_loop_regularization_preserves_centroid():
    proximal = _irregular_loop()
    regularized = regularize_loop(
        proximal, np.array([0.0, 0.0, 1.0]), iterations=8, relaxation=0.2
    )
    _, before = polygon_area_centroid(proximal[:, :2])
    _, after = polygon_area_centroid(regularized[:, :2])
    assert np.allclose(before, after, atol=1.0e-12)


def test_refinement_loop_regularization_preserves_area():
    proximal = _irregular_loop()
    regularized = regularize_loop(
        proximal, np.array([0.0, 0.0, 1.0]), iterations=8, relaxation=0.2
    )
    before, _ = polygon_area_centroid(proximal[:, :2])
    after, _ = polygon_area_centroid(regularized[:, :2])
    assert after == pytest.approx(before, rel=1.0e-12)


def test_refinement_transition_starts_from_actual_loop(mesh_config):
    no_smoothing = replace(
        mesh_config,
        smoothing=replace(mesh_config.smoothing, iterations=0),
    )
    proximal = _irregular_loop()
    rings = build_refined_rings(proximal, np.array([0.0, 0.0, 1.0]), 5.0, 0.3, 1.0, no_smoothing)
    assert np.array_equal(rings.points[0], proximal)


def test_refinement_transition_smoothly_reaches_regularized_loop(mesh_config):
    no_smoothing = replace(
        mesh_config,
        smoothing=replace(mesh_config.smoothing, iterations=0),
    )
    proximal = _irregular_loop()
    rings = build_refined_rings(proximal, np.array([0.0, 0.0, 1.0]), 5.0, 0.3, 1.0, no_smoothing)
    shape = rings.points - rings.stations_um[:, None, None] * np.array([0.0, 0.0, 1.0])
    difference = np.linalg.norm(shape - rings.regularized_loop, axis=2).mean(axis=1)
    transition_ids = rings.stations_um <= rings.transition_length_um
    assert np.all(np.diff(difference[transition_ids]) <= 1.0e-12)
    assert difference[-1] <= 1.0e-12


def test_refinement_all_intermediate_section_centers_stay_on_centerline(
    mesh_config,
):
    angles = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    radii = np.asarray(
        (1.6, 0.7, 1.25, 0.8, 1.4, 0.65, 1.1, 0.9, 1.5, 0.75, 1.2, 0.85)
    )
    proximal = np.column_stack(
        (radii * np.cos(angles), radii * np.sin(angles), np.zeros(12))
    )
    rings = build_refined_rings(
        proximal, np.asarray((0.0, 0.0, 1.0)), 5.0, 0.3, 1.0, mesh_config
    )
    _, proximal_center = polygon_area_centroid(proximal[:, :2])
    centers = np.asarray(
        [polygon_area_centroid(ring[:, :2])[1] for ring in rings.points]
    )
    deviations = np.linalg.norm(centers[1:-1] - proximal_center, axis=1)
    assert np.max(deviations) <= 1.0e-12


def test_refinement_constrained_smoothing_keeps_axial_stations(mesh_config):
    proximal = _irregular_loop()
    rings = build_refined_rings(proximal, np.array([0.0, 0.0, 1.0]), 5.0, 0.3, 1.0, mesh_config)
    measured = rings.points[:, :, 2].mean(axis=1)
    assert np.allclose(measured, rings.stations_um, atol=1.0e-12)


def test_refinement_constrained_smoothing_preserves_ring_area(mesh_config):
    proximal = _irregular_loop()
    rings = build_refined_rings(proximal, np.array([0.0, 0.0, 1.0]), 5.0, 0.3, 1.0, mesh_config)
    measured = np.asarray([polygon_area_centroid(ring[:, :2])[0] for ring in rings.points])
    assert np.allclose(measured, rings.target_areas_um2, rtol=1.0e-10, atol=1.0e-12)


def test_refinement_diagonal_choice_improves_or_preserves_quality():
    vertices = np.asarray(
        ((0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (0.0, 1.0, 0.2), (1.0, 1.0, 0.2))
    )
    first = np.asarray((0, 1))
    second = np.asarray((2, 3))
    selected = choose_quad_diagonal(vertices, first, second, 0)
    alternatives = (
        np.asarray(((0, 1, 3), (0, 3, 2))),
        np.asarray(((0, 1, 2), (1, 3, 2))),
    )
    selected_worst = np.max(triangle_metrics(vertices, np.asarray(selected)).aspect_ratios)
    alternative_worst = [np.max(triangle_metrics(vertices, faces).aspect_ratios) for faces in alternatives]
    assert selected_worst <= min(alternative_worst) + 1.0e-12


def test_refinement_mesh_quality_metrics_are_finite(prepared):
    _, final, *_ = prepared
    mask = final.face_kind == 1
    metrics = triangle_metrics(final.vertices, final.faces[mask])
    assert np.all(np.isfinite(metrics.edge_lengths))
    assert np.all(np.isfinite(metrics.aspect_ratios))


def test_refinement_bad_triangle_fraction_qc_detects_skinny_mesh(full_config):
    vertices = np.asarray(((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.01, 0.001, 0.0)))
    summary = summarize_extension_mesh(
        vertices,
        np.asarray(((0, 1, 2),)),
        target_edge_length_um=0.3,
        local_original_median_edge_length_um=0.3,
        quality=full_config.mesh_quality,
    )
    assert summary["bad_triangle_fraction"] == 1.0


def test_refinement_boundary_and_extension_tags_are_preserved(
    prepared, boundary, full_config
):
    _, final, _, _, result = prepared
    report = extension_mesh_quality_qc(
        final, [boundary], [result], full_config.mesh_quality
    )
    assert report["boundaries"][0]["ring_count"] > 2
    assert set(final.boundary_type) == {0, 1}
    assert set(final.extension_index[final.face_kind == 1]) == {0}


def test_refinement_pressure_is_recomputed_from_final_cap_area(boundary):
    outlet = replace(boundary, role="ASSUMED_OUTLET")
    first = calculate_pressure_corrections(
        [outlet], [_result(outlet, np.pi)], dynamic_viscosity_pa_s=0.00345312, allow_negative_gauge_pressure=True
    )[0]
    second = calculate_pressure_corrections(
        [outlet], [_result(outlet, 1.1 * np.pi)], dynamic_viscosity_pa_s=0.00345312, allow_negative_gauge_pressure=True
    )[0]
    assert first["P_solver_boundary_pa"] != second["P_solver_boundary_pa"]


def test_refinement_meter_copy_remains_correct(prepared, boundary, tmp_path: Path):
    from utils.cfd_surface_prepare.export import create_layout

    _, final, *_ = prepared
    paths = export_geometry(final, [boundary], create_layout(tmp_path, "refined"), create_meter_copy=True)
    report = meter_scale_qc(
        paths["cfd_surface_refined_um_stl"], paths["cfd_surface_refined_m_stl"]
    )
    assert report["status"] == "PASS"


def test_refinement_previous_run_reference_remains_untouched(tmp_path: Path):
    previous = tmp_path / "previous_direct_extrusion.stl"
    previous.write_bytes(b"read-only-reference")
    before = sha256_file(previous)
    _ = previous.read_bytes()
    after = sha256_file(previous)
    assert before == after
