from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import trimesh

from utils.cfd_lumen.branch_local_qc import (
    _BranchDistanceIndex,
    evaluate_branch_local_cross_section_qc,
)
from utils.cfd_lumen.collision_qc import detect_nonadjacent_collisions
from utils.cfd_lumen.config import CFDLumenConfig
from utils.cfd_lumen.context_domain import _DomainBuilder, _core_signature
from utils.cfd_lumen.export import create_run_layout, create_roi_layout, write_geometry_exports
from utils.cfd_lumen.geometry_preprocess import resample_branches, validate_and_extract_branches
from utils.cfd_lumen.hybrid_qc import evaluate_hybrid_surface_qc
from utils.cfd_lumen.hybrid_validation import (
    run_hybrid_synthetic_controls,
    run_v5_synthetic_controls,
)
from utils.cfd_lumen.lumen_builder import build_lumen_surface
from utils.cfd_lumen.mesh_defects import diagnose_mesh_defects
from utils.cfd_lumen.port_geometry import construct_port_geometry
from utils.cfd_lumen.surface_qc import (
    evaluate_radius_fidelity,
    evaluate_surface_qc,
    identify_port_patches,
)
from utils.cfd_lumen.synthetic_controls import run_synthetic_controls
from utils.cfd_lumen.types import BranchGeometry
from utils.cfd_lumen.unified_polyball import (
    JunctionBlendSpec,
    _candidate_polyball_values,
    _stable_smooth_min_reduce,
    build_unified_polyball_surface,
    prepare_polyball_raster,
    release_prepared_polyball_raster,
)
from utils.cfd_lumen.v6_validation import run_v6_synthetic_controls
from utils.sampling.sampling_types import CutPort, GlobalEdge, GlobalVascularModel, ROIRecord


def _config(*, tube_sides: int = 24) -> CFDLumenConfig:
    config = CFDLumenConfig()
    config.geometry.tube_sides = tube_sides
    config.geometry.min_resample_spacing_um = 0.5
    config.geometry.max_resample_spacing_um = 1.0
    config.boolean.allow_implicit_fallback = False
    config.convergence.enabled = False
    config.output.visualizations = False
    return config


def _roi(
    roi_id: str,
    positions: np.ndarray,
    radii: np.ndarray,
    edges: np.ndarray,
    cut_nodes_and_faces: tuple[tuple[int, str], ...],
) -> ROIRecord:
    positions = np.asarray(positions, dtype=float)
    radii = np.asarray(radii, dtype=float)
    edges = np.asarray(edges, dtype=np.int64).reshape((-1, 2))
    minimum = positions.min(axis=0)
    maximum = positions.max(axis=0)
    span = maximum - minimum
    minimum = minimum - np.where(span == 0, 5.0, 0.0)
    maximum = maximum + np.where(span == 0, 5.0, 0.0)
    ports = tuple(
        CutPort(
            cut_port_id=f"{roi_id}__cut_{index:03d}",
            local_node_id=node,
            global_edge_id=100 + next(
                edge_index for edge_index, edge in enumerate(edges) if node in edge
            ),
            intersection_position_um=tuple(map(float, positions[node])),
            radius_at_cut_um=float(radii[node]),
            boundary_face=face,
        )
        for index, (node, face) in enumerate(cut_nodes_and_faces)
    )
    edge_points = positions[edges]
    edge_radii = radii[edges]
    lengths = np.linalg.norm(edge_points[:, 1] - edge_points[:, 0], axis=1)
    return ROIRecord(
        roi_id=roi_id,
        source_model_id="synthetic",
        source_mouse_id="synthetic",
        anchor_id=0,
        anchor_position_um=tuple(map(float, positions[0])),
        bbox_min_um=tuple(map(float, minimum)),
        bbox_max_um=tuple(map(float, maximum)),
        bbox_center_um=tuple(map(float, (minimum + maximum) * 0.5)),
        bbox_size_um=tuple(map(float, maximum - minimum)),
        global_node_ids=tuple(range(1000, 1000 + len(positions))),
        global_edge_ids=tuple(range(100, 100 + len(edges))),
        local_node_ids=np.arange(len(positions), dtype=np.int64),
        local_node_global_ids=np.arange(1000, 1000 + len(positions), dtype=np.int64),
        local_node_positions_um=positions,
        local_node_radius_um=radii,
        local_edges=edges,
        local_edge_ids=np.arange(len(edges), dtype=np.int64),
        local_edge_global_ids=np.arange(100, 100 + len(edges), dtype=np.int64),
        local_edge_points_um=edge_points,
        local_edge_radius_um=edge_radii,
        true_terminal_local_ids=(),
        true_terminal_global_ids=(),
        cut_ports=ports,
        raw_component_count=1,
        raw_total_vessel_length_um=float(lengths.sum()),
        retained_component_length_um=float(lengths.sum()),
    )


def _build(roi: ROIRecord, config: CFDLumenConfig):
    branches, pre_qc = validate_and_extract_branches(roi, config)
    resample_branches(branches, config)
    ports = construct_port_geometry(roi, config)
    build = build_lumen_surface(branches, roi, ports, config)
    patch = identify_port_patches(build.mesh, ports, config)
    surface = evaluate_surface_qc(build.mesh, patch, roi, branches, config)
    samples, fidelity = evaluate_radius_fidelity(build.mesh, branches, roi, config)
    return branches, ports, build, patch, surface, samples, fidelity, pre_qc


def test_straight_constant_radius_vessel_is_watertight_and_faithful() -> None:
    roi = _roi(
        "straight",
        np.asarray(((0, 0, 0), (20, 0, 0), (40, 0, 0))),
        np.full(3, 2.0),
        np.asarray(((0, 1), (1, 2))),
        ((0, "x_min"), (2, "x_max")),
    )
    _, ports, _, patch, surface, samples, fidelity, pre_qc = _build(roi, _config())
    assert pre_qc["branch_count"] == 1
    assert len(ports) == patch.detected_port_count == 2
    assert surface["status"] == "PASS"
    assert surface["checks"]["watertight_trimesh"]
    assert surface["surface_component_count"] == 1
    assert samples
    assert fidelity["p95_absolute_relative_error"] < 0.03
    assert max(row["area_relative_error"] for row in patch.port_rows) < 0.03


def test_tapered_vessel_preserves_radius_change() -> None:
    roi = _roi(
        "tapered",
        np.asarray(((0, 0, 0), (15, 0, 0), (30, 0, 0), (45, 0, 0))),
        np.asarray((1.2, 1.8, 2.4, 3.0)),
        np.asarray(((0, 1), (1, 2), (2, 3))),
        ((0, "x_min"), (3, "x_max")),
    )
    _, _, _, _, surface, samples, fidelity, _ = _build(roi, _config(tube_sides=32))
    assert surface["status"] == "PASS"
    assert len(samples) >= 2
    assert np.corrcoef(
        [sample.source_radius_um for sample in samples],
        [sample.reconstructed_radius_um for sample in samples],
    )[0, 1] > 0.98
    assert fidelity["p95_absolute_relative_error"] < 0.08


def test_curved_vessel_has_continuous_manifold_surface() -> None:
    angles = np.linspace(0.0, np.pi / 2.0, 7)
    positions = np.column_stack((20.0 * np.cos(angles), 20.0 * np.sin(angles), np.zeros_like(angles)))
    roi = _roi(
        "curved",
        positions,
        np.full(len(positions), 1.5),
        np.column_stack((np.arange(len(positions) - 1), np.arange(1, len(positions)))),
        ((0, "x_max"), (len(positions) - 1, "y_max")),
    )
    _, _, _, _, surface, _, _, _ = _build(roi, _config())
    assert surface["status"] == "PASS"
    assert surface["checks"]["manifold"]
    assert surface["checks"]["winding_consistent"]


def test_y_bifurcation_unions_three_branches_and_one_junction() -> None:
    roi = _roi(
        "y_branch",
        np.asarray(((0, 0, 0), (15, 0, 0), (30, 10, 0), (30, -10, 0))),
        np.asarray((2.2, 2.0, 1.5, 1.5)),
        np.asarray(((0, 1), (1, 2), (1, 3))),
        ((0, "x_min"), (2, "y_max"), (3, "y_min")),
    )
    branches, ports, build, patch, surface, _, _, pre_qc = _build(roi, _config(tube_sides=32))
    assert len(branches) == pre_qc["branch_count"] == 3
    assert build.junction_primitive_count == 1
    assert len(ports) == patch.detected_port_count == 3
    assert surface["status"] == "PASS"
    assert surface["surface_component_count"] == 1


def test_nonadjacent_collision_is_rejected_before_boolean() -> None:
    first = BranchGeometry(
        0, (0, 1), (10, 11), (100,),
        np.asarray(((-5, 0, 0), (5, 0, 0)), dtype=float), np.asarray((1.0, 1.0)),
        np.asarray(((-5, 0, 0), (5, 0, 0)), dtype=float), np.asarray((1.0, 1.0)), np.asarray((0.0, 10.0)),
    )
    second = BranchGeometry(
        1, (2, 3), (12, 13), (101,),
        np.asarray(((0, -5, 0), (0, 5, 0)), dtype=float), np.asarray((1.0, 1.0)),
        np.asarray(((0, -5, 0), (0, 5, 0)), dtype=float), np.asarray((1.0, 1.0)), np.asarray((0.0, 10.0)),
    )
    events, report = detect_nonadjacent_collisions([first, second], _config())
    assert report["status"] == "FAIL"
    assert report["hard_collision_count"] == 1
    assert events[0].classification == "HARD_COLLISION"
    assert events[0].clearance_um == pytest.approx(-2.0)


def test_multiple_cut_ports_have_unique_ids_and_outward_normals() -> None:
    roi = _roi(
        "multi_port",
        np.asarray(((0, 0, 0), (15, 0, 0), (30, 10, 0), (30, -10, 0))),
        np.asarray((2.0, 2.0, 1.5, 1.5)),
        np.asarray(((0, 1), (1, 2), (1, 3))),
        ((0, "x_min"), (2, "y_max"), (3, "y_min")),
    )
    _, ports, _, patch, surface, _, _, _ = _build(roi, _config(tube_sides=32))
    assert [port.port_id for port in ports] == [0, 1, 2]
    assert len(set(patch.patch_id[patch.patch_id > 0])) == 3
    assert all(row["normal_alignment"] > 0.9 for row in patch.port_rows)
    assert surface["detected_port_patch_count"] == 3


def test_um_to_m_export_scale_is_exact(tmp_path: Path) -> None:
    roi = _roi(
        "unit_scale",
        np.asarray(((0, 0, 0), (20, 0, 0), (40, 0, 0))),
        np.full(3, 2.0),
        np.asarray(((0, 1), (1, 2))),
        ((0, "x_min"), (2, "x_max")),
    )
    _, ports, build, patch, _, _, _, _ = _build(roi, _config())
    layout = create_run_layout(tmp_path, "units")
    directories = create_roi_layout(layout, roi.roi_id)
    write_geometry_exports(build.mesh, patch, ports, directories)
    exported = trimesh.load_mesh(directories["geometry"] / "lumen_surface_m.stl", process=False)
    assert np.asarray(exported.bounds) == pytest.approx(np.asarray(build.mesh.bounds) * 1.0e-6)


def test_manifold_reconstruction_is_deterministic() -> None:
    roi = _roi(
        "deterministic",
        np.asarray(((0, 0, 0), (20, 0, 0), (40, 0, 0))),
        np.full(3, 2.0),
        np.asarray(((0, 1), (1, 2))),
        ((0, "x_min"), (2, "x_max")),
    )
    first = _build(roi, _config())
    second = _build(roi, _config())
    first_build, first_patch, first_surface = first[2], first[3], first[4]
    second_build, second_patch, second_surface = second[2], second[3], second[4]
    assert first_surface["triangle_count"] == second_surface["triangle_count"]
    assert first_surface["enclosed_volume_um3"] == pytest.approx(second_surface["enclosed_volume_um3"], rel=0, abs=1e-10)
    assert first_patch.patch_id.tolist() == second_patch.patch_id.tolist()
    assert first_build.backend_used == second_build.backend_used == "hybrid_loop_stitch"


def test_legacy_sphere_backend_remains_available_only_when_requested() -> None:
    roi = _roi(
        "legacy_y",
        np.asarray(((0, 0, 0), (15, 0, 0), (30, 10, 0), (30, -10, 0))),
        np.asarray((2.2, 2.0, 1.5, 1.5)),
        np.asarray(((0, 1), (1, 2), (1, 3))),
        ((0, "x_min"), (2, "y_max"), (3, "y_min")),
    )
    config = _config(tube_sides=32)
    config.junction.backend = "legacy_sphere"
    config.reconstruction.junction_backend = "legacy_sphere"
    build = _build(roi, config)[2]
    assert build.backend_used == "manifold"
    assert build.hybrid_details is None


def test_implicit_backend_preserves_flat_detectable_ports() -> None:
    roi = _roi(
        "implicit",
        np.asarray(((0, 0, 0), (10, 0, 0), (20, 0, 0))),
        np.full(3, 2.0),
        np.asarray(((0, 1), (1, 2))),
        ((0, "x_min"), (2, "x_max")),
    )
    config = _config()
    config.boolean.backend = "implicit"
    config.implicit_fallback.min_spacing_um = 0.4
    config.implicit_fallback.max_spacing_um = 0.4
    branches, _ = validate_and_extract_branches(roi, config)
    resample_branches(branches, config)
    ports = construct_port_geometry(roi, config)
    build = build_lumen_surface(branches, roi, ports, config)
    patch = identify_port_patches(
        build.mesh,
        ports,
        config,
        plane_tolerance_um=0.55 * float(build.implicit_grid["spacing_um"]),
        face_normal_alignment=0.0,
    )
    surface = evaluate_surface_qc(build.mesh, patch, roi, branches, config)
    assert build.backend_used == "implicit"
    assert surface["status"] == "PASS"
    assert patch.detected_port_count == 2


def test_v2_synthetic_controls_pass() -> None:
    result = run_synthetic_controls(_config(tube_sides=32))
    assert result["status"] == "PASS"
    controls = {row["name"]: row for row in result["controls"]}
    assert controls["vtk_absolute_radius_scalar_semantics"]["measured_mesh_radius_um"] == pytest.approx(
        2.0, rel=0.01
    )
    assert controls["straight_tube_straight_extension"]["step_radius_relative_error"] < 0.01
    assert controls["simple_symmetric_y_bifurcation"]["self_intersection_count"] == 0


def test_v3_hybrid_synthetic_controls_pass() -> None:
    result = run_hybrid_synthetic_controls(_config(tube_sides=24))
    assert result["status"] == "PASS"
    assert len(result["cases"]) == 6
    assert all(row["self_intersections"] == 0 for row in result["cases"])
    assert all(row["internal_faces"] == 0 for row in result["cases"])
    assert all(row["internal_caps"] == 0 for row in result["cases"])


def test_v5_surface_continuity_synthetic_controls_pass() -> None:
    result = run_v5_synthetic_controls(_config(tube_sides=24))
    assert result["status"] == "PASS"
    assert result["case_count"] == 8
    cases = {row["name"]: row for row in result["cases"]}
    port = cases["curved_tube_to_straight_continuous_extension"]
    assert port["separate_cylinder_primitive_count"] == 0
    assert port["derived_extension_point_count"] > 0
    curved = cases["curved_daughter_loop_stitch_transition"]
    assert curved["transition_backend"] == "loop_stitch"
    assert curved["fallback_reason"] is None
    assert curved["self_intersections"] == 0
    assert curved["internal_faces"] == 0
    assert curved["internal_caps"] == 0


def test_v6_continuous_field_synthetic_controls_pass() -> None:
    result = run_v6_synthetic_controls(_config(tube_sides=24))
    assert result["status"] == "PASS"
    assert result["case_count"] == 6
    cases = {row["name"]: row for row in result["cases"]}
    taper = cases["straight_branch_with_source_taper"]
    assert taper["source_points_unchanged"]
    assert taper["source_radius_unchanged"]
    assert taper["v6_radius_slope_jump"] < taper["v5_radius_slope_jump"]
    curved_port = cases["curved_branch_to_cfd_extension"]
    assert curved_port["v6_tangent_jump_deg"] < curved_port["v5_tangent_jump_deg"]
    for name in (
        "y_junction_continuous_implicit_transition",
        "acute_y_junction",
        "unequal_radius_y_junction",
        "curved_daughter_branches",
    ):
        row = cases[name]
        assert row["transition_backend"] == "continuous_implicit_field"
        assert row["self_intersections"] == 0
        assert row["internal_faces"] == 0
        assert row["internal_caps"] == 0
        assert row["boundary_edges"] == 0
        assert row["nonmanifold_edges"] == 0
        assert row["components"] == 1


def test_v7_exact_polyballline_uses_continuous_radius_interpolation() -> None:
    query = np.asarray(((5.0, 1.5, 0.0), (5.0, 2.0, 0.0)))
    candidates = np.zeros((2, 1), dtype=np.int64)
    phi, gradient = _candidate_polyball_values(
        query,
        candidates,
        np.asarray(((0.0, 0.0, 0.0),)),
        np.asarray(((10.0, 0.0, 0.0),)),
        np.asarray((1.0,)),
        np.asarray((2.0,)),
        gradients=True,
    )
    t = np.linspace(0.0, 1.0, 1_000_001)
    expected = [
        np.min(np.sqrt((5.0 - 10.0 * t) ** 2 + radial**2) - (1.0 + t))
        for radial in (1.5, 2.0)
    ]
    assert phi == pytest.approx(expected, abs=2.0e-10)
    assert gradient is not None
    assert np.linalg.norm(gradient[0]) == pytest.approx(1.0, abs=1.0e-12)
    assert gradient[0, 1] > 0.99


def test_v7_unified_y_has_one_field_no_interface_and_exact_flat_ports() -> None:
    roi = _roi(
        "v7_y",
        np.asarray(((0, 0, 0), (10, 0, 0), (20, 8, 0), (20, -8, 0))),
        np.asarray((2.2, 2.0, 1.5, 1.5)),
        np.asarray(((0, 1), (1, 2), (1, 3))),
        ((0, "x_min"), (2, "y_max"), (3, "y_min")),
    )
    config = _config(tube_sides=24)
    config.v7.cells_across_min_diameter = 12
    config.v7.compare_marching_cubes = True
    config.v7.remesh_minimum_clusters = 100
    branches, _ = validate_and_extract_branches(roi, config)
    branches = resample_branches(branches, config)
    ports = construct_port_geometry(roi, config)
    build = build_unified_polyball_surface(branches, ports, config)
    metadata = build.metadata
    assert build.mesh.is_watertight
    assert len(build.mesh.split(only_watertight=False)) == 1
    assert build.patch.all_ports_pass
    assert metadata["single_continuous_implicit_field"]
    assert metadata["one_iso_surface_extraction"]
    assert metadata["vtk_tube_filter_surface_count"] == 0
    assert metadata["surface_stitch_count"] == 0
    assert metadata["surface_boolean_count"] == 0
    assert metadata["collar_surface_count"] == 0
    assert metadata["hybrid_interface_edge_count"] == 0
    assert all(row["tail_passes_final_plane"] for row in metadata["port_tail_rows"])
    assert all(row["exact_plane_clip"] for row in metadata["port_clip_rows"])
    assert {row["backend"] for row in metadata["extractor_comparison"]} == {
        "flying_edges",
        "marching_cubes",
    }
    assert metadata["projection_after_remesh"][
        "post_projection_phi_abs_p95_um"
    ] <= config.v7.projection_tolerance_um
    assert tuple(build.stage_meshes) == (
        "S0_raw_flying_edges",
        "S1_newton_projected",
        "S2_pyacvd_before_second_projection",
        "S3_final_projected_before_port_clip",
    )


def test_v8_polynomial_smin_and_quintic_support_are_exact() -> None:
    values = np.asarray(((0.0, 0.0), (0.0, 1.0), (1.0, 0.0)))
    spatial_k = np.asarray((0.2, 0.2, 0.0))
    reduced, _ = _stable_smooth_min_reduce(
        values,
        gradients=None,
        spatial_k=spatial_k,
        spatial_k_gradient=np.zeros((3, 3), dtype=float),
    )
    assert reduced[0] == pytest.approx(-0.05)
    assert reduced[1] == pytest.approx(0.0)
    assert reduced[2] == pytest.approx(0.0)

    spec = JunctionBlendSpec(
        junction_node_id=1,
        center_world_um=np.zeros(3),
        radius_um=2.0,
        blend_length_um=4.0,
        incident_branch_ids=(0, 1, 2),
    )
    weight, gradient = spec.compact_weight_and_gradient(
        np.asarray(((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (5.0, 0.0, 0.0)))
    )
    assert weight == pytest.approx((1.0, 0.0, 0.0), abs=1.0e-14)
    assert np.linalg.norm(gradient[0]) == pytest.approx(0.0, abs=1.0e-14)
    assert np.linalg.norm(gradient[1]) == pytest.approx(0.0, abs=1.0e-14)
    assert np.linalg.norm(gradient[2]) == pytest.approx(0.0, abs=1.0e-14)


def test_v8_smooth_junction_uses_same_field_for_extraction_and_projection() -> None:
    roi = _roi(
        "v8_y",
        np.asarray(((0, 0, 0), (10, 0, 0), (20, 8, 0), (20, -8, 0))),
        np.asarray((2.2, 2.0, 1.5, 1.5)),
        np.asarray(((0, 1), (1, 2), (1, 3))),
        ((0, "x_min"), (2, "y_max"), (3, "y_min")),
    )
    config = _config(tube_sides=24)
    config.v7.cells_across_min_diameter = 12
    config.v7.compare_marching_cubes = False
    config.v7.remesh_minimum_clusters = 100
    branches, _ = validate_and_extract_branches(roi, config)
    branches = resample_branches(branches, config)
    ports = construct_port_geometry(roi, config)
    spec = JunctionBlendSpec(
        junction_node_id=1,
        center_world_um=np.asarray(roi.local_node_positions_um[1]),
        radius_um=float(roi.local_node_radius_um[1]),
        blend_length_um=5.0,
        incident_branch_ids=tuple(sorted(branch.branch_id for branch in branches)),
    )
    prepared = prepare_polyball_raster(branches, ports, config)
    original_field = np.asarray(prepared.field).copy()
    try:
        build = build_unified_polyball_surface(
            branches,
            ports,
            config,
            compare_extractors=False,
            junction_specs=(spec,),
            smooth_k_radius_ratio=0.2,
            prepared_raster=prepared,
        )
        assert np.array_equal(np.asarray(prepared.field), original_field)
    finally:
        release_prepared_polyball_raster(prepared)
    assert build.mesh.is_watertight
    assert build.patch.all_ports_pass
    assert build.metadata["backend"] == "unified_polyball_smooth_junction"
    assert build.metadata["projection_field"] == "Phi_v8"
    assert build.metadata["same_field_for_flying_edges_and_both_newton_projections"]
    assert build.metadata["projection_after_remesh"][
        "post_projection_phi_abs_p95_um"
    ] <= config.v7.projection_tolerance_um


def test_v4_segment_index_matches_exhaustive_radius_normalized_distance() -> None:
    rng = np.random.default_rng(17)
    points = np.cumsum(rng.normal(size=(121, 3)), axis=0)
    radii = np.linspace(0.5, 2.0, len(points))
    arc = np.concatenate(
        ([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
    )
    branch = BranchGeometry(
        7,
        tuple(range(len(points))),
        tuple(range(len(points))),
        tuple(range(len(points) - 1)),
        points,
        radii,
        points,
        radii,
        arc,
    )
    queries = rng.normal(size=(500, 3)) * 8.0
    indexed, _ = _BranchDistanceIndex(branch).query(queries)
    starts = points[:-1]
    vectors = np.diff(points, axis=0)
    squared_lengths = np.einsum("ij,ij->i", vectors, vectors)
    relative = queries[:, None, :] - starts[None, :, :]
    projection = np.clip(
        np.einsum("qsi,si->qs", relative, vectors) / squared_lengths[None, :],
        0.0,
        1.0,
    )
    closest = starts[None, :, :] + projection[:, :, None] * vectors[None, :, :]
    local_radius = radii[:-1][None, :] + projection * np.diff(radii)[None, :]
    exhaustive = np.min(
        np.linalg.norm(queries[:, None, :] - closest, axis=2) / local_radius,
        axis=1,
    )
    assert indexed == pytest.approx(exhaustive, rel=0.0, abs=1.0e-12)


def test_v4_context_builder_preserves_core_and_marks_new_cfd_port() -> None:
    core = _roi(
        "context",
        np.asarray(((0, 0, 0), (2, 0, 0))),
        np.asarray((1.0, 1.0)),
        np.asarray(((0, 1),)),
        ((1, "x_max"),),
    )
    core.global_node_ids = (0,)
    core.global_edge_ids = (0,)
    core.local_node_global_ids = np.asarray((0, -1), dtype=np.int64)
    core.local_edge_global_ids = np.asarray((0,), dtype=np.int64)
    core.cut_ports = (
        CutPort(
            cut_port_id="context__cut_000",
            local_node_id=1,
            global_edge_id=0,
            intersection_position_um=(2.0, 0.0, 0.0),
            radius_at_cut_um=1.0,
            boundary_face="x_max",
        ),
    )
    positions = np.asarray(((0, 0, 0), (10, 0, 0), (20, 0, 0)), dtype=float)
    global_edges = (
        GlobalEdge(0, 0, 1, positions[0], positions[1], 1.0, 1.0),
        GlobalEdge(1, 1, 2, positions[1], positions[2], 1.0, 1.0),
    )
    model = GlobalVascularModel(
        source_model_id="synthetic",
        source_mouse_id="synthetic",
        node_ids=np.asarray((0, 1, 2), dtype=np.int64),
        node_positions_um=positions,
        node_radius_um=np.ones(3),
        parent_ids=np.asarray((-1, 0, 1), dtype=np.int64),
        edges=global_edges,
        model_bounds_xyz_um=(0.0, 20.0, 0.0, 0.0, 0.0, 0.0),
        node_index_by_id={0: 0, 1: 1, 2: 2},
        incident_edge_ids_by_node={0: (0,), 1: (0, 1), 2: (1,)},
        global_degree_by_node={0: 1, 1: 2, 2: 1},
    )
    signature = _core_signature(core)
    builder = _DomainBuilder(core, model, _config())
    builder.start_extension(0, "synthetic port conflict")
    domain, _, _ = builder.build_roi()
    assert _core_signature(core) == signature
    assert core.cut_ports[0].boundary_role == "CORE_ROI_BOUNDARY"
    assert domain.node_count == core.node_count + 1
    assert domain.cut_ports[0].boundary_role == "CFD_BOUNDARY_PORT"
    assert builder.owner_added_length(0) == pytest.approx(8.0)


def test_v4_controlled_field_and_branch_local_qc_are_available_on_demand() -> None:
    roi = _roi(
        "controlled_y",
        np.asarray(((0, 0, 0), (15, 0, 0), (30, 10, 0), (30, -10, 0))),
        np.asarray((2.2, 2.0, 1.5, 1.5)),
        np.asarray(((0, 1), (1, 2), (1, 3))),
        ((0, "x_min"), (2, "y_max"), (3, "y_min")),
    )
    config = _config(tube_sides=24)
    config.branch_local_qc.cross_section_cells_across_source_diameter = 32
    config.branch_local_qc.maximum_cross_section_grid_points = 32_768
    branches, _ = validate_and_extract_branches(roi, config)
    resample_branches(branches, config)
    ports = construct_port_geometry(roi, config)
    build = build_lumen_surface(
        branches,
        roi,
        ports,
        config,
        controlled_local_implicit=True,
        transition_backend="manifold_boolean",
        continuous_port_extensions=False,
    )
    assert build.backend_used == "hybrid_controlled_local_implicit"
    assert build.hybrid_details is not None
    assert all(
        patch.metadata["controlled_local_implicit"]
        and patch.metadata["junction_core_field"] == "norm(x - x_J) - r_J"
        for patch in build.hybrid_details.patches.values()
    )
    ownership, rows, report = evaluate_branch_local_cross_section_qc(
        build.mesh, build.hybrid_details, branches, config
    )
    hybrid_qc, _ = evaluate_hybrid_surface_qc(
        build.mesh, build.hybrid_details, roi, config
    )
    assert len(ownership.points) == len(build.mesh.vertices)
    assert rows
    assert report["method"] == "BRANCH_LOCAL_CROSS_SECTION_QC"
    assert hybrid_qc["status"] == "PASS"
    assert hybrid_qc["self_intersection_pairs"] == 0
    assert hybrid_qc["internal_face_count"] == 0
    assert hybrid_qc["internal_cap_face_count"] == 0
    assert hybrid_qc["boundary_edge_count"] == 0
    assert hybrid_qc["nonmanifold_edge_count"] == 0
    assert hybrid_qc["surface_component_count"] == 1


def test_v2_mesh_defects_separates_boundary_and_nonmanifold_edges() -> None:
    closed = trimesh.creation.icosphere(subdivisions=2, radius=2.0)
    closed_report, _ = diagnose_mesh_defects(closed, [])
    assert closed_report["boundary_edge_count"] == 0
    assert closed_report["non_manifold_edge_count"] == 0
    assert closed_report["surface_connected_component_count"] == 1

    open_mesh = closed.copy()
    open_mesh.update_faces(np.arange(len(open_mesh.faces)) != 0)
    open_mesh.remove_unreferenced_vertices()
    open_report, _ = diagnose_mesh_defects(open_mesh, [])
    assert open_report["boundary_edge_count"] == 3
    assert open_report["non_manifold_edge_count"] == 0
