"""Integration quality controls for official VMTK-derived CFD surfaces."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import numpy as np
import pyvista as pv
import trimesh
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

from utils.cfd_lumen.ultraliser_qc import _section_polygon, _triangle_intersections

from .config import LocalCutConfig, MeshQualityConfig
from .io import BoundaryInput, SurfacePrepareError
from .local_cut import orthogonal_basis
from .mesh_quality import summarize_extension_mesh, triangle_metrics


@dataclass(frozen=True, slots=True)
class BoundaryLoop:
    point_ids: np.ndarray
    points: np.ndarray
    center_um: np.ndarray
    area_um2: float
    equivalent_radius_um: float
    normal: np.ndarray
    planarity_error_um: float


def _polygon_area_centroid(points: np.ndarray) -> tuple[float, np.ndarray]:
    following = np.roll(points, -1, axis=0)
    cross = points[:, 0] * following[:, 1] - following[:, 0] * points[:, 1]
    area = 0.5 * float(np.sum(cross))
    if abs(area) <= np.finfo(float).eps:
        raise SurfacePrepareError("VMTK_BOUNDARY_HAS_ZERO_AREA")
    centroid = np.sum((points + following) * cross[:, None], axis=0) / (6.0 * area)
    return area, centroid


def polydata_mesh(path: Path) -> tuple[pv.PolyData, trimesh.Trimesh]:
    data = pv.read(path).triangulate()
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    mesh = trimesh.Trimesh(
        vertices=np.asarray(data.points, dtype=float), faces=faces, process=False
    )
    return data, mesh


def _faces(data: pv.PolyData) -> np.ndarray:
    return np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]


def _output_precision_tolerance(points: np.ndarray) -> float:
    values = np.asarray(points)
    spacing = np.abs(np.spacing(values))
    maximum = float(np.max(spacing)) if spacing.size else 0.0
    return max(np.sqrt(3.0) * maximum, np.finfo(float).eps)


def raw_core_exact_copy_qc(
    input_vtp: Path,
    raw_vtp: Path,
    *,
    tag_regions: bool = True,
) -> dict[str, Any]:
    """Prove that VMTK only appends geometry, allowing its float output ULP."""

    source = pv.read(input_vtp).triangulate()
    raw = pv.read(raw_vtp).triangulate()
    source_faces = _faces(source)
    raw_faces = _faces(raw)
    point_count = source.n_points
    cell_count = source.n_cells
    counts_valid = raw.n_points >= point_count and raw.n_cells >= cell_count
    if counts_valid:
        source_points = np.asarray(source.points, dtype=float)
        retained = np.asarray(raw.points[:point_count])
        motion = np.linalg.norm(retained.astype(float) - source_points, axis=1)
        cast_source = np.asarray(source.points).astype(retained.dtype, copy=False)
        exact_after_cast = bool(np.array_equal(retained, cast_source))
        connectivity_changed = int(
            np.count_nonzero(np.any(raw_faces[:cell_count] != source_faces, axis=1))
        )
        tolerance = _output_precision_tolerance(retained)
        maximum = float(np.max(motion)) if len(motion) else float("inf")
        p95 = float(np.percentile(motion, 95)) if len(motion) else float("inf")
    else:
        exact_after_cast = False
        connectivity_changed = max(cell_count, 1)
        tolerance = 0.0
        maximum = float("inf")
        p95 = float("inf")
    checks = {
        "raw_appends_points_and_cells": counts_valid,
        "retained_points_exact_after_output_dtype_cast": exact_after_cast,
        "retained_point_motion_within_output_machine_precision": maximum <= tolerance,
        "original_input_cell_connectivity_unchanged": connectivity_changed == 0,
    }
    passed = all(checks.values())
    if tag_regions and passed:
        region = np.ones(raw.n_cells, dtype=np.uint8)
        region[:cell_count] = 0
        raw.cell_data["SurfaceRegionId"] = region
        raw.cell_data["SurfaceRegion"] = np.where(
            region == 0, "CORE", "EXTENSION"
        )
        raw.save(raw_vtp, binary=True)
    return {
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "input_point_count": int(point_count),
        "raw_point_count": int(raw.n_points),
        "input_cell_count": int(cell_count),
        "raw_cell_count": int(raw.n_cells),
        "retained_input_point_max_motion_um": maximum,
        "retained_input_point_P95_motion_um": p95,
        "original_input_cell_connectivity_changed_count": connectivity_changed,
        "new_extension_point_count": int(raw.n_points - point_count),
        "new_extension_cell_count": int(raw.n_cells - cell_count),
        "input_point_dtype": str(np.asarray(source.points).dtype),
        "raw_point_dtype": str(np.asarray(raw.points).dtype),
        "machine_precision_tolerance_um": tolerance,
        "retained_points_exact_after_output_dtype_cast": exact_after_cast,
        "classification_method": "validated input cell/point prefix ordering",
        "surface_region_codes": {"CORE": 0, "EXTENSION": 1},
    }


def _ordered_component(graph: nx.Graph, component: set[int]) -> np.ndarray:
    if len(component) < 3 or any(graph.degree(node) != 2 for node in component):
        raise SurfacePrepareError("VMTK_BOUNDARY_NOT_SIMPLE_CLOSED_LOOP")
    start = min(component)
    ordered = [start]
    previous = -1
    current = start
    while True:
        candidates = [node for node in graph.neighbors(current) if node != previous]
        if not candidates:
            raise SurfacePrepareError("VMTK_BOUNDARY_NOT_SIMPLE_CLOSED_LOOP")
        following = candidates[0]
        if following == start:
            break
        if following in ordered:
            raise SurfacePrepareError("VMTK_BOUNDARY_NOT_SIMPLE_CLOSED_LOOP")
        ordered.append(following)
        previous, current = current, following
    if len(ordered) != len(component):
        raise SurfacePrepareError("VMTK_BOUNDARY_NOT_SIMPLE_CLOSED_LOOP")
    return np.asarray(ordered, dtype=np.int64)


def extract_boundary_loops(mesh: trimesh.Trimesh) -> tuple[BoundaryLoop, ...]:
    faces = np.asarray(mesh.faces, dtype=np.int64)
    directed = np.vstack((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    sorted_edges = np.sort(directed, axis=1)
    unique, counts = np.unique(sorted_edges, axis=0, return_counts=True)
    boundary_edges = unique[(counts == 1) & (unique[:, 0] != unique[:, 1])]
    if len(boundary_edges) == 0:
        return ()
    graph = nx.Graph()
    graph.add_edges_from((int(a), int(b)) for a, b in boundary_edges)
    loops: list[BoundaryLoop] = []
    vertices = np.asarray(mesh.vertices, dtype=float)
    for component in nx.connected_components(graph):
        ids = _ordered_component(graph, set(component))
        points = vertices[ids]
        centered = points - points.mean(axis=0)
        _, _, right = np.linalg.svd(centered, full_matrices=False)
        normal = right[-1]
        first, second = orthogonal_basis(normal)
        projected = np.column_stack((points @ first, points @ second))
        signed_area, center_2d = _polygon_area_centroid(projected)
        if signed_area < 0:
            ids = ids[::-1]
            points = points[::-1]
            projected = projected[::-1]
            signed_area, center_2d = _polygon_area_centroid(projected)
            normal *= -1.0
        center = center_2d[0] * first + center_2d[1] * second
        plane_offset = float(np.mean(points @ normal))
        center += plane_offset * normal
        residual = float(np.max(np.abs((points - center) @ normal)))
        area = float(abs(signed_area))
        loops.append(
            BoundaryLoop(
                point_ids=ids,
                points=points,
                center_um=center,
                area_um2=area,
                equivalent_radius_um=float(np.sqrt(area / np.pi)),
                normal=normal / np.linalg.norm(normal),
                planarity_error_um=residual,
            )
        )
    return tuple(loops)


def map_loops_to_boundaries(
    loops: Iterable[BoundaryLoop],
    boundaries: Iterable[BoundaryInput],
    *,
    distal: bool,
) -> dict[int, BoundaryLoop]:
    loop_list = list(loops)
    boundary_list = list(boundaries)
    if len(loop_list) != len(boundary_list):
        raise SurfacePrepareError("VMTK_BOUNDARY_COUNT_MISMATCH")
    targets = np.asarray(
        [boundary.extension_end_um if distal else boundary.center_um for boundary in boundary_list]
    )
    centers = np.asarray([loop.center_um for loop in loop_list])
    distances = np.linalg.norm(targets[:, None, :] - centers[None, :, :], axis=2)
    rows, columns = linear_sum_assignment(distances)
    if len(set(map(int, columns))) != len(boundary_list):
        raise SurfacePrepareError("VMTK_BOUNDARY_MAPPING_FAILED")
    return {
        boundary_list[int(row)].index: loop_list[int(column)]
        for row, column in zip(rows, columns, strict=True)
    }


def open_profile_qc(
    mesh: trimesh.Trimesh, boundaries: Iterable[BoundaryInput]
) -> tuple[dict[str, Any], dict[int, BoundaryLoop]]:
    boundary_list = list(boundaries)
    loops = extract_boundary_loops(mesh)
    mapping = map_loops_to_boundaries(loops, boundary_list, distal=False)
    rows: list[dict[str, Any]] = []
    for boundary in boundary_list:
        loop = mapping[boundary.index]
        normal_dot = float(abs(np.dot(loop.normal, boundary.outward_normal)))
        center_distance = float(np.linalg.norm(loop.center_um - boundary.center_um))
        checks = {
            "simple_closed_loop": len(loop.point_ids) >= 3,
            "positive_area": np.isfinite(loop.area_um2) and loop.area_um2 > 0.0,
            "center_matches_preprocess_plane": center_distance
            <= boundary.source_radius_um,
            "normal_matches_preprocess_outward": normal_dot >= 0.999,
        }
        rows.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "role": boundary.role,
                "boundary_origin": boundary.boundary_origin,
                "point_count": len(loop.point_ids),
                "area_um2": loop.area_um2,
                "equivalent_radius_um": loop.equivalent_radius_um,
                "center_um": loop.center_um.tolist(),
                "center_distance_um": center_distance,
                "center_distance_to_expected_boundary_um": center_distance,
                "normal_abs_dot_outward": normal_dot,
                "boundary_plane_normal_abs_dot_expected_outward": normal_dot,
                "planarity_error_um": loop.planarity_error_um,
                "checks": checks,
                "status": "PASS" if all(checks.values()) else "FAIL",
            }
        )
    count_checks = {
        "exactly_four_open_profiles": len(loops) == 4,
        "one_inlet": sum(item.role == "ASSUMED_INLET" for item in boundary_list) == 1,
        "three_outlets": sum(item.role == "ASSUMED_OUTLET" for item in boundary_list) == 3,
    }
    passed = all(count_checks.values()) and all(row["status"] == "PASS" for row in rows)
    return {
        "status": "PASS" if passed else "FAIL",
        "profile_count": len(loops),
        "checks": count_checks,
        "boundaries": rows,
    }, mapping


def topology_qc(
    mesh: trimesh.Trimesh,
    *,
    expected_open_profile_count: int,
    allow_degenerate: bool = False,
    require_winding_consistent: bool = False,
) -> tuple[dict[str, Any], list[tuple[int, int]]]:
    sorted_edges = np.sort(np.asarray(mesh.edges, dtype=np.int64), axis=1)
    _, edge_counts = np.unique(sorted_edges, axis=0, return_counts=True)
    boundary_edges = int(np.count_nonzero(edge_counts == 1))
    nonmanifold = int(np.count_nonzero(edge_counts > 2))
    loops = extract_boundary_loops(mesh)
    diagonal = float(np.linalg.norm(np.ptp(np.asarray(mesh.vertices), axis=0)))
    area_tolerance = max(np.finfo(float).eps * diagonal**2 * 100.0, 1.0e-18)
    areas = np.asarray(mesh.area_faces, dtype=float)
    repeated = np.asarray([len(set(map(int, face))) < 3 for face in mesh.faces])
    degenerate = int(np.count_nonzero((areas <= area_tolerance) | repeated))
    valid_face_mask = (areas > area_tolerance) & ~repeated
    assessment_mesh = mesh
    if allow_degenerate and not np.all(valid_face_mask):
        assessment_mesh = trimesh.Trimesh(
            vertices=np.asarray(mesh.vertices, dtype=float),
            faces=np.asarray(mesh.faces, dtype=np.int64)[valid_face_mask],
            process=False,
        )
    assessment_edges = np.sort(np.asarray(assessment_mesh.edges, dtype=np.int64), axis=1)
    _, assessment_counts = np.unique(assessment_edges, axis=0, return_counts=True)
    assessment_nonmanifold = int(np.count_nonzero(assessment_counts > 2))
    assessment_components = len(assessment_mesh.split(only_watertight=False))
    intersections, candidate_count = _triangle_intersections(
        mesh, np.arange(len(mesh.faces), dtype=np.int64)
    )
    components = len(mesh.split(only_watertight=False))
    checks = {
        "single_component": assessment_components == 1,
        "expected_open_profiles": len(loops) == expected_open_profile_count,
        "watertight_when_closed": expected_open_profile_count != 0 or mesh.is_watertight,
        "zero_nonmanifold_edges": assessment_nonmanifold == 0,
        "zero_self_intersections": len(intersections) == 0,
        "zero_degenerate_triangles": degenerate == 0 or allow_degenerate,
        "winding_consistent": (
            not require_winding_consistent or bool(mesh.is_winding_consistent)
        ),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "vertex_count": len(mesh.vertices),
        "triangle_count": len(mesh.faces),
        "component_count": components,
        "assessment_component_count_excluding_official_raw_degenerates": assessment_components,
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "boundary_edge_count": boundary_edges,
        "open_profile_count": len(loops),
        "nonmanifold_edge_count": nonmanifold,
        "assessment_nonmanifold_edge_count_excluding_official_raw_degenerates": assessment_nonmanifold,
        "self_intersection_count": len(intersections),
        "self_intersection_candidate_pairs_checked": candidate_count,
        "degenerate_triangle_count": degenerate,
        "official_raw_degenerate_triangles_allowed_until_vmtk_remesh": allow_degenerate,
    }, intersections


def _interface_angles(
    mesh: trimesh.Trimesh,
    core_face_mask: np.ndarray,
    extension_face_mask: np.ndarray,
    boundaries: list[BoundaryInput],
) -> dict[int, dict[str, Any]]:
    edge_faces: dict[tuple[int, int], list[int]] = {}
    for face_id, face in enumerate(np.asarray(mesh.faces, dtype=np.int64)):
        for first, second in zip(face, np.roll(face, -1)):
            edge_faces.setdefault(tuple(sorted((int(first), int(second)))), []).append(face_id)
    values: dict[int, list[float]] = {boundary.index: [] for boundary in boundaries}
    edge_count = {boundary.index: 0 for boundary in boundaries}
    vertices = np.asarray(mesh.vertices, dtype=float)
    normals = np.asarray(mesh.face_normals, dtype=float)
    for edge, linked in edge_faces.items():
        if len(linked) != 2:
            continue
        first, second = linked
        if not (
            (core_face_mask[first] and extension_face_mask[second])
            or (core_face_mask[second] and extension_face_mask[first])
        ):
            continue
        midpoint = vertices[np.asarray(edge)].mean(axis=0)
        distances = [np.linalg.norm(midpoint - boundary.center_um) for boundary in boundaries]
        boundary = boundaries[int(np.argmin(distances))]
        dot = float(np.clip(np.dot(normals[first], normals[second]), -1.0, 1.0))
        values[boundary.index].append(float(np.degrees(np.arccos(dot))))
        edge_count[boundary.index] += 1
    output: dict[int, dict[str, Any]] = {}
    for boundary in boundaries:
        angles = np.asarray(values[boundary.index], dtype=float)
        if len(angles) == 0:
            output[boundary.index] = {
                "interface_edge_count": 0,
                "normal_jump_P50_deg": None,
                "normal_jump_P95_deg": None,
                "normal_jump_P99_deg": None,
                "normal_jump_max_deg": None,
            }
        else:
            output[boundary.index] = {
                "interface_edge_count": edge_count[boundary.index],
                "normal_jump_P50_deg": float(np.percentile(angles, 50)),
                "normal_jump_P95_deg": float(np.percentile(angles, 95)),
                "normal_jump_P99_deg": float(np.percentile(angles, 99)),
                "normal_jump_max_deg": float(np.max(angles)),
            }
    return output


def interface_smoothness_from_raw(
    mesh: trimesh.Trimesh,
    *,
    input_point_count: int,
    boundaries: Iterable[BoundaryInput],
) -> dict[str, Any]:
    faces = np.asarray(mesh.faces, dtype=np.int64)
    extension = np.any(faces >= input_point_count, axis=1)
    core = ~extension
    boundary_list = list(boundaries)
    rows = _interface_angles(mesh, core, extension, boundary_list)
    return {
        "method": "normal angle across shared original/VMTK interface edges",
        "boundaries": [
            {"port_id": boundary.port_id, **rows[boundary.index]}
            for boundary in boundary_list
        ],
    }


def interface_smoothness_from_old_custom(
    old_vtp: Path, boundaries: Iterable[BoundaryInput]
) -> dict[str, Any]:
    from .visualization import load_previous_tagged_surface

    boundary_list = list(boundaries)
    tagged = load_previous_tagged_surface(old_vtp, boundary_list)
    mesh = tagged.mesh()
    rows = _interface_angles(
        mesh, tagged.face_kind == 0, tagged.face_kind == 1, boundary_list
    )
    return {
        "method": "normal angle across shared old-custom core/extension edges",
        "boundaries": [
            {"port_id": boundary.port_id, **rows[boundary.index]}
            for boundary in boundary_list
        ],
    }


def extension_vector_measurements(
    proximal_center_um: np.ndarray,
    distal_center_um: np.ndarray,
    expected_outward_normal: np.ndarray,
) -> dict[str, float]:
    """Measure extension direction with a signed outward dot product."""

    vector = np.asarray(distal_center_um) - np.asarray(proximal_center_um)
    norm = float(np.linalg.norm(vector))
    direction_dot = (
        -1.0
        if not np.isfinite(norm) or norm <= np.finfo(float).eps
        else float(np.dot(vector / norm, expected_outward_normal))
    )
    return {
        "actual_center_to_center_norm_um": norm,
        "actual_axial_length_um": float(np.dot(vector, expected_outward_normal)),
        "extension_direction_dot": direction_dot,
    }


def extension_geometry_qc(
    raw_mesh: trimesh.Trimesh,
    boundaries: Iterable[BoundaryInput],
    proximal_loops: dict[int, BoundaryLoop],
    interface_report: dict[str, Any],
) -> tuple[dict[str, Any], dict[int, BoundaryLoop]]:
    boundary_list = list(boundaries)
    distal_loops = extract_boundary_loops(raw_mesh)
    distal_mapping = map_loops_to_boundaries(distal_loops, boundary_list, distal=True)
    interfaces = {row["port_id"]: row for row in interface_report["boundaries"]}
    rows: list[dict[str, Any]] = []
    for boundary in boundary_list:
        proximal = proximal_loops[boundary.index]
        distal = distal_mapping[boundary.index]
        vector_measurements = extension_vector_measurements(
            proximal.center_um, distal.center_um, boundary.outward_normal
        )
        norm = vector_measurements["actual_center_to_center_norm_um"]
        direction_dot = vector_measurements["extension_direction_dot"]
        actual_axial_length = vector_measurements["actual_axial_length_um"]
        length_error = (
            abs(actual_axial_length - boundary.extension_length_um)
            / boundary.extension_length_um
        )
        area_error = abs(distal.area_um2 - proximal.area_um2) / proximal.area_um2
        checks = {
            "direction_dot_at_least_0_999": direction_dot >= 0.999,
            "relative_length_error_at_most_0_02": length_error <= 0.02,
            "distal_area_relative_error_at_most_0_05": area_error <= 0.05,
            "positive_distal_area": distal.area_um2 > 0.0,
        }
        rows.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "role": boundary.role,
                "boundary_origin": boundary.boundary_origin,
                "planned_extension_length_um": boundary.extension_length_um,
                "actual_extension_length_um": actual_axial_length,
                "actual_axial_length_um": actual_axial_length,
                "actual_center_to_center_norm_um": norm,
                "extension_length_relative_error": length_error,
                "extension_direction_dot": direction_dot,
                "original_area_um2": proximal.area_um2,
                "proximal_area_um2": proximal.area_um2,
                "distal_area_um2": distal.area_um2,
                "distal_area_relative_error": area_error,
                "equivalent_radius_original_um": proximal.equivalent_radius_um,
                "equivalent_radius_distal_um": distal.equivalent_radius_um,
                "distal_center_um": distal.center_um.tolist(),
                "distal_planarity_error_um": distal.planarity_error_um,
                **interfaces[boundary.port_id],
                "interface_diagnostic_available": (
                    interfaces[boundary.port_id]["interface_edge_count"] > 0
                ),
                "checks": checks,
                "status": "PASS" if all(checks.values()) else "FAIL",
            }
        )
    passed = len(rows) == 4 and all(row["status"] == "PASS" for row in rows)
    return {
        "status": "PASS" if passed else "FAIL",
        "boundary_count": len(rows),
        "thresholds": {
            "minimum_direction_dot": 0.999,
            "maximum_length_relative_error": 0.02,
            "maximum_distal_area_relative_error": 0.05,
        },
        "boundaries": rows,
    }, distal_mapping


def extension_mesh_metrics(
    mesh: trimesh.Trimesh,
    boundaries: Iterable[BoundaryInput],
    *,
    added_face_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    boundary_list = list(boundaries)
    centers = np.asarray(mesh.triangles_center, dtype=float)
    rows: list[dict[str, Any]] = []
    for boundary in boundary_list:
        relative = centers - boundary.center_um
        axial = relative @ boundary.outward_normal
        radial = np.linalg.norm(relative - np.outer(axial, boundary.outward_normal), axis=1)
        mask = (
            (axial >= -0.1 * boundary.source_radius_um)
            & (axial <= 1.05 * boundary.extension_length_um)
            & (radial <= 2.5 * boundary.source_radius_um)
        )
        if added_face_mask is not None:
            mask &= added_face_mask
        metrics = triangle_metrics(mesh.vertices, np.asarray(mesh.faces)[mask])
        rows.append(
            {
                "port_id": boundary.port_id,
                "triangle_count": int(np.count_nonzero(mask)),
                "minimum_angle_deg": float(np.min(metrics.minimum_angles_deg)),
                "aspect_ratio_P95": float(np.percentile(metrics.aspect_ratios, 95)),
                "aspect_ratio_max": float(np.max(metrics.aspect_ratios)),
                "edge_length_median_um": float(np.median(metrics.edge_lengths)),
            }
        )
    return {"boundaries": rows}


def extension_mesh_quality_from_raw(
    raw_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    *,
    local_target_edge_um: dict[int, float],
    quality: MeshQualityConfig,
    input_point_count: int,
) -> dict[str, Any]:
    """Apply the existing project mesh-quality rules to RAW extension cells only."""

    data, mesh = polydata_mesh(raw_vtp)
    if "SurfaceRegionId" not in data.cell_data:
        raise SurfacePrepareError("VMTK_RAW_CORE_NOT_EXACT_COPY:regions_missing")
    region = np.asarray(data.cell_data["SurfaceRegionId"], dtype=np.uint8)
    extension_ids = np.flatnonzero(region == 1)
    if len(extension_ids) == 0:
        raise SurfacePrepareError("VMTK_RAW_EXTENSION_MESH_QUALITY_FAILED")
    boundary_list = list(boundaries)
    centers = np.asarray(mesh.triangles_center, dtype=float)[extension_ids]
    scores: list[np.ndarray] = []
    for boundary in boundary_list:
        relative = centers - boundary.center_um
        axial = relative @ boundary.outward_normal
        radial = np.linalg.norm(
            relative - np.outer(axial, boundary.outward_normal), axis=1
        )
        before = np.maximum(-axial, 0.0)
        after = np.maximum(axial - boundary.extension_length_um, 0.0)
        scores.append(
            (radial / boundary.source_radius_um) ** 2
            + ((before + after) / boundary.source_radius_um) ** 2
        )
    assignment = np.argmin(np.column_stack(scores), axis=1)
    rows: list[dict[str, Any]] = []
    faces = np.asarray(mesh.faces, dtype=np.int64)
    for local_index, boundary in enumerate(boundary_list):
        selected_ids = extension_ids[assignment == local_index]
        selected_faces = faces[selected_ids]
        target = float(local_target_edge_um[boundary.index])
        summary = summarize_extension_mesh(
            np.asarray(mesh.vertices, dtype=float),
            selected_faces,
            target_edge_length_um=target,
            local_original_median_edge_length_um=target,
            quality=quality,
        )
        interface_mask = np.any(selected_faces < input_point_count, axis=1)
        if np.any(interface_mask):
            interface = triangle_metrics(
                np.asarray(mesh.vertices, dtype=float), selected_faces[interface_mask]
            )
            interface_ratio = float(np.median(interface.edge_lengths) / target)
        else:
            interface_ratio = float("inf")
        finite = all(
            np.isfinite(value)
            for value in summary.values()
            if isinstance(value, float)
        ) and np.isfinite(interface_ratio)
        checks = {
            "finite_metrics": bool(finite),
            "bad_triangle_fraction": (
                summary["bad_triangle_fraction"]
                <= quality.maximum_bad_triangle_fraction
            ),
            "neighbor_area_ratio_p95": (
                summary["neighbor_area_ratio_p95"]
                <= quality.maximum_neighbor_area_ratio
            ),
            "interface_edge_length_ratio": (
                interface_ratio <= quality.maximum_interface_edge_length_ratio
            ),
        }
        rows.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "role": boundary.role,
                "triangle_count": int(len(selected_faces)),
                "minimum_angle_deg": summary["minimum_triangle_angle_deg"],
                "angle_P05_deg": summary["triangle_angle_p05_deg"],
                "aspect_ratio_median": summary["aspect_ratio_median"],
                "aspect_ratio_P95": summary["aspect_ratio_p95"],
                "aspect_ratio_max": summary["aspect_ratio_max"],
                "edge_length_median_um": summary["edge_length_median_um"],
                "edge_length_P95_um": summary["edge_length_p95_um"],
                "bad_triangle_fraction": summary["bad_triangle_fraction"],
                "neighbor_area_ratio_P95": summary["neighbor_area_ratio_p95"],
                "interface_edge_length_ratio": interface_ratio,
                "local_original_target_edge_um": target,
                "checks": checks,
                "status": "PASS" if all(checks.values()) else "FAIL",
            }
        )
    return {
        "status": (
            "PASS"
            if len(rows) == len(boundary_list)
            and all(row["status"] == "PASS" for row in rows)
            else "FAIL"
        ),
        "scope": "SurfaceRegionId == EXTENSION only",
        "thresholds": {
            "minimum_triangle_angle_deg": quality.minimum_triangle_angle_deg,
            "maximum_aspect_ratio": quality.maximum_aspect_ratio,
            "maximum_edge_length_to_local_target_ratio": (
                quality.maximum_edge_length_to_local_target_ratio
            ),
            "maximum_neighbor_area_ratio": quality.maximum_neighbor_area_ratio,
            "maximum_interface_edge_length_ratio": (
                quality.maximum_interface_edge_length_ratio
            ),
            "maximum_bad_triangle_fraction": quality.maximum_bad_triangle_fraction,
        },
        "boundaries": rows,
    }


def tag_and_export_final_surface(
    capped_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    geometry_directory: Path,
    boundary_directory: Path,
    *,
    output_stem: str = "cfd_surface_vmtk_tps_boundarynormal",
    raw_vtp: Path | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    geometry_directory.mkdir(parents=True, exist_ok=True)
    data = pv.read(capped_vtp).triangulate()
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    entity = np.asarray(data.cell_data["CellEntityIds"], dtype=np.int32)
    unique_ids, entity_counts = np.unique(entity, return_counts=True)
    wall_entity_id = int(unique_ids[int(np.argmax(entity_counts))])
    cap_ids = sorted(int(value) for value in unique_ids if int(value) != wall_entity_id)
    boundary_list = list(boundaries)
    if len(cap_ids) != len(boundary_list):
        raise SurfacePrepareError("VMTK_BOUNDARY_COUNT_MISMATCH")
    centers = np.asarray(data.points)[faces].mean(axis=1)
    cap_centers = np.asarray([centers[entity == cap_id].mean(axis=0) for cap_id in cap_ids])
    targets = np.asarray([boundary.extension_end_um for boundary in boundary_list])
    rows, columns = linear_sum_assignment(
        np.linalg.norm(targets[:, None, :] - cap_centers[None, :, :], axis=2)
    )
    cap_to_boundary = {
        cap_ids[int(column)]: boundary_list[int(row)]
        for row, column in zip(rows, columns, strict=True)
    }
    boundary_type = np.zeros(len(faces), dtype=np.uint8)
    boundary_index = np.full(len(faces), -1, dtype=np.int32)
    origin_code = np.zeros(len(faces), dtype=np.uint8)
    origin = np.full(len(faces), "WALL", dtype="<U16")
    port_width = max(len(boundary.port_id) for boundary in boundary_list)
    port_id = np.full(len(faces), "", dtype=f"<U{port_width}")
    surface_region_id = np.full(len(faces), 2, dtype=np.uint8)
    preserved_raw_triangle_count = 0
    missing_raw_triangle_count = 0
    if raw_vtp is not None:
        raw = pv.read(raw_vtp).triangulate()
        if "SurfaceRegionId" not in raw.cell_data:
            raise SurfacePrepareError("VMTK_RAW_CORE_NOT_EXACT_COPY:regions_missing")
        raw_faces = _faces(raw)
        raw_regions = np.asarray(raw.cell_data["SurfaceRegionId"], dtype=np.uint8)

        def key(points: np.ndarray, face: np.ndarray) -> tuple[tuple[float, ...], ...]:
            return tuple(
                sorted(tuple(float(value) for value in point) for point in points[face])
            )

        available: dict[tuple[tuple[float, ...], ...], list[int]] = defaultdict(list)
        raw_points = np.asarray(raw.points)
        for raw_id, raw_face in enumerate(raw_faces):
            available[key(raw_points, raw_face)].append(int(raw_regions[raw_id]))
        capped_points = np.asarray(data.points)
        for face_id, face in enumerate(faces):
            values = available.get(key(capped_points, face))
            if values:
                surface_region_id[face_id] = values.pop()
                preserved_raw_triangle_count += 1
        missing_raw_triangle_count = sum(len(values) for values in available.values())
        if missing_raw_triangle_count or preserved_raw_triangle_count != len(raw_faces):
            raise SurfacePrepareError("VMTK_CAPONLY_CORE_PRESERVATION_FAILED")
    else:
        surface_region_id[entity == wall_entity_id] = 0
    mapping_rows: list[dict[str, Any]] = []
    boundary_paths: list[str] = []
    boundary_directory.mkdir(parents=True, exist_ok=True)
    for cap_id, boundary in cap_to_boundary.items():
        mask = entity == cap_id
        boundary_type[mask] = 1 if boundary.role == "ASSUMED_INLET" else 2
        boundary_index[mask] = boundary.index
        origin_code[mask] = 1 if boundary.boundary_origin == "CUT_PORT" else 2
        origin[mask] = boundary.boundary_origin
        port_id[mask] = boundary.port_id
        cap_mesh = trimesh.Trimesh(
            vertices=np.asarray(data.points, dtype=float), faces=faces[mask], process=False
        )
        name = "inlet.stl" if boundary.role == "ASSUMED_INLET" else f"outlet_{boundary.index:02d}.stl"
        path = boundary_directory / name
        cap_mesh.export(path)
        boundary_paths.append(str(path.resolve()))
        cap_area = float(cap_mesh.area)
        centroid = centers[mask].mean(axis=0)
        mapping_rows.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "role": boundary.role,
                "boundary_origin": boundary.boundary_origin,
                "vmtk_cap_entity_id": cap_id,
                "triangle_count": int(np.count_nonzero(mask)),
                "area_um2": cap_area,
                "centroid_um": centroid.tolist(),
                "predicted_end_distance_um": float(np.linalg.norm(centroid - boundary.extension_end_um)),
                "stl_path": str(path.resolve()),
            }
        )
    data.cell_data["boundary_type_code"] = boundary_type
    data.cell_data["boundary_index"] = boundary_index
    data.cell_data["boundary_origin_code"] = origin_code
    data.cell_data["boundary_origin"] = origin
    data.cell_data["port_id"] = port_id
    data.cell_data["SurfaceRegionId"] = surface_region_id
    data.cell_data["SurfaceRegion"] = np.asarray(
        ["CORE", "EXTENSION", "CAP"], dtype="<U9"
    )[surface_region_id]
    tagged_vtp = geometry_directory / f"{output_stem}_um.vtp"
    tagged_stl = geometry_directory / f"{output_stem}_um.stl"
    meter_stl = geometry_directory / f"{output_stem}_m.stl"
    data.save(tagged_vtp, binary=True)
    final_mesh = trimesh.Trimesh(
        vertices=np.asarray(data.points, dtype=float), faces=faces, process=False
    )
    final_mesh.export(tagged_stl)
    meter_mesh = final_mesh.copy()
    meter_mesh.vertices = np.asarray(meter_mesh.vertices) * 1.0e-6
    meter_mesh.export(meter_stl)
    manifest = boundary_directory / "boundary_manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "boundary_index",
                "port_id",
                "role",
                "boundary_origin",
                "vmtk_cap_entity_id",
                "triangle_count",
                "area_um2",
                "predicted_end_distance_um",
                "stl_path",
            ),
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(sorted(mapping_rows, key=lambda row: int(row["boundary_index"])))
    outputs = {
        "manual_review_stl": str(tagged_stl.resolve()),
        "tagged_vtp": str(tagged_vtp.resolve()),
        "meter_stl": str(meter_stl.resolve()),
        "boundary_stl_paths": boundary_paths,
        "boundary_manifest_csv": str(manifest.resolve()),
    }
    one_inlet = sum(row["role"] == "ASSUMED_INLET" for row in mapping_rows) == 1
    expected_outlets = (
        sum(row["role"] == "ASSUMED_OUTLET" for row in mapping_rows)
        == len(boundary_list) - 1
    )
    report = {
        "status": (
            "PASS"
            if len(cap_ids) == len(boundary_list)
            and (raw_vtp is None or missing_raw_triangle_count == 0)
            and one_inlet
            and expected_outlets
            else "FAIL"
        ),
        "mapping_method": "one-to-one nearest predicted extension end",
        "wall_entity_id": wall_entity_id,
        "distal_boundary_count": len(cap_ids),
        "one_inlet": one_inlet,
        "three_outlets": sum(
            row["role"] == "ASSUMED_OUTLET" for row in mapping_rows
        )
        == 3,
        "expected_outlet_count": expected_outlets,
        "raw_noncap_triangle_count": preserved_raw_triangle_count,
        "raw_noncap_triangle_missing_count": missing_raw_triangle_count,
        "surface_region_codes": {"CORE": 0, "EXTENSION": 1, "CAP": 2},
        "boundaries": sorted(mapping_rows, key=lambda row: int(row["boundary_index"])),
    }
    return outputs, report


def _outside_local_zones(
    points: np.ndarray,
    boundaries: Iterable[BoundaryInput],
    local: LocalCutConfig,
    *,
    exclude_full_extensions: bool = False,
) -> np.ndarray:
    keep = np.ones(len(points), dtype=bool)
    for boundary in boundaries:
        relative = np.asarray(points, dtype=float) - boundary.center_um
        axial = relative @ boundary.outward_normal
        radial = np.linalg.norm(
            relative - np.outer(axial, boundary.outward_normal), axis=1
        )
        upper = (
            1.1 * boundary.extension_length_um
            if exclude_full_extensions
            else local.local_axial_forward_radius_factor * boundary.source_radius_um
        )
        in_zone = (
            (radial <= local.local_radial_radius_factor * boundary.source_radius_um)
            & (
                axial
                >= -local.local_axial_back_radius_factor
                * boundary.source_radius_um
            )
            & (axial <= upper)
        )
        keep &= ~in_zone
    return keep


def _distance_summary(values: np.ndarray) -> dict[str, float]:
    distances = np.asarray(values, dtype=float)
    if len(distances) == 0:
        return {key: float("inf") for key in ("P50_um", "P95_um", "P99_um", "max_um")}
    return {
        "P50_um": float(np.percentile(distances, 50)),
        "P95_um": float(np.percentile(distances, 95)),
        "P99_um": float(np.percentile(distances, 99)),
        "max_um": float(np.max(distances)),
    }


def _angle_summary(values: np.ndarray) -> dict[str, float]:
    angles = np.asarray(values, dtype=float)
    if len(angles) == 0:
        return {key: float("inf") for key in ("P50_deg", "P95_deg", "P99_deg", "max_deg")}
    return {
        "P50_deg": float(np.percentile(angles, 50)),
        "P95_deg": float(np.percentile(angles, 95)),
        "P99_deg": float(np.percentile(angles, 99)),
        "max_deg": float(np.max(angles)),
    }


def _triangle_key(
    points: np.ndarray, face: np.ndarray
) -> tuple[tuple[float, ...], ...]:
    return tuple(
        sorted(tuple(float(value) for value in point) for point in points[face])
    )


def _face_lookup(
    points: np.ndarray, faces: np.ndarray, face_ids: np.ndarray
) -> dict[tuple[tuple[float, ...], ...], list[int]]:
    lookup: dict[tuple[tuple[float, ...], ...], list[int]] = defaultdict(list)
    for face_id in face_ids:
        lookup[_triangle_key(points, faces[int(face_id)])].append(int(face_id))
    return lookup


def core_exact_preservation_qc(
    original: trimesh.Trimesh,
    final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    local: LocalCutConfig,
) -> dict[str, Any]:
    """Verify original vertices and triangles outside surgery zones exactly."""

    data, final_mesh = polydata_mesh(final_vtp)
    region = np.asarray(data.cell_data["SurfaceRegionId"], dtype=np.uint8)
    final_faces = np.asarray(final_mesh.faces, dtype=np.int64)
    core_face_ids = np.flatnonzero(region == 0)
    final_points = np.asarray(data.points)
    original_points = np.asarray(original.vertices, dtype=float)
    original_faces = np.asarray(original.faces, dtype=np.int64)
    keep_vertices = _outside_local_zones(
        original_points, boundaries, local, exclude_full_extensions=False
    )
    keep_faces = _outside_local_zones(
        original_points[original_faces].mean(axis=1),
        boundaries,
        local,
        exclude_full_extensions=False,
    )
    cast_original = original_points.astype(final_points.dtype, copy=False)
    lookup = _face_lookup(final_points, final_faces, core_face_ids)
    missing = 0
    for face_id in np.flatnonzero(keep_faces):
        values = lookup.get(_triangle_key(cast_original, original_faces[face_id]))
        if values:
            values.pop()
        else:
            missing += 1
    used_core_vertices = np.unique(final_faces[core_face_ids])
    tree = cKDTree(final_points[used_core_vertices].astype(float))
    motion, _ = tree.query(original_points[keep_vertices], k=1)
    cast_motion, _ = tree.query(cast_original[keep_vertices].astype(float), k=1)
    tolerance = _output_precision_tolerance(final_points)
    maximum = float(np.max(motion)) if len(motion) else float("inf")
    p95 = float(np.percentile(motion, 95)) if len(motion) else float("inf")
    cast_maximum = float(np.max(cast_motion)) if len(cast_motion) else float("inf")
    checks = {
        "retained_original_vertices_available": bool(np.any(keep_vertices)),
        "retained_original_vertex_motion_within_output_machine_precision": (
            maximum <= tolerance
        ),
        "retained_vertices_exact_after_output_dtype_cast": cast_maximum == 0.0,
        "core_triangle_connectivity_unchanged": missing == 0,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "retained_original_vertex_count": int(np.count_nonzero(keep_vertices)),
        "retained_original_vertex_max_motion_um": maximum,
        "retained_original_vertex_P95_motion_um": p95,
        "retained_original_vertex_max_motion_after_output_dtype_cast_um": cast_maximum,
        "machine_precision_tolerance_um": tolerance,
        "core_triangle_count_checked": int(np.count_nonzero(keep_faces)),
        "core_triangle_connectivity_changed_count": int(missing),
        "final_core_triangle_count": int(len(core_face_ids)),
        "final_point_dtype": str(final_points.dtype),
    }


def core_symmetric_distance_qc(
    original: trimesh.Trimesh,
    final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    local: LocalCutConfig,
) -> dict[str, Any]:
    """Measure both directions so final-only spikes cannot be hidden."""

    data, final_mesh = polydata_mesh(final_vtp)
    region = np.asarray(data.cell_data["SurfaceRegionId"], dtype=np.uint8)
    original_points = np.asarray(original.vertices, dtype=float)
    original_faces = np.asarray(original.faces, dtype=np.int64)
    original_face_keep = _outside_local_zones(
        original_points[original_faces].mean(axis=1), boundaries, local
    )
    original_core = trimesh.Trimesh(
        vertices=original_points,
        faces=original_faces[original_face_keep],
        process=False,
    )
    final_faces = np.asarray(final_mesh.faces, dtype=np.int64)
    final_centers = np.asarray(final_mesh.triangles_center, dtype=float)
    final_face_keep = (region == 0) & _outside_local_zones(
        final_centers, boundaries, local
    )
    final_core = trimesh.Trimesh(
        vertices=np.asarray(final_mesh.vertices, dtype=float),
        faces=final_faces[final_face_keep],
        process=False,
    )
    original_vertex_ids = np.unique(original_core.faces)
    final_vertex_ids = np.unique(final_core.faces)
    original_samples = np.vstack(
        (original_core.vertices[original_vertex_ids], original_core.triangles_center)
    )
    final_samples = np.vstack(
        (final_core.vertices[final_vertex_ids], final_core.triangles_center)
    )
    _, forward, _ = trimesh.proximity.closest_point(final_core, original_samples)
    _, reverse, _ = trimesh.proximity.closest_point(original_core, final_samples)
    forward_report = _distance_summary(forward)
    reverse_report = _distance_summary(reverse)
    tolerance = _output_precision_tolerance(np.asarray(data.points))
    symmetric_max = max(forward_report["max_um"], reverse_report["max_um"])
    return {
        "status": "PASS" if symmetric_max <= tolerance else "FAIL",
        "method": "bidirectional vertex-plus-face-center closest surface distance",
        "original_core_to_caponly_final_core": forward_report,
        "caponly_final_core_to_original_core": reverse_report,
        "symmetric_max_um": symmetric_max,
        "machine_precision_tolerance_um": tolerance,
        "original_sample_count": int(len(original_samples)),
        "final_sample_count": int(len(final_samples)),
    }


def _inconsistent_adjacent_orientation_count(faces: np.ndarray) -> int:
    edge_faces: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for face_id, face in enumerate(np.asarray(faces, dtype=np.int64)):
        for first, second in zip(face, np.roll(face, -1)):
            key = tuple(sorted((int(first), int(second))))
            direction = 1 if (int(first), int(second)) == key else -1
            edge_faces[key].append((face_id, direction))
    return sum(
        1
        for linked in edge_faces.values()
        if len(linked) == 2 and linked[0][1] == linked[1][1]
    )


def normal_consistency_qc(
    original: trimesh.Trimesh,
    final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    local: LocalCutConfig,
) -> dict[str, Any]:
    data, final_mesh = polydata_mesh(final_vtp)
    faces = np.asarray(final_mesh.faces, dtype=np.int64)
    adjacency = np.asarray(final_mesh.face_adjacency, dtype=np.int64)
    adjacent_dots = np.sum(
        final_mesh.face_normals[adjacency[:, 0]]
        * final_mesh.face_normals[adjacency[:, 1]],
        axis=1,
    )
    adjacent_angles = np.degrees(
        np.arccos(np.clip(adjacent_dots, -1.0, 1.0))
    )
    region = np.asarray(data.cell_data["SurfaceRegionId"], dtype=np.uint8)
    core_ids = np.flatnonzero(region == 0)
    original_points = np.asarray(original.vertices, dtype=float)
    original_faces = np.asarray(original.faces, dtype=np.int64)
    keep = _outside_local_zones(
        original_points[original_faces].mean(axis=1), boundaries, local
    )
    final_points = np.asarray(data.points)
    cast_original = original_points.astype(final_points.dtype, copy=False)
    lookup = _face_lookup(final_points, faces, core_ids)
    cast_mesh = trimesh.Trimesh(
        vertices=cast_original.astype(float), faces=original_faces, process=False
    )
    source_angles: list[float] = []
    cast_angles: list[float] = []
    cast_normal_differences: list[float] = []
    missing = 0
    for source_id in np.flatnonzero(keep):
        values = lookup.get(_triangle_key(cast_original, original_faces[source_id]))
        if not values:
            missing += 1
            continue
        final_id = values.pop()
        for output, source_normal in (
            (source_angles, original.face_normals[source_id]),
            (cast_angles, cast_mesh.face_normals[source_id]),
        ):
            dot = float(
                np.clip(np.dot(source_normal, final_mesh.face_normals[final_id]), -1.0, 1.0)
            )
            output.append(float(np.degrees(np.arccos(dot))))
        cast_normal_differences.append(
            float(
                np.max(
                    np.abs(
                        cast_mesh.face_normals[source_id]
                        - final_mesh.face_normals[final_id]
                    )
                )
            )
        )
    source_array = np.asarray(source_angles, dtype=float)
    cast_array = np.asarray(cast_angles, dtype=float)
    source_report = _angle_summary(source_array)
    cast_report = _angle_summary(cast_array)
    cast_normal_max_difference = (
        float(np.max(cast_normal_differences))
        if cast_normal_differences
        else float("inf")
    )
    cast_exact = cast_normal_max_difference <= 1.0e-12
    inconsistent = _inconsistent_adjacent_orientation_count(faces)
    checks = {
        "winding_consistent": bool(final_mesh.is_winding_consistent),
        "zero_flipped_or_opposing_adjacent_faces": inconsistent == 0,
        "core_face_correspondence_complete": missing == 0,
        "core_normals_exact_after_output_dtype_cast": cast_exact,
    }
    adjacent_report = _angle_summary(adjacent_angles)
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "winding_consistent": bool(final_mesh.is_winding_consistent),
        "flipped_or_opposing_adjacent_face_count": int(inconsistent),
        "adjacent_normal_jump_P50_deg": adjacent_report["P50_deg"],
        "adjacent_normal_jump_P95_deg": adjacent_report["P95_deg"],
        "adjacent_normal_jump_P99_deg": adjacent_report["P99_deg"],
        "adjacent_normal_jump_max_deg": adjacent_report["max_deg"],
        "core_normal_deviation_P50_deg": source_report["P50_deg"],
        "core_normal_deviation_P95_deg": source_report["P95_deg"],
        "core_normal_deviation_P99_deg": source_report["P99_deg"],
        "core_normal_deviation_max_deg": source_report["max_deg"],
        "core_normal_after_output_dtype_cast_P50_deg": cast_report["P50_deg"],
        "core_normal_after_output_dtype_cast_P95_deg": cast_report["P95_deg"],
        "core_normal_after_output_dtype_cast_P99_deg": cast_report["P99_deg"],
        "core_normal_after_output_dtype_cast_max_deg": cast_report["max_deg"],
        "core_normal_after_output_dtype_cast_max_vector_difference": (
            cast_normal_max_difference
        ),
        "core_face_correspondence_missing_count": int(missing),
    }


def previous_global_remesh_diagnostics(
    original: trimesh.Trimesh,
    previous_final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    local: LocalCutConfig,
    *,
    hotspot_count: int = 8,
) -> dict[str, Any]:
    """Diagnose prior global-remesh changes outside all boundary surgery zones."""

    data, previous = polydata_mesh(previous_final_vtp)
    faces = np.asarray(previous.faces, dtype=np.int64)
    centers = np.asarray(previous.triangles_center, dtype=float)
    wall = (
        np.asarray(data.cell_data["boundary_type_code"], dtype=np.uint8) == 0
        if "boundary_type_code" in data.cell_data
        else np.ones(len(faces), dtype=bool)
    )
    previous_keep = wall & _outside_local_zones(
        centers,
        boundaries,
        local,
        exclude_full_extensions=True,
    )
    previous_ids = np.flatnonzero(previous_keep)
    previous_core = trimesh.Trimesh(
        vertices=np.asarray(previous.vertices, dtype=float),
        faces=faces[previous_keep],
        process=False,
    )
    original_points = np.asarray(original.vertices, dtype=float)
    original_faces = np.asarray(original.faces, dtype=np.int64)
    original_keep = _outside_local_zones(
        original_points[original_faces].mean(axis=1), boundaries, local
    )
    original_core = trimesh.Trimesh(
        vertices=original_points, faces=original_faces[original_keep], process=False
    )
    original_vertex_ids = np.unique(original_core.faces)
    previous_vertex_ids = np.unique(previous_core.faces)
    original_samples = np.vstack(
        (original_core.vertices[original_vertex_ids], original_core.triangles_center)
    )
    previous_samples = np.vstack(
        (previous_core.vertices[previous_vertex_ids], previous_core.triangles_center)
    )
    _, forward, _ = trimesh.proximity.closest_point(previous_core, original_samples)
    _, reverse, _ = trimesh.proximity.closest_point(original_core, previous_samples)
    forward_report = _distance_summary(forward)
    reverse_report = _distance_summary(reverse)

    selected_centers = centers[previous_keep]
    _, center_distances, nearest_original = trimesh.proximity.closest_point(
        original, selected_centers
    )
    normal_dots = np.sum(
        previous.face_normals[previous_ids]
        * original.face_normals[np.asarray(nearest_original, dtype=np.int64)],
        axis=1,
    )
    normal_angles = np.degrees(np.arccos(np.clip(normal_dots, -1.0, 1.0)))
    normal_report = _angle_summary(normal_angles)
    distance_scale = max(float(np.percentile(center_distances, 99)), 1.0e-15)
    normal_scale = max(float(np.percentile(normal_angles, 99)), 1.0e-15)
    score = center_distances / distance_scale + normal_angles / normal_scale
    order = np.argsort(score)[::-1]
    separation = max(float(np.median(original.edges_unique_length)) * 5.0, 0.5)
    chosen: list[int] = []
    for local_id in order:
        center = selected_centers[local_id]
        if all(np.linalg.norm(center - selected_centers[other]) >= separation for other in chosen):
            chosen.append(int(local_id))
        if len(chosen) == hotspot_count:
            break
    boundary_list = list(boundaries)
    hotspots = []
    for hotspot_id, local_id in enumerate(chosen):
        center = selected_centers[local_id]
        hotspots.append(
            {
                "hotspot_id": hotspot_id,
                "center_x_um": float(center[0]),
                "center_y_um": float(center[1]),
                "center_z_um": float(center[2]),
                "distance_to_original_um": float(center_distances[local_id]),
                "local_normal_deviation_deg": float(normal_angles[local_id]),
                "nearest_original_face_id": int(nearest_original[local_id]),
                "previous_final_face_id": int(previous_ids[local_id]),
                "distance_from_nearest_CFD_boundary_um": float(
                    min(np.linalg.norm(center - boundary.center_um) for boundary in boundary_list)
                ),
            }
        )
    tolerance = _output_precision_tolerance(np.asarray(data.points))
    symmetric_max = max(forward_report["max_um"], reverse_report["max_um"])
    return {
        "status": "PASS",
        "method": "bidirectional core surface distance outside boundary surgery cylinders",
        "original_core_to_previous_final": forward_report,
        "previous_final_core_to_original": reverse_report,
        "symmetric_max_um": symmetric_max,
        "core_normal_deviation": normal_report,
        "winding_consistent": bool(previous.is_winding_consistent),
        "global_remesh_artifact_detected": symmetric_max > 10.0 * tolerance,
        "machine_precision_reference_um": tolerance,
        "hotspot_count": len(hotspots),
        "hotspots": hotspots,
    }


def meter_scale_qc(um_stl: Path, meter_stl: Path) -> dict[str, Any]:
    um = trimesh.load_mesh(um_stl, process=False)
    meter = trimesh.load_mesh(meter_stl, process=False)
    um_extent = np.asarray(um.extents, dtype=float)
    meter_extent = np.asarray(meter.extents, dtype=float)
    passed = bool(np.allclose(meter_extent, um_extent * 1.0e-6, rtol=1.0e-7, atol=1.0e-14))
    return {
        "status": "PASS" if passed else "FAIL",
        "scale_factor": 1.0e-6,
        "um_extent": um_extent.tolist(),
        "meter_extent": meter_extent.tolist(),
    }


def geometry_pressure_correction(
    final_mesh: trimesh.Trimesh,
    boundaries: Iterable[BoundaryInput],
    proximal_loops: dict[int, BoundaryLoop],
    distal_loops: dict[int, BoundaryLoop],
    *,
    dynamic_viscosity_pa_s: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for boundary in boundaries:
        start = proximal_loops[boundary.index].center_um
        end = distal_loops[boundary.index].center_um
        vector = end - start
        length_um = float(np.linalg.norm(vector))
        direction = vector / length_um
        fractions = np.linspace(0.025, 0.975, 20)
        areas: list[float] = []
        for fraction in fractions:
            section = _section_polygon(final_mesh, start + fraction * vector, direction)
            if section is None:
                raise SurfacePrepareError(
                    f"VMTK_EXTENSION_GEOMETRY_FAILED:{boundary.port_id}:cross_section"
                )
            areas.append(float(section[0]))
        area_array = np.asarray(areas)
        radii_m = np.sqrt(area_array / np.pi) * 1.0e-6
        positions_m = fractions * length_um * 1.0e-6
        integrand = 8.0 * dynamic_viscosity_pa_s / (np.pi * radii_m**4)
        resistance = float(np.trapz(integrand, positions_m))
        pressure_drop = float(boundary.expected_flow_m3_s * resistance)
        is_outlet = boundary.role == "ASSUMED_OUTLET"
        rows.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "role": boundary.role,
                "P_original_1D_pa": boundary.pressure_original_pa,
                "Q_expected_1D_m3_s": boundary.expected_flow_m3_s,
                "station_count": 20,
                "station_fractions": fractions.tolist(),
                "cross_section_area_um2": areas,
                "equivalent_radius_um": (np.sqrt(area_array / np.pi)).tolist(),
                "actual_extension_length_um": length_um,
                "extension_resistance_pa_s_m3": resistance,
                "predicted_extension_pressure_drop_pa": pressure_drop if is_outlet else None,
                "P_solver_boundary_pa": boundary.pressure_original_pa - pressure_drop if is_outlet else None,
                "Q_solver_m3_s": boundary.expected_flow_m3_s if not is_outlet else None,
                "profile": "PARABOLIC" if not is_outlet else None,
                "pressure_correction_role": "NUMERICAL_ARTIFICIAL_EXTENSION_CORRECTION",
                "geometry_source": "final VMTK direct-cap surface; 20 normal cross sections",
            }
        )
    return rows, {
        "method": "integral 8*mu/(pi*r_eq(s)^4) ds",
        "station_count_per_extension": 20,
        "dynamic_viscosity_pa_s": dynamic_viscosity_pa_s,
        "rows": rows,
    }


def write_json(path: Path, value: Any) -> Path:
    path.write_text(json.dumps(value, indent=2, default=str), encoding="utf-8")
    return path
