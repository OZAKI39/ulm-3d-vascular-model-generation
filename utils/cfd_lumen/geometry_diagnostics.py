"""Port and junction measurements required by the v2 diagnostic protocol."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh
from shapely.geometry import Polygon

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .lumen_builder import balanced_manifold_union
from .surface_qc import _section_polygon
from .types import BranchGeometry, LumenPrimitives, PortGeometry


def _mesh_polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces, dtype=np.int64))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def save_diagnostic_primitives(
    primitives: LumenPrimitives,
    post_boolean_mesh: trimesh.Trimesh,
    root: Path,
) -> list[Path]:
    branch_dir = root / "branch_tubes"
    junction_dir = root / "junction_solids"
    port_dir = root / "port_extensions"
    for folder in (branch_dir, junction_dir, port_dir):
        folder.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    blocks: list[pv.PolyData] = []
    for branch_id, mesh in primitives.branch_tubes.items():
        path = branch_dir / f"branch_tube_{branch_id:03d}.vtp"
        polydata = _mesh_polydata(mesh)
        polydata.cell_data["primitive_type"] = np.full(polydata.n_cells, 1, dtype=np.int8)
        polydata.cell_data["primitive_id"] = np.full(polydata.n_cells, branch_id, dtype=np.int32)
        polydata.save(path)
        paths.append(path)
        blocks.append(polydata)
    for node_id, mesh in primitives.junction_solids.items():
        path = junction_dir / f"junction_solid_node_{node_id:03d}.vtp"
        polydata = _mesh_polydata(mesh)
        polydata.cell_data["primitive_type"] = np.full(polydata.n_cells, 2, dtype=np.int8)
        polydata.cell_data["primitive_id"] = np.full(polydata.n_cells, node_id, dtype=np.int32)
        polydata.save(path)
        paths.append(path)
        blocks.append(polydata)
    for port_id, mesh in primitives.port_extensions.items():
        path = port_dir / f"port_extension_{port_id:03d}.vtp"
        polydata = _mesh_polydata(mesh)
        polydata.cell_data["primitive_type"] = np.full(polydata.n_cells, 3, dtype=np.int8)
        polydata.cell_data["primitive_id"] = np.full(polydata.n_cells, port_id, dtype=np.int32)
        polydata.save(path)
        paths.append(path)
        blocks.append(polydata)
    combined = (
        pv.MultiBlock(blocks)
        .combine(merge_points=False)
        .extract_surface(algorithm="dataset_surface")
        .triangulate()
    )
    combined_path = root / "pre_boolean_combined.vtp"
    combined.save(combined_path)
    post_path = root / "post_boolean_surface.vtp"
    _mesh_polydata(post_boolean_mesh).save(post_path)
    return [*paths, combined_path, post_path]


def _intersection_volume(
    first: trimesh.Trimesh,
    second: trimesh.Trimesh,
) -> tuple[float | None, str | None, float]:
    started = time.perf_counter()
    try:
        result = trimesh.boolean.intersection(
            [first, second], engine="manifold", check_volume=True
        )
        if not isinstance(result, trimesh.Trimesh) or len(result.faces) == 0:
            return 0.0, None, time.perf_counter() - started
        return float(abs(result.volume)), None, time.perf_counter() - started
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}", time.perf_counter() - started


def _endpoint_branch(
    branches: list[BranchGeometry],
    local_node_id: int,
) -> tuple[BranchGeometry, int]:
    matches: list[tuple[BranchGeometry, int]] = []
    for branch in branches:
        if branch.local_node_ids[0] == local_node_id:
            matches.append((branch, 0))
        if branch.local_node_ids[-1] == local_node_id:
            matches.append((branch, -1))
    if len(matches) != 1:
        raise RuntimeError(
            f"CUT_PORT local node {local_node_id} maps to {len(matches)} branch endpoints"
        )
    return matches[0]


def diagnose_ports(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    primitives: LumenPrimitives,
    config: CFDLumenConfig,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Measure source, resampled, tube, extension, and overlap continuity per port."""

    edge_ids = np.asarray(roi.local_edge_global_ids, dtype=np.int64)
    edge_points = np.asarray(roi.local_edge_points_um, dtype=float)
    edge_radius = np.asarray(roi.local_edge_radius_um, dtype=float)
    rows: list[dict[str, Any]] = []
    for cut, port in zip(roi.cut_ports, ports):
        edge_match = np.flatnonzero(edge_ids == int(cut.global_edge_id))
        if len(edge_match) != 1:
            raise RuntimeError(
                f"CUT_PORT {cut.cut_port_id} source edge maps to {len(edge_match)} ROI edges"
            )
        edge_index = int(edge_match[0])
        edge_start, edge_end = edge_points[edge_index]
        radius_start, radius_end = edge_radius[edge_index]
        vector = edge_end - edge_start
        projection = float(np.dot(np.asarray(cut.intersection_position_um) - edge_start, vector))
        projection /= float(np.dot(vector, vector))
        recomputed_radius = float((1.0 - projection) * radius_start + projection * radius_end)
        branch, endpoint_index = _endpoint_branch(branches, int(cut.local_node_id))
        branch_position = np.asarray(branch.points_um[endpoint_index], dtype=float)
        branch_radius = float(branch.radius_um[endpoint_index])
        position_error = float(
            np.linalg.norm(branch_position - np.asarray(cut.intersection_position_um, dtype=float))
        )
        endpoint_radius_error = float(abs(branch_radius - float(cut.radius_at_cut_um)))
        adjacent_index = 1 if endpoint_index == 0 else len(branch.points_um) - 2
        available = float(np.linalg.norm(branch.points_um[adjacent_index] - branch_position))
        epsilon = max(min(0.25 * port.radius_um, 0.25 * available), 1.0e-6)
        branch_plane = branch_position - epsilon * port.outward_tangent
        extension_plane = branch_position + epsilon * port.outward_tangent
        branch_section = _section_polygon(
            primitives.branch_tubes[branch.branch_id], branch_plane, port.outward_tangent
        )
        extension_section = _section_polygon(
            primitives.port_extensions[port.port_id], extension_plane, port.outward_tangent
        )
        branch_area = float(branch_section[0]) if branch_section else None
        extension_area = float(extension_section[0]) if extension_section else None
        branch_mesh_radius = float(np.sqrt(branch_area / np.pi)) if branch_area else None
        extension_mesh_radius = (
            float(np.sqrt(extension_area / np.pi)) if extension_area else None
        )
        radius_step = (
            abs(branch_mesh_radius - extension_mesh_radius)
            / ((branch_mesh_radius + extension_mesh_radius) * 0.5)
            if branch_mesh_radius is not None and extension_mesh_radius is not None
            else None
        )
        area_step = (
            abs(branch_area - extension_area) / ((branch_area + extension_area) * 0.5)
            if branch_area is not None and extension_area is not None
            else None
        )
        polygon_shape_error = None
        if branch_section and extension_section:
            branch_polygon = Polygon(branch_section[1])
            extension_polygon = Polygon(extension_section[1])
            if branch_polygon.is_valid and extension_polygon.is_valid:
                mean_area = 0.5 * (branch_polygon.area + extension_polygon.area)
                polygon_shape_error = (
                    float(branch_polygon.symmetric_difference(extension_polygon).area / mean_area)
                    if mean_area > 0
                    else None
                )
        intersection_volume, intersection_error, intersection_runtime = _intersection_volume(
            primitives.branch_tubes[branch.branch_id],
            primitives.port_extensions[port.port_id],
        )
        expected_overlap = float(port.overlap_length_um)
        actual_overlap = float(
            np.dot(branch_position - port.cylinder_start_um, port.outward_tangent)
        )
        root_candidates: list[str] = []
        tolerance = max(config.geometry.minimum_edge_length_um * 10.0, 1.0e-8)
        if position_error > tolerance:
            root_candidates.append("PORT_RESAMPLING_ENDPOINT_MISMATCH")
        if endpoint_radius_error > tolerance:
            root_candidates.append("PORT_RADIUS_MISMATCH")
        if intersection_volume is not None and intersection_volume <= 1.0e-12:
            root_candidates.append("PORT_INSUFFICIENT_OVERLAP")
        branch_vertices = len(branch_section[1]) - 1 if branch_section else None
        extension_vertices = len(extension_section[1]) - 1 if extension_section else None
        if branch_vertices is not None and extension_vertices is not None and branch_vertices != extension_vertices:
            root_candidates.append("PORT_POLYGON_RESOLUTION_MISMATCH")
        elif polygon_shape_error is not None and polygon_shape_error > 0.01:
            root_candidates.append("PORT_POLYGON_RESOLUTION_MISMATCH")
        rows.append(
            {
                "port_id": port.port_id,
                "cut_port_id": cut.cut_port_id,
                "exact_cut_x_um": float(cut.intersection_position_um[0]),
                "exact_cut_y_um": float(cut.intersection_position_um[1]),
                "exact_cut_z_um": float(cut.intersection_position_um[2]),
                "cut_radius_um": float(cut.radius_at_cut_um),
                "source_global_edge": int(cut.global_edge_id),
                "source_edge_start_x_um": float(edge_start[0]),
                "source_edge_start_y_um": float(edge_start[1]),
                "source_edge_start_z_um": float(edge_start[2]),
                "source_edge_end_x_um": float(edge_end[0]),
                "source_edge_end_y_um": float(edge_end[1]),
                "source_edge_end_z_um": float(edge_end[2]),
                "source_edge_start_radius_um": float(radius_start),
                "source_edge_end_radius_um": float(radius_end),
                "cut_projection_t": projection,
                "radius_recomputed_um": recomputed_radius,
                "radius_interpolation_error_um": abs(float(cut.radius_at_cut_um) - recomputed_radius),
                "branch_id": branch.branch_id,
                "endpoint_position_error_um": position_error,
                "endpoint_radius_error_um": endpoint_radius_error,
                "endpoint_radius_relative_error": endpoint_radius_error / float(cut.radius_at_cut_um),
                "section_offset_epsilon_um": epsilon,
                "source_radius_um": float(cut.radius_at_cut_um),
                "branch_mesh_area_um2": branch_area,
                "branch_mesh_radius_um": branch_mesh_radius,
                "extension_radius_input_um": port.radius_um,
                "extension_mesh_area_um2": extension_area,
                "extension_mesh_radius_um": extension_mesh_radius,
                "step_radius_relative_error": radius_step,
                "step_area_relative_error": area_step,
                "branch_tube_sides": config.geometry.tube_sides,
                "extension_cylinder_sides": config.geometry.tube_sides,
                "cross_section_vertex_count_branch": branch_vertices,
                "cross_section_vertex_count_extension": extension_vertices,
                "cross_section_polygon_symmetric_difference_relative": polygon_shape_error,
                "expected_overlap_um": expected_overlap,
                "actual_axial_overlap_um": actual_overlap,
                "intersection_volume_um3": intersection_volume,
                "intersection_boolean_error": intersection_error,
                "intersection_boolean_runtime_s": intersection_runtime,
                "root_cause_candidates": ";".join(root_candidates),
            }
        )
    valid_step = [row["step_area_relative_error"] for row in rows if row["step_area_relative_error"] is not None]
    valid_overlap = [row["intersection_volume_um3"] for row in rows if row["intersection_volume_um3"] is not None]
    summary = {
        "port_count": len(rows),
        "max_endpoint_position_error_um": max((row["endpoint_position_error_um"] for row in rows), default=None),
        "max_endpoint_radius_error_um": max((row["endpoint_radius_error_um"] for row in rows), default=None),
        "max_radius_interpolation_error_um": max((row["radius_interpolation_error_um"] for row in rows), default=None),
        "max_step_area_relative_error": max(valid_step, default=None),
        "min_intersection_volume_um3": min(valid_overlap, default=None),
    }
    return rows, summary


def _sample_away_from_endpoint(
    branch: BranchGeometry,
    endpoint_index: int,
    distance: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    if endpoint_index == 0:
        points = branch.points_um
        radii = branch.radius_um
    else:
        points = branch.points_um[::-1]
        radii = branch.radius_um[::-1]
    cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))))
    actual = min(float(distance), float(cumulative[-1]) * 0.9)
    index = min(max(int(np.searchsorted(cumulative, actual)), 1), len(points) - 1)
    span = cumulative[index] - cumulative[index - 1]
    fraction = 0.0 if span <= 0 else (actual - cumulative[index - 1]) / span
    point = (1.0 - fraction) * points[index - 1] + fraction * points[index]
    radius = float((1.0 - fraction) * radii[index - 1] + fraction * radii[index])
    tangent = points[index] - points[index - 1]
    tangent /= np.linalg.norm(tangent)
    return point, tangent, radius, actual


def _surviving_cap_faces(
    mesh: trimesh.Trimesh,
    center: np.ndarray,
    tangent: np.ndarray,
    radius: float,
    confirmed_internal_face_ids: set[int],
) -> tuple[int, int]:
    centroids = np.asarray(mesh.triangles_center)
    relative = centroids - center
    axial = np.abs(relative @ tangent)
    radial_vector = relative - (relative @ tangent)[:, None] * tangent
    radial = np.linalg.norm(radial_vector, axis=1)
    alignment = np.abs(np.asarray(mesh.face_normals) @ tangent)
    tolerance = max(radius * 1.0e-4, 1.0e-8)
    candidates = np.flatnonzero(
        (axial <= tolerance) & (radial <= radius * 1.01) & (alignment >= 0.95)
    )
    confirmed = sum(int(face_id) in confirmed_internal_face_ids for face_id in candidates)
    return int(len(candidates)), int(confirmed)


def diagnose_junctions(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    primitives: LumenPrimitives,
    post_boolean_mesh: trimesh.Trimesh,
    config: CFDLumenConfig,
    *,
    boolean_runtime_s: float,
    confirmed_internal_face_ids: set[int] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    degree = np.bincount(np.asarray(roi.local_edges).ravel(), minlength=roi.node_count)
    rows: list[dict[str, Any]] = []
    internal_ids = confirmed_internal_face_ids or set()
    for node_id in np.flatnonzero(degree >= 3):
        center = np.asarray(roi.local_node_positions_um[node_id], dtype=float)
        source_radius = float(roi.local_node_radius_um[node_id])
        solid = primitives.junction_solids[int(node_id)]
        adjacent: list[tuple[BranchGeometry, int]] = []
        for branch in branches:
            if branch.local_node_ids[0] == int(node_id):
                adjacent.append((branch, 0))
            elif branch.local_node_ids[-1] == int(node_id):
                adjacent.append((branch, -1))
        local_radius = max([source_radius, *[float(branch.radius_um[index]) for branch, index in adjacent]])
        local_face_count = int(
            np.count_nonzero(np.linalg.norm(post_boolean_mesh.triangles_center - center, axis=1) <= 3.0 * local_radius)
        )
        for branch, endpoint_index in adjacent:
            endpoint_radius = float(branch.radius_um[endpoint_index])
            next_index = 1 if endpoint_index == 0 else len(branch.points_um) - 2
            away = np.asarray(branch.points_um[next_index] - branch.points_um[endpoint_index], dtype=float)
            away /= np.linalg.norm(away)
            intersection_volume, intersection_error, intersection_runtime = _intersection_volume(
                solid, primitives.branch_tubes[branch.branch_id]
            )
            sections: dict[str, float | None] = {}
            section_radii: dict[str, float | None] = {}
            actual_distances: dict[str, float] = {}
            diameter = 2.0 * endpoint_radius
            for label, distance in (("1D", diameter), ("2D", 2.0 * diameter)):
                point, tangent, _, actual = _sample_away_from_endpoint(
                    branch, endpoint_index, distance
                )
                section = _section_polygon(post_boolean_mesh, point, tangent)
                area = float(section[0]) if section else None
                sections[label] = area
                section_radii[label] = float(np.sqrt(area / np.pi)) if area else None
                actual_distances[label] = actual
            cap_plane_candidates, cap_survivors = _surviving_cap_faces(
                post_boolean_mesh, center, away, endpoint_radius, internal_ids
            )
            rows.append(
                {
                    "junction_node_id": int(node_id),
                    "junction_x_um": float(center[0]),
                    "junction_y_um": float(center[1]),
                    "junction_z_um": float(center[2]),
                    "branch_id": branch.branch_id,
                    "adjacent_role": "undirected_adjacent_branch",
                    "junction_swc_radius_um": source_radius,
                    "junction_primitive_radius_um": source_radius * config.junction.radius_scale,
                    "adjacent_branch_radius_um": endpoint_radius,
                    "junction_to_adjacent_radius_ratio": source_radius / endpoint_radius,
                    "intersection_volume_um3": intersection_volume,
                    "intersection_over_junction_volume": (
                        intersection_volume / abs(float(solid.volume))
                        if intersection_volume is not None and abs(float(solid.volume)) > 0
                        else None
                    ),
                    "intersection_boolean_error": intersection_error,
                    "intersection_boolean_runtime_s": intersection_runtime,
                    "branch_internal_cap_present_before_boolean": True,
                    "cap_plane_candidate_face_count_after_boolean": cap_plane_candidates,
                    "suspected_internal_cap_face_count_after_boolean": cap_survivors,
                    "pre_boolean_vertices": int(len(solid.vertices) + len(primitives.branch_tubes[branch.branch_id].vertices)),
                    "pre_boolean_triangles": int(len(solid.faces) + len(primitives.branch_tubes[branch.branch_id].faces)),
                    "post_boolean_local_triangles": local_face_count,
                    "number_of_solids_unioned": len(primitives.all_meshes),
                    "boolean_backend": "manifold",
                    "boolean_runtime_s": boolean_runtime_s,
                    "boolean_warnings_or_exceptions": None,
                    "section_1D_actual_distance_um": actual_distances["1D"],
                    "section_1D_area_um2": sections["1D"],
                    "section_1D_equivalent_radius_um": section_radii["1D"],
                    "section_2D_actual_distance_um": actual_distances["2D"],
                    "section_2D_area_um2": sections["2D"],
                    "section_2D_equivalent_radius_um": section_radii["2D"],
                }
            )
    volumes = [row["intersection_volume_um3"] for row in rows if row["intersection_volume_um3"] is not None]
    ratios = [row["junction_to_adjacent_radius_ratio"] for row in rows]
    areas = [row[key] for row in rows for key in ("section_1D_area_um2", "section_2D_area_um2") if row[key] is not None]
    summary = {
        "junction_count": int(np.count_nonzero(degree >= 3)),
        "junction_branch_pair_count": len(rows),
        "min_intersection_volume_um3": min(volumes, default=None),
        "min_junction_to_adjacent_radius_ratio": min(ratios, default=None),
        "max_junction_to_adjacent_radius_ratio": max(ratios, default=None),
        "suspected_internal_cap_face_count": sum(row["suspected_internal_cap_face_count_after_boolean"] for row in rows),
        "A_min_local_um2": min(areas, default=None),
        "A_max_local_um2": max(areas, default=None),
    }
    return rows, summary


def explicit_union_for_diagnostics(primitives: LumenPrimitives) -> tuple[trimesh.Trimesh, float]:
    started = time.perf_counter()
    mesh = balanced_manifold_union(primitives.all_meshes)
    return mesh, time.perf_counter() - started
