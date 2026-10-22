"""Output conversion, source validation, and QC for Ultraliser surfaces."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pyvista as pv
import trimesh
import vtk
from shapely.geometry import Point, Polygon

from utils.sampling.sampling_types import ROIRecord
from utils.sampling.structural_features import _branch_paths

from .config import CFDLumenConfig
from .export import write_json
from .types import BranchGeometry, GeometryValidationError, RadiusFidelitySample


def load_mesh(path: Path) -> trimesh.Trimesh:
    loaded = trimesh.load_mesh(path, process=True, validate=False)
    if isinstance(loaded, trimesh.Scene):
        loaded = loaded.to_mesh()
    if not isinstance(loaded, trimesh.Trimesh) or len(loaded.faces) == 0:
        raise ValueError(f"No triangular surface was loaded from {path}")
    return loaded


def _polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces, dtype=np.int64))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def write_geometry_outputs(raw_surface: Path, geometry_directory: Path) -> trimesh.Trimesh:
    """Copy official geometry without smoothing, then create unit-explicit mirrors."""

    geometry_directory.mkdir(parents=True, exist_ok=True)
    surface_um_stl = geometry_directory / "lumen_surface_um.stl"
    surface_um_vtp = geometry_directory / "lumen_surface_um.vtp"
    surface_m_stl = geometry_directory / "lumen_surface_m.stl"
    shutil.copy2(raw_surface, surface_um_stl)
    mesh = load_mesh(surface_um_stl)
    _polydata(mesh).save(surface_um_vtp, binary=True)
    mesh_m = trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float) * 1.0e-6,
        faces=np.asarray(mesh.faces, dtype=np.int64),
        process=False,
    )
    mesh_m.export(surface_m_stl)
    return mesh


def validate_source_roi(roi: ROIRecord) -> tuple[list[BranchGeometry], dict[str, Any]]:
    """Validate and decompose a saved ROI without resampling or modifying it."""

    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    edges = np.asarray(roi.local_edges, dtype=np.int64).reshape((-1, 2))
    node_ids = np.asarray(roi.local_node_ids, dtype=np.int64)
    failures: list[dict[str, Any]] = []
    if positions.shape != (len(node_ids), 3):
        failures.append({"reason": "node_position_shape", "value": list(positions.shape)})
    if radii.shape != (len(node_ids),):
        failures.append({"reason": "node_radius_shape", "value": list(radii.shape)})
    if not np.array_equal(node_ids, np.arange(len(node_ids), dtype=np.int64)):
        failures.append({"reason": "local_node_ids_are_not_contiguous_indices"})
    invalid_coordinates = np.flatnonzero(~np.all(np.isfinite(positions), axis=1))
    invalid_radii = np.flatnonzero(~np.isfinite(radii) | (radii <= 0.0))
    failures.extend(
        {"reason": "nonfinite_coordinate", "node_id": int(index)}
        for index in invalid_coordinates
    )
    failures.extend(
        {
            "reason": "invalid_radius",
            "node_id": int(index),
            "radius_um": float(radii[index]),
        }
        for index in invalid_radii
    )
    if len(edges) == 0:
        failures.append({"reason": "roi_has_no_edges"})
        lengths = np.empty(0, dtype=float)
    elif int(edges.min()) < 0 or int(edges.max()) >= len(node_ids):
        failures.append({"reason": "edge_node_index_out_of_range"})
        lengths = np.empty(0, dtype=float)
    else:
        lengths = np.linalg.norm(positions[edges[:, 1]] - positions[edges[:, 0]], axis=1)
        for edge_index in np.flatnonzero(~np.isfinite(lengths) | (lengths <= 1.0e-12)):
            failures.append(
                {
                    "reason": "zero_or_near_zero_edge",
                    "edge_index": int(edge_index),
                    "length_um": float(lengths[edge_index]),
                }
            )

    graph = nx.Graph()
    graph.add_nodes_from(map(int, node_ids))
    graph.add_edges_from((int(first), int(second)) for first, second in edges)
    component_count = nx.number_connected_components(graph) if graph else 0
    duplicate_edges = len(edges) - len({tuple(sorted(map(int, edge))) for edge in edges})
    if duplicate_edges:
        failures.append({"reason": "duplicate_local_edges", "count": duplicate_edges})
    if component_count != 1:
        failures.append({"reason": "roi_not_connected", "component_count": component_count})
    if failures:
        error = GeometryValidationError(
            f"ROI {roi.roi_id} failed source validation: {failures[0]['reason']}"
        )
        error.failures = failures  # type: ignore[attr-defined]
        raise error

    edge_index_by_nodes = {
        tuple(sorted((int(first), int(second)))): edge_index
        for edge_index, (first, second) in enumerate(edges)
    }
    branches: list[BranchGeometry] = []
    for branch_id, path in enumerate(_branch_paths(graph)):
        source_edge_ids = tuple(
            int(roi.local_edge_global_ids[edge_index_by_nodes[tuple(sorted((first, second)))]] )
            for first, second in zip(path[:-1], path[1:])
        )
        raw_points = positions[path].copy()
        raw_radii = radii[path].copy()
        arc_length = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(raw_points, axis=0), axis=1)))
        )
        branches.append(
            BranchGeometry(
                branch_id=branch_id,
                local_node_ids=tuple(map(int, path)),
                source_global_nodes=tuple(int(roi.local_node_global_ids[node]) for node in path),
                source_global_edges=source_edge_ids,
                raw_points_um=raw_points,
                raw_radius_um=raw_radii,
                points_um=raw_points.copy(),
                radius_um=raw_radii.copy(),
                arc_length_um=arc_length,
            )
        )
    covered_edges = [edge for branch in branches for edge in branch.source_global_edges]
    expected_edges = list(map(int, roi.local_edge_global_ids))
    if sorted(covered_edges) != sorted(expected_edges) or len(covered_edges) != len(expected_edges):
        raise GeometryValidationError(
            f"ROI {roi.roi_id} branch extraction did not preserve every source edge exactly once"
        )
    cycle_rank = graph.number_of_edges() - graph.number_of_nodes() + component_count
    report = {
        "roi_id": roi.roi_id,
        "node_count": int(len(node_ids)),
        "edge_count": int(len(edges)),
        "branch_count": len(branches),
        "connected_component_count": int(component_count),
        "cycle_rank": int(cycle_rank),
        "cut_port_count": len(roi.cut_ports),
        "radius_min_um": float(radii.min()),
        "radius_median_um": float(np.median(radii)),
        "radius_max_um": float(radii.max()),
        "centerline_total_length_um": float(lengths.sum()),
        "source_geometry_modified": False,
        "source_radius_modified": False,
        "source_edge_coverage_exactly_once": True,
    }
    return branches, report


def _triangle_intersections(
    mesh: trimesh.Trimesh,
    face_ids: np.ndarray,
) -> tuple[list[tuple[int, int]], int]:
    """Find non-adjacent triangle intersections using an R-tree AABB index."""

    if len(face_ids) < 2:
        return [], 0
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces[face_ids], dtype=np.int64)
    triangles = vertices[faces]
    triangle_minimum = triangles.min(axis=1)
    triangle_maximum = triangles.max(axis=1)
    bounds = np.column_stack((triangle_minimum, triangle_maximum))
    tree = trimesh.util.bounds_tree(bounds)
    pairs: list[tuple[int, int]] = []
    candidates_checked = 0
    for local_id, current in enumerate(triangles):
        face_vertices = set(map(int, faces[local_id]))
        for other_id in tree.intersection(bounds[local_id]):
            other_id = int(other_id)
            if other_id <= local_id or face_vertices.intersection(map(int, faces[other_id])):
                continue
            if np.any(triangle_maximum[local_id] < triangle_minimum[other_id]) or np.any(
                triangle_maximum[other_id] < triangle_minimum[local_id]
            ):
                continue
            candidates_checked += 1
            if vtk.vtkTriangle.TrianglesIntersect(*current, *triangles[other_id]):
                pairs.append((int(face_ids[local_id]), int(face_ids[other_id])))
    return pairs, candidates_checked


def evaluate_surface_topology(
    mesh: trimesh.Trimesh,
    config: CFDLumenConfig,
) -> dict[str, Any]:
    sorted_edges = np.sort(np.asarray(mesh.edges, dtype=np.int64), axis=1)
    _, edge_counts = np.unique(sorted_edges, axis=0, return_counts=True)
    boundary_edges = int(np.count_nonzero(edge_counts == 1))
    nonmanifold_edges = int(np.count_nonzero(edge_counts > 2))
    areas = np.asarray(mesh.area_faces, dtype=float)
    diagonal = float(np.linalg.norm(np.ptp(np.asarray(mesh.vertices), axis=0)))
    area_tolerance = max(np.finfo(float).eps * diagonal**2 * 100.0, 1.0e-18)
    repeated_index = np.asarray(
        [len(set(map(int, face))) < 3 for face in np.asarray(mesh.faces)], dtype=bool
    )
    degenerate = int(np.count_nonzero((areas <= area_tolerance) | repeated_index))
    intersections, candidate_pairs = _triangle_intersections(
        mesh, np.arange(len(mesh.faces), dtype=np.int64)
    )
    components = mesh.split(only_watertight=False)
    qc = config.surface_qc
    checks = {
        "watertight": not qc.require_watertight or bool(mesh.is_watertight),
        "single_component": not qc.require_single_component or len(components) == 1,
        "zero_boundary_edges": not qc.require_zero_boundary_edges or boundary_edges == 0,
        "zero_nonmanifold_edges": not qc.require_zero_nonmanifold_edges
        or nonmanifold_edges == 0,
        "zero_self_intersections": not qc.require_zero_self_intersections
        or len(intersections) == 0,
        "zero_degenerate_triangles": not qc.require_zero_degenerate_triangles
        or degenerate == 0,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "vertex_count": int(len(mesh.vertices)),
        "triangle_count": int(len(mesh.faces)),
        "component_count": int(len(components)),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "boundary_edge_count": boundary_edges,
        "nonmanifold_edge_count": nonmanifold_edges,
        "self_intersection_count": int(len(intersections)),
        "self_intersection_candidate_pairs_checked": int(candidate_pairs),
        "degenerate_triangle_count": degenerate,
        "degenerate_area_tolerance_um2": area_tolerance,
        "surface_area_um2": float(mesh.area),
        "volume_um3": float(abs(mesh.volume)),
        "bounds_um": np.asarray(mesh.bounds, dtype=float).tolist(),
    }


def _orthogonal_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    normal = normal / np.linalg.norm(normal)
    helper = (
        np.asarray((1.0, 0.0, 0.0))
        if abs(normal[0]) < 0.8
        else np.asarray((0.0, 1.0, 0.0))
    )
    first = np.cross(normal, helper)
    first /= np.linalg.norm(first)
    return first, np.cross(normal, first)


def _section_polygon(
    mesh: trimesh.Trimesh,
    center: np.ndarray,
    tangent: np.ndarray,
) -> tuple[float, tuple[tuple[float, float], ...]] | None:
    section = mesh.section(plane_origin=center, plane_normal=tangent)
    if section is None:
        return None
    basis_first, basis_second = _orthogonal_basis(tangent)
    candidates: list[tuple[bool, float, float, np.ndarray]] = []
    for discrete in section.discrete:
        points = np.asarray(discrete, dtype=float)
        if len(points) < 4:
            continue
        relative = points - center
        projected = np.column_stack((relative @ basis_first, relative @ basis_second))
        polygon = Polygon(projected)
        if not polygon.is_valid:
            polygon = polygon.buffer(0)
        if polygon.is_empty or polygon.area <= 0:
            continue
        polygons = polygon.geoms if polygon.geom_type == "MultiPolygon" else (polygon,)
        for component in polygons:
            if component.is_empty or component.area <= 0:
                continue
            origin = Point(0.0, 0.0)
            candidates.append(
                (
                    bool(component.buffer(1.0e-9).contains(origin)),
                    float(component.distance(origin)),
                    float(component.area),
                    np.asarray(component.exterior.coords, dtype=float),
                )
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (not item[0], item[1], -item[2]))
    _, _, area, coordinates = candidates[0]
    return area, tuple((float(point[0]), float(point[1])) for point in coordinates)


def _tangent(points: np.ndarray, index: int) -> np.ndarray:
    if index == 0:
        vector = points[1] - points[0]
    elif index == len(points) - 1:
        vector = points[-1] - points[-2]
    else:
        vector = points[index + 1] - points[index - 1]
    return vector / np.linalg.norm(vector)


def evaluate_radius_fidelity(
    mesh: trimesh.Trimesh,
    branches: list[BranchGeometry],
    roi: ROIRecord,
    config: CFDLumenConfig,
) -> tuple[list[RadiusFidelitySample], dict[str, Any]]:
    degree = np.bincount(np.asarray(roi.local_edges).ravel(), minlength=roi.node_count)
    cut_ids = {int(port.local_node_id) for port in roi.cut_ports}
    samples: list[RadiusFidelitySample] = []
    attempted = 0
    missed = 0
    for branch in branches:
        length = float(branch.arc_length_um[-1])
        first_local = branch.local_node_ids[0]
        last_local = branch.local_node_ids[-1]
        start_exclusion = (
            2.0 * branch.radius_um[0] * config.surface_qc.radius_fidelity_skip_diameters
            if degree[first_local] >= 3 or first_local in cut_ids
            else 0.0
        )
        end_exclusion = (
            2.0 * branch.radius_um[-1] * config.surface_qc.radius_fidelity_skip_diameters
            if degree[last_local] >= 3 or last_local in cut_ids
            else 0.0
        )
        valid = np.flatnonzero(
            (branch.arc_length_um > start_exclusion)
            & ((length - branch.arc_length_um) > end_exclusion)
        )
        if len(valid) == 0:
            continue
        count = min(config.surface_qc.radius_fidelity_samples_per_branch, len(valid))
        chosen = np.unique(valid[np.linspace(0, len(valid) - 1, count).round().astype(int)])
        for index in chosen:
            attempted += 1
            center = branch.points_um[index]
            tangent = _tangent(branch.points_um, int(index))
            section = _section_polygon(mesh, center, tangent)
            if section is None:
                missed += 1
                continue
            area, coordinates = section
            reconstructed = float(np.sqrt(area / np.pi))
            source = float(branch.radius_um[index])
            samples.append(
                RadiusFidelitySample(
                    branch_id=branch.branch_id,
                    sample_index=int(index),
                    arc_length_um=float(branch.arc_length_um[index]),
                    center_um=tuple(map(float, center)),
                    tangent=tuple(map(float, tangent)),
                    source_radius_um=source,
                    reconstructed_radius_um=reconstructed,
                    relative_error=float((reconstructed - source) / source),
                    section_xy_um=coordinates,
                )
            )
    signed = np.asarray([sample.relative_error for sample in samples], dtype=float)
    absolute = np.abs(signed)
    if len(samples) == 0:
        metrics: dict[str, Any] = {
            "median_signed_relative_error": None,
            "mean_signed_relative_error": None,
            "median_absolute_relative_error": None,
            "mean_absolute_relative_error": None,
            "p95_absolute_relative_error": None,
            "max_absolute_relative_error": None,
            "positive_error_fraction": None,
        }
    else:
        metrics = {
            "median_signed_relative_error": float(np.median(signed)),
            "mean_signed_relative_error": float(np.mean(signed)),
            "median_absolute_relative_error": float(np.median(absolute)),
            "mean_absolute_relative_error": float(np.mean(absolute)),
            "p95_absolute_relative_error": float(np.percentile(absolute, 95)),
            "max_absolute_relative_error": float(np.max(absolute)),
            "positive_error_fraction": float(np.mean(signed > 0.0)),
        }
    p95 = metrics["p95_absolute_relative_error"]
    threshold = config.surface_qc.max_radius_p95_error
    return samples, {
        "status": "PASS" if p95 is not None and float(p95) <= threshold else "FAIL",
        "max_radius_p95_error": threshold,
        "attempted_sample_count": attempted,
        "successful_sample_count": len(samples),
        "missed_section_count": missed,
        **metrics,
    }


def finalize_ultraliser_outputs(
    roi: ROIRecord,
    config: CFDLumenConfig,
    *,
    raw_surface: Path,
    geometry_directory: Path,
    qc_directory: Path,
) -> dict[str, Any]:
    mesh = write_geometry_outputs(raw_surface, geometry_directory)
    surface_report = evaluate_surface_topology(mesh, config)
    branches, source_report = validate_source_roi(roi)
    samples, radius_report = evaluate_radius_fidelity(mesh, branches, roi, config)
    radius_report["samples"] = [sample.report() for sample in samples]
    write_json(qc_directory / "surface_qc.json", surface_report)
    write_json(qc_directory / "radius_fidelity.json", radius_report)
    return {
        "source_qc": source_report,
        "surface_qc": surface_report,
        "radius_fidelity": radius_report,
    }


def write_reconstruction_report(path: Path, summary: dict[str, Any]) -> Path:
    surface = summary["surface_qc"]
    radius = summary["radius_fidelity"]
    lines = [
        "# Ultraliser vascular surface reconstruction",
        "",
        f"- ROI: `{summary['roi_id']}`",
        f"- Source SWC modified: `{summary['source_swc_modified']}`",
        f"- Source radius modified: `{summary['source_radius_modified']}`",
        f"- H5 feed radius scale: `{summary['radius_scale']}`",
        f"- Ultraliser invocation count: `{summary['ultraliser_invocation_count']}`",
        f"- Watertight: `{surface['watertight']}`",
        f"- Components: `{surface['component_count']}`",
        f"- Self intersections: `{surface['self_intersection_count']}`",
        f"- Nonmanifold edges: `{surface['nonmanifold_edge_count']}`",
        f"- Degenerate triangles: `{surface['degenerate_triangle_count']}`",
        f"- Median signed radius error: `{radius['median_signed_relative_error']}`",
        f"- P95 absolute radius error: `{radius['p95_absolute_relative_error']}`",
        f"- Status: `{summary['status']}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
