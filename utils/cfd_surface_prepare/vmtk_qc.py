"""Integration quality controls for official VMTK-derived CFD surfaces."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import numpy as np
import pyvista as pv
import trimesh
from scipy.optimize import linear_sum_assignment

from utils.cfd_lumen.ultraliser_qc import _section_polygon, _triangle_intersections

from .io import BoundaryInput, SurfacePrepareError
from .local_cut import orthogonal_basis
from .mesh_quality import triangle_metrics


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
                "normal_abs_dot_outward": normal_dot,
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
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "vertex_count": len(mesh.vertices),
        "triangle_count": len(mesh.faces),
        "component_count": components,
        "assessment_component_count_excluding_official_raw_degenerates": assessment_components,
        "watertight": bool(mesh.is_watertight),
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
        vector = distal.center_um - proximal.center_um
        norm = float(np.linalg.norm(vector))
        actual_length = float(np.dot(vector, boundary.outward_normal))
        direction_dot = float(np.dot(vector / norm, boundary.outward_normal))
        length_error = abs(actual_length - boundary.extension_length_um) / boundary.extension_length_um
        area_error = abs(distal.area_um2 - proximal.area_um2) / proximal.area_um2
        checks = {
            "direction_dot_at_least_0_999": direction_dot >= 0.999,
            "relative_length_error_at_most_0_02": length_error <= 0.02,
            "distal_area_relative_error_at_most_0_05": area_error <= 0.05,
            "positive_distal_area": distal.area_um2 > 0.0,
            "interface_edges_detected": interfaces[boundary.port_id]["interface_edge_count"] > 0,
        }
        rows.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "role": boundary.role,
                "boundary_origin": boundary.boundary_origin,
                "planned_extension_length_um": boundary.extension_length_um,
                "actual_extension_length_um": actual_length,
                "extension_length_relative_error": length_error,
                "extension_direction_dot": direction_dot,
                "original_area_um2": proximal.area_um2,
                "distal_area_um2": distal.area_um2,
                "distal_area_relative_error": area_error,
                "equivalent_radius_original_um": proximal.equivalent_radius_um,
                "equivalent_radius_distal_um": distal.equivalent_radius_um,
                "distal_center_um": distal.center_um.tolist(),
                "distal_planarity_error_um": distal.planarity_error_um,
                **interfaces[boundary.port_id],
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


def tag_and_export_final_surface(
    capped_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    geometry_directory: Path,
    boundary_directory: Path,
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
    tagged_vtp = geometry_directory / "cfd_surface_vmtk_tps_remeshed_um.vtp"
    tagged_stl = geometry_directory / "cfd_surface_vmtk_tps_remeshed_um.stl"
    meter_stl = geometry_directory / "cfd_surface_vmtk_tps_remeshed_m.stl"
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
    report = {
        "status": "PASS",
        "mapping_method": "one-to-one nearest predicted extension end",
        "wall_entity_id": wall_entity_id,
        "distal_boundary_count": len(cap_ids),
        "boundaries": sorted(mapping_rows, key=lambda row: int(row["boundary_index"])),
    }
    return outputs, report


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
                "geometry_source": "final VMTK remeshed surface; 20 normal cross sections",
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
