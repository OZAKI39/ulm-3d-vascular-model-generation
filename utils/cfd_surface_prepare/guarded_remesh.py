"""Topological guard assignment and exact diagnostics for guarded VMTK remeshing."""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv
import trimesh
from scipy.spatial import cKDTree

from utils.cfd_lumen.ultraliser_qc import _triangle_intersections

from .config import MeshQualityConfig
from .io import BoundaryInput, SurfacePrepareError
from .mesh_quality import summarize_extension_mesh, triangle_metrics
from .vmtk_qc import polydata_mesh, symmetric_mesh_size_mismatch


ENTITY_NAMES = {0: "CAP", 1: "CORE", 2: "PROXIMAL_GUARD", 3: "EXTENSION_BODY"}


def _faces(data: pv.PolyData) -> np.ndarray:
    return np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]


def _edge_face_map(faces: np.ndarray) -> dict[tuple[int, int], list[int]]:
    mapping: dict[tuple[int, int], list[int]] = defaultdict(list)
    for face_id, face in enumerate(np.asarray(faces, dtype=np.int64)):
        for first, second in zip(face, np.roll(face, -1)):
            mapping[tuple(sorted((int(first), int(second))))].append(face_id)
    return mapping


def extension_face_adjacency_layers(
    faces: np.ndarray, regions: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return extension BFS distance and layer-zero faces at the CORE interface."""

    faces = np.asarray(faces, dtype=np.int64)
    regions = np.asarray(regions, dtype=np.uint8)
    if len(faces) != len(regions):
        raise ValueError("faces and regions must have equal length")
    edge_faces = _edge_face_map(faces)
    extension_neighbors: dict[int, set[int]] = defaultdict(set)
    seeds: set[int] = set()
    for linked in edge_faces.values():
        if len(linked) != 2:
            continue
        first, second = linked
        values = {int(regions[first]), int(regions[second])}
        if values == {0, 1}:
            seeds.add(first if regions[first] == 1 else second)
        elif values == {1}:
            extension_neighbors[first].add(second)
            extension_neighbors[second].add(first)
    distances = np.full(len(faces), -1, dtype=np.int64)
    queue: deque[int] = deque()
    for face_id in sorted(seeds):
        distances[face_id] = 0
        queue.append(face_id)
    while queue:
        current = queue.popleft()
        for neighbor in extension_neighbors[current]:
            if distances[neighbor] >= 0:
                continue
            distances[neighbor] = distances[current] + 1
            queue.append(neighbor)
    return distances, np.asarray(sorted(seeds), dtype=np.int64)


def _boundary_assignment(
    centers: np.ndarray, boundaries: list[BoundaryInput]
) -> np.ndarray:
    scores: list[np.ndarray] = []
    for boundary in boundaries:
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
    return np.argmin(np.column_stack(scores), axis=1)


def assign_guarded_remesh_entities(
    raw_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    *,
    face_layers: int = 2,
    entity_array_name: str = "RemeshEntityId",
    core_entity_id: int = 1,
    guard_entity_id: int = 2,
    body_entity_id: int = 3,
) -> tuple[dict[str, Any], np.ndarray]:
    """Assign CORE/GUARD/BODY using exactly two extension adjacency layers."""

    if face_layers != 2:
        raise SurfacePrepareError("VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED:face_layers")
    data = pv.read(raw_vtp).triangulate()
    if "SurfaceRegionId" not in data.cell_data:
        raise SurfacePrepareError("VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED:regions_missing")
    faces = _faces(data)
    regions = np.asarray(data.cell_data["SurfaceRegionId"], dtype=np.uint8)
    if sorted(int(value) for value in np.unique(regions)) != [0, 1]:
        raise SurfacePrepareError("VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED:regions")
    distances, seeds = extension_face_adjacency_layers(faces, regions)
    core = regions == 0
    extension = regions == 1
    guard = extension & (distances >= 0) & (distances < face_layers)
    body = extension & ~guard
    entities = np.zeros(len(faces), dtype=np.int32)
    entities[core] = core_entity_id
    entities[guard] = guard_entity_id
    entities[body] = body_entity_id
    points_before = np.asarray(data.points).copy()
    faces_before = faces.copy()
    data.cell_data[entity_array_name] = entities
    data.save(raw_vtp, binary=True)
    saved = pv.read(raw_vtp).triangulate()
    saved_entities = np.asarray(saved.cell_data[entity_array_name], dtype=np.int32)

    boundary_list = list(boundaries)
    centers = np.asarray(saved.points)[faces].mean(axis=1)
    extension_ids = np.flatnonzero(extension)
    port_assignment = _boundary_assignment(centers[extension_ids], boundary_list)
    face_port = np.full(len(faces), -1, dtype=np.int32)
    face_port[extension_ids] = port_assignment
    edge_faces = _edge_face_map(faces)
    distal_face_ids = {
        linked[0]
        for linked in edge_faces.values()
        if len(linked) == 1 and extension[linked[0]]
    }
    per_port: list[dict[str, Any]] = []
    suspicious = False
    for local_index, boundary in enumerate(boundary_list):
        guard_ids = np.flatnonzero(guard & (face_port == local_index))
        body_ids = np.flatnonzero(body & (face_port == local_index))
        port_extension_count = len(guard_ids) + len(body_ids)
        guard_vertices = (
            np.unique(faces[guard_ids]) if len(guard_ids) else np.empty(0, dtype=np.int64)
        )
        if len(guard_ids):
            axial = (
                np.asarray(saved.points)[faces[guard_ids]].reshape((-1, 3))
                - boundary.center_um
            ) @ boundary.outward_normal
            axial_min = float(np.min(axial))
            axial_max = float(np.max(axial))
            axial_width = float(axial_max - axial_min)
        else:
            axial_min = axial_max = axial_width = float("nan")
        contains_distal = bool(set(map(int, guard_ids)) & distal_face_ids)
        guard_fraction = (
            float(len(guard_ids) / port_extension_count)
            if port_extension_count
            else 1.0
        )
        port_suspicious = bool(
            len(guard_ids) == 0
            or len(body_ids) == 0
            or contains_distal
            or guard_fraction >= 0.25
            or (np.isfinite(axial_width) and axial_width > 2.0 * boundary.source_radius_um)
        )
        suspicious |= port_suspicious
        per_port.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "guard_face_count": int(len(guard_ids)),
                "guard_vertex_count": int(len(guard_vertices)),
                "extension_body_face_count": int(len(body_ids)),
                "guard_axial_extent_min_um": axial_min,
                "guard_axial_extent_max_um": axial_max,
                "guard_maximum_axial_width_um": axial_width,
                "guard_width_in_source_radius": axial_width / boundary.source_radius_um,
                "guard_width_in_diameter": axial_width / (2.0 * boundary.source_radius_um),
                "guard_face_fraction_of_extension": guard_fraction,
                "guard_contains_distal_boundary": contains_distal,
                "classification_suspicious": port_suspicious,
            }
        )
    unknown = int(np.count_nonzero(entities == 0))
    checks = {
        "geometry_points_unchanged": bool(
            np.array_equal(np.asarray(saved.points), points_before)
        ),
        "geometry_connectivity_unchanged": bool(
            np.array_equal(_faces(saved), faces_before)
        ),
        "guard_layer_zero_detected": bool(len(seeds)),
        "guard_triangle_count_positive": bool(np.count_nonzero(guard)),
        "body_triangle_count_positive": bool(np.count_nonzero(body)),
        "unknown_face_count_zero": unknown == 0,
        "entity_ids_exact": sorted(int(value) for value in np.unique(saved_entities))
        == [core_entity_id, guard_entity_id, body_entity_id],
        "guard_excludes_core": not bool(np.any(guard & core)),
        "guard_excludes_distal_boundary": not any(
            row["guard_contains_distal_boundary"] for row in per_port
        ),
    }
    if not all(checks.values()):
        status = "VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED"
    elif suspicious:
        status = "GUARD_REGION_CLASSIFICATION_SUSPICIOUS"
    else:
        status = "PASS"
    return {
        "status": status,
        "checks": checks,
        "mode": "extension_face_adjacency_layers",
        "guard_face_layers": face_layers,
        "guard_layers_included": list(range(face_layers)),
        "core_entity_id": core_entity_id,
        "guard_entity_id": guard_entity_id,
        "extension_body_entity_id": body_entity_id,
        "entity_ids": sorted(int(value) for value in np.unique(saved_entities)),
        "core_face_count": int(np.count_nonzero(core)),
        "guard_face_count": int(np.count_nonzero(guard)),
        "body_face_count": int(np.count_nonzero(body)),
        "unknown_face_count": unknown,
        "interface_seed_face_count": int(len(seeds)),
        "per_port": per_port,
    }, distances


def _triangle_key(points: np.ndarray, face: np.ndarray) -> tuple[tuple[float, ...], ...]:
    return tuple(sorted(tuple(float(value) for value in point) for point in points[face]))


def _entity_vertex_ids(faces: np.ndarray, entities: np.ndarray, value: int) -> np.ndarray:
    selected = faces[entities == value]
    return np.unique(selected) if len(selected) else np.empty(0, dtype=np.int64)


def _vertex_motion(
    source: np.ndarray, target: np.ndarray, target_dtype: np.dtype[Any]
) -> dict[str, Any]:
    if not len(source) or not len(target):
        return {"max_um": float("inf"), "P95_um": float("inf"), "exact": False}
    distances, _ = cKDTree(np.asarray(target, dtype=float)).query(
        np.asarray(source, dtype=float), k=1
    )
    source_cast = np.asarray(source).astype(target_dtype, copy=False)
    source_set = Counter(tuple(float(v) for v in point) for point in source_cast)
    target_set = Counter(tuple(float(v) for v in point) for point in target)
    return {
        "max_um": float(np.max(distances)),
        "P95_um": float(np.percentile(distances, 95)),
        "exact": source_set == target_set,
    }


def _locked_entity_report(
    raw: pv.PolyData,
    remeshed: pv.PolyData,
    *,
    entity_array_name: str,
    entity_id: int,
    label: str,
) -> dict[str, Any]:
    raw_faces = _faces(raw)
    remeshed_faces = _faces(remeshed)
    raw_entities = np.asarray(raw.cell_data[entity_array_name], dtype=np.int32)
    remeshed_entities = np.asarray(
        remeshed.cell_data[entity_array_name], dtype=np.int32
    )
    raw_ids = np.flatnonzero(raw_entities == entity_id)
    remeshed_ids = np.flatnonzero(remeshed_entities == entity_id)
    raw_points = np.asarray(raw.points)
    remeshed_points = np.asarray(remeshed.points)
    cast_points = raw_points.astype(remeshed_points.dtype, copy=False)
    raw_counter = Counter(_triangle_key(cast_points, raw_faces[index]) for index in raw_ids)
    remeshed_counter = Counter(
        _triangle_key(remeshed_points, remeshed_faces[index]) for index in remeshed_ids
    )
    missing = int(sum((raw_counter - remeshed_counter).values()))
    added = int(sum((remeshed_counter - raw_counter).values()))
    raw_vertex_ids = _entity_vertex_ids(raw_faces, raw_entities, entity_id)
    remeshed_vertex_ids = _entity_vertex_ids(
        remeshed_faces, remeshed_entities, entity_id
    )
    motion = _vertex_motion(
        raw_points[raw_vertex_ids],
        remeshed_points[remeshed_vertex_ids],
        remeshed_points.dtype,
    )
    tolerance = max(
        float(np.max(np.abs(np.spacing(remeshed_points)))) * np.sqrt(3.0),
        np.finfo(float).eps,
    )
    checks = {
        f"{label}_face_count_unchanged": len(raw_ids) == len(remeshed_ids),
        f"{label}_vertex_count_unchanged": len(raw_vertex_ids)
        == len(remeshed_vertex_ids),
        f"{label}_motion_within_output_machine_precision": motion["max_um"]
        <= tolerance,
        f"{label}_vertices_exact_after_output_dtype_cast": motion["exact"],
        f"{label}_connectivity_unchanged": missing == 0 and added == 0,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        f"raw_{label}_face_count": int(len(raw_ids)),
        f"remeshed_{label}_face_count": int(len(remeshed_ids)),
        f"raw_{label}_vertex_count": int(len(raw_vertex_ids)),
        f"remeshed_{label}_vertex_count": int(len(remeshed_vertex_ids)),
        f"{label}_vertex_max_motion_um": motion["max_um"],
        f"{label}_vertex_P95_motion_um": motion["P95_um"],
        f"{label}_connectivity_changed_count": missing + added,
        f"{label}_missing_face_count": missing,
        f"{label}_added_face_count": added,
        "machine_precision_tolerance_um": tolerance,
    }


def _entity_boundary_report(
    raw: pv.PolyData,
    remeshed: pv.PolyData,
    *,
    entity_array_name: str,
    first_id: int,
    second_id: int,
    label: str,
) -> dict[str, Any]:
    raw_faces = _faces(raw)
    remeshed_faces = _faces(remeshed)
    raw_entities = np.asarray(raw.cell_data[entity_array_name], dtype=np.int32)
    remeshed_entities = np.asarray(
        remeshed.cell_data[entity_array_name], dtype=np.int32
    )
    raw_shared = np.intersect1d(
        _entity_vertex_ids(raw_faces, raw_entities, first_id),
        _entity_vertex_ids(raw_faces, raw_entities, second_id),
    )
    remeshed_shared = np.intersect1d(
        _entity_vertex_ids(remeshed_faces, remeshed_entities, first_id),
        _entity_vertex_ids(remeshed_faces, remeshed_entities, second_id),
    )
    motion = _vertex_motion(
        np.asarray(raw.points)[raw_shared],
        np.asarray(remeshed.points)[remeshed_shared],
        np.asarray(remeshed.points).dtype,
    )
    tolerance = max(
        float(np.max(np.abs(np.spacing(np.asarray(remeshed.points))))) * np.sqrt(3.0),
        np.finfo(float).eps,
    )
    passed = (
        len(raw_shared) == len(remeshed_shared)
        and motion["exact"]
        and motion["max_um"] <= tolerance
    )
    return {
        "status": "PASS" if passed else "FAIL",
        f"{label}_shared_vertex_count": int(len(raw_shared)),
        f"remeshed_{label}_shared_vertex_count": int(len(remeshed_shared)),
        f"{label}_max_motion_um": motion["max_um"],
        f"{label}_P95_motion_um": motion["P95_um"],
        f"{label}_exact_after_output_dtype_cast": motion["exact"],
        "machine_precision_tolerance_um": tolerance,
    }


def guarded_entity_preservation_qc(
    raw_vtp: Path,
    remeshed_vtp: Path,
    *,
    entity_array_name: str = "RemeshEntityId",
    core_entity_id: int = 1,
    guard_entity_id: int = 2,
    body_entity_id: int = 3,
    restore_region_arrays: bool = True,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Prove CORE/GUARD locks, both interfaces, and BODY remesh effect."""

    raw = pv.read(raw_vtp).triangulate()
    remeshed = pv.read(remeshed_vtp).triangulate()
    expected = [core_entity_id, guard_entity_id, body_entity_id]
    for data in (raw, remeshed):
        if entity_array_name not in data.cell_data:
            raise SurfacePrepareError("VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED:array_missing")
        if sorted(int(v) for v in np.unique(data.cell_data[entity_array_name])) != expected:
            raise SurfacePrepareError("VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED:ids_changed")
    core = _locked_entity_report(
        raw,
        remeshed,
        entity_array_name=entity_array_name,
        entity_id=core_entity_id,
        label="core",
    )
    guard = _locked_entity_report(
        raw,
        remeshed,
        entity_array_name=entity_array_name,
        entity_id=guard_entity_id,
        label="guard",
    )
    core_guard = _entity_boundary_report(
        raw,
        remeshed,
        entity_array_name=entity_array_name,
        first_id=core_entity_id,
        second_id=guard_entity_id,
        label="core_guard",
    )
    guard_body = _entity_boundary_report(
        raw,
        remeshed,
        entity_array_name=entity_array_name,
        first_id=guard_entity_id,
        second_id=body_entity_id,
        label="guard_body",
    )
    boundary = {
        "status": (
            "PASS"
            if core_guard["status"] == guard_body["status"] == "PASS"
            else "FAIL"
        ),
        **{key: value for key, value in core_guard.items() if key != "status"},
        **{key: value for key, value in guard_body.items() if key != "status"},
    }
    raw_faces = _faces(raw)
    remeshed_faces = _faces(remeshed)
    raw_entities = np.asarray(raw.cell_data[entity_array_name], dtype=np.int32)
    remeshed_entities = np.asarray(
        remeshed.cell_data[entity_array_name], dtype=np.int32
    )
    raw_body_ids = np.flatnonzero(raw_entities == body_entity_id)
    remeshed_body_ids = np.flatnonzero(remeshed_entities == body_entity_id)
    raw_points = np.asarray(raw.points).astype(np.asarray(remeshed.points).dtype, copy=False)
    remeshed_points = np.asarray(remeshed.points)
    raw_counter = Counter(
        _triangle_key(raw_points, raw_faces[index]) for index in raw_body_ids
    )
    remeshed_counter = Counter(
        _triangle_key(remeshed_points, remeshed_faces[index])
        for index in remeshed_body_ids
    )
    raw_vertices = _entity_vertex_ids(raw_faces, raw_entities, body_entity_id)
    remeshed_vertices = _entity_vertex_ids(
        remeshed_faces, remeshed_entities, body_entity_id
    )
    changed = raw_counter != remeshed_counter
    effect = bool(
        len(raw_body_ids) != len(remeshed_body_ids)
        or len(raw_vertices) != len(remeshed_vertices)
        or changed
    )
    body = {
        "status": "PASS" if effect else "FAIL",
        "raw_body_face_count": int(len(raw_body_ids)),
        "remeshed_body_face_count": int(len(remeshed_body_ids)),
        "raw_body_vertex_count": int(len(raw_vertices)),
        "remeshed_body_vertex_count": int(len(remeshed_vertices)),
        "body_connectivity_changed": bool(changed),
        "body_remesh_effect_detected": effect,
    }
    if restore_region_arrays:
        output_entities = np.asarray(
            remeshed.cell_data[entity_array_name], dtype=np.int32
        )
        remeshed.cell_data["SurfaceRegionId"] = np.where(
            output_entities == core_entity_id, 0, 1
        ).astype(np.uint8)
        remeshed.cell_data["SurfaceRegion"] = np.where(
            output_entities == core_entity_id, "CORE", "EXTENSION"
        )
        remeshed.save(remeshed_vtp, binary=True)
    return core, guard, boundary, body


def _plane_section_points(
    triangle: np.ndarray, plane_point: np.ndarray, plane_normal: np.ndarray
) -> np.ndarray:
    signed = (triangle - plane_point) @ plane_normal
    points: list[np.ndarray] = []
    for index in range(3):
        following = (index + 1) % 3
        first = float(signed[index])
        second = float(signed[following])
        if first == 0.0:
            points.append(triangle[index])
        if first * second < 0.0:
            fraction = first / (first - second)
            points.append(
                triangle[index]
                + fraction * (triangle[following] - triangle[index])
            )
    unique: list[np.ndarray] = []
    for point in points:
        if not any(np.linalg.norm(point - prior) <= 1.0e-12 for prior in unique):
            unique.append(point)
    return np.asarray(unique, dtype=float)


def triangle_pair_intersection_diagnosis(
    vertices: np.ndarray,
    faces: np.ndarray,
    first_face_id: int,
    second_face_id: int,
) -> dict[str, Any]:
    """Distinguish legal mesh adjacency, point contact, and true penetration."""

    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    first_face = faces[first_face_id]
    second_face = faces[second_face_id]
    shared_ids = np.intersect1d(first_face, second_face)
    first_edges = {
        tuple(sorted((int(a), int(b))))
        for a, b in zip(first_face, np.roll(first_face, -1))
    }
    second_edges = {
        tuple(sorted((int(a), int(b))))
        for a, b in zip(second_face, np.roll(second_face, -1))
    }
    shared_edges = sorted(first_edges & second_edges)
    first = vertices[first_face]
    second = vertices[second_face]
    first_cross = np.cross(first[1] - first[0], first[2] - first[0])
    second_cross = np.cross(second[1] - second[0], second[2] - second[0])
    first_normal = first_cross / np.linalg.norm(first_cross)
    second_normal = second_cross / np.linalg.norm(second_cross)
    direction = np.cross(first_normal, second_normal)
    direction_norm = float(np.linalg.norm(direction))
    segment_start: np.ndarray | None = None
    segment_end: np.ndarray | None = None
    segment_length = 0.0
    intersection_type = "NO_INTERSECTION"
    if shared_edges:
        intersection_type = "LEGAL_SHARED_EDGE_CONTACT"
    elif len(shared_ids):
        intersection_type = "LEGAL_SHARED_VERTEX_CONTACT"
    elif direction_norm <= np.finfo(float).eps:
        intersection_type = "COPLANAR_REVIEW_REQUIRED"
    else:
        unit_direction = direction / direction_norm
        first_section = _plane_section_points(first, second[0], second_normal)
        second_section = _plane_section_points(second, first[0], first_normal)
        if len(first_section) >= 2 and len(second_section) >= 2:
            first_interval = np.sort(first_section @ unit_direction)[[0, -1]]
            second_interval = np.sort(second_section @ unit_direction)[[0, -1]]
            lower = float(max(first_interval[0], second_interval[0]))
            upper = float(min(first_interval[1], second_interval[1]))
            if upper >= lower:
                first_offset = float(first_normal @ first[0])
                second_offset = float(second_normal @ second[0])
                line_direction = np.cross(first_normal, second_normal)
                line_point = np.cross(
                    first_offset * second_normal - second_offset * first_normal,
                    line_direction,
                ) / float(line_direction @ line_direction)
                line_projection = float(line_point @ unit_direction)
                segment_start = line_point + (lower - line_projection) * unit_direction
                segment_end = line_point + (upper - line_projection) * unit_direction
                segment_length = max(0.0, upper - lower)
                intersection_type = (
                    "NONCOPLANAR_SEGMENT_PENETRATION"
                    if segment_length > 0.0
                    else "POINT_CONTACT"
                )
    true_intersection = bool(
        not shared_edges
        and len(shared_ids) == 0
        and intersection_type == "NONCOPLANAR_SEGMENT_PENETRATION"
        and segment_length > 0.0
    )
    return {
        "first_face_id": int(first_face_id),
        "second_face_id": int(second_face_id),
        "first_triangle_vertices_xyz": first.tolist(),
        "second_triangle_vertices_xyz": second.tolist(),
        "shared_vertex_count": int(len(shared_ids)),
        "shared_edge_count": int(len(shared_edges)),
        "shared_vertex_ids": [int(value) for value in shared_ids],
        "shared_edges": [list(edge) for edge in shared_edges],
        "first_normal": first_normal.tolist(),
        "second_normal": second_normal.tolist(),
        "first_triangle_centroid": first.mean(axis=0).tolist(),
        "second_triangle_centroid": second.mean(axis=0).tolist(),
        "minimum_triangle_distance_um": 0.0 if true_intersection else None,
        "true_triangle_triangle_intersection": true_intersection,
        "intersection_type": intersection_type,
        "intersection_segment_start_xyz": (
            segment_start.tolist() if segment_start is not None else None
        ),
        "intersection_segment_end_xyz": (
            segment_end.tolist() if segment_end is not None else None
        ),
        "intersection_segment_length_um": segment_length,
    }


def guarded_intersection_qc(
    surface_vtp: Path,
    *,
    entity_array_name: str = "RemeshEntityId",
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Classify every true non-adjacent intersection by guarded entity pair."""

    data, mesh = polydata_mesh(surface_vtp)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    entities = np.asarray(data.cell_data[entity_array_name], dtype=np.int32)
    pairs, candidates = _triangle_intersections(
        mesh, np.arange(len(faces), dtype=np.int64)
    )
    records: list[dict[str, Any]] = []
    counts = {
        "CORE_CORE": 0,
        "CORE_PROXIMAL_GUARD": 0,
        "CORE_EXTENSION_BODY": 0,
        "PROXIMAL_GUARD_PROXIMAL_GUARD": 0,
        "PROXIMAL_GUARD_EXTENSION_BODY": 0,
        "EXTENSION_BODY_EXTENSION_BODY": 0,
    }
    for first, second in pairs:
        diagnosis = triangle_pair_intersection_diagnosis(
            np.asarray(mesh.vertices), faces, first, second
        )
        first_name = ENTITY_NAMES[int(entities[first])]
        second_name = ENTITY_NAMES[int(entities[second])]
        order = {"CAP": 0, "CORE": 1, "PROXIMAL_GUARD": 2, "EXTENSION_BODY": 3}
        ordered_names = sorted((first_name, second_name), key=order.__getitem__)
        canonical = "_".join(ordered_names)
        diagnosis.update(
            {
                "first_entity_id": int(entities[first]),
                "second_entity_id": int(entities[second]),
                "entity_pair": canonical,
            }
        )
        if diagnosis["true_triangle_triangle_intersection"]:
            counts.setdefault(canonical, 0)
            counts[canonical] += 1
            records.append(diagnosis)
    total = int(sum(counts.values()))
    return {
        "status": "PASS" if total == 0 else "FAIL",
        "true_self_intersection_count": total,
        "vtk_candidate_intersection_count": int(len(pairs)),
        "candidate_pairs_checked": int(candidates),
        "classification_counts": counts,
        "intersections": records,
    }, records


def _point_to_segments_distance(point: np.ndarray, segments: np.ndarray) -> float:
    starts = segments[:, 0]
    vectors = segments[:, 1] - starts
    denominator = np.sum(vectors * vectors, axis=1)
    fraction = np.divide(
        np.sum((point - starts) * vectors, axis=1),
        denominator,
        out=np.zeros(len(segments), dtype=float),
        where=denominator > 0.0,
    )
    fraction = np.clip(fraction, 0.0, 1.0)
    closest = starts + fraction[:, None] * vectors
    return float(np.min(np.linalg.norm(closest - point, axis=1)))


def diagnose_previous_entityremesh(
    raw_vtp: Path,
    remeshed_vtp: Path,
    boundaries: Iterable[BoundaryInput],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Recompute the saved 180737 collision and locate its RAW BFS layer."""

    raw_data, raw_mesh = polydata_mesh(raw_vtp)
    remeshed_data, remeshed_mesh = polydata_mesh(remeshed_vtp)
    pairs, candidates = _triangle_intersections(
        remeshed_mesh,
        np.arange(len(remeshed_mesh.faces), dtype=np.int64),
    )
    entities = np.asarray(remeshed_data.cell_data["RemeshEntityId"], dtype=np.int32)
    target_pair: tuple[int, int] | None = None
    target_diagnosis: dict[str, Any] | None = None
    for first, second in pairs:
        if {int(entities[first]), int(entities[second])} != {1, 2}:
            continue
        diagnosis = triangle_pair_intersection_diagnosis(
            np.asarray(remeshed_mesh.vertices),
            np.asarray(remeshed_mesh.faces),
            first,
            second,
        )
        if diagnosis["true_triangle_triangle_intersection"]:
            target_pair = (first, second)
            target_diagnosis = diagnosis
            break
    if target_pair is None or target_diagnosis is None:
        return {
            "status": "ENTITY_REMESH_INTERSECTION_DETECTOR_REVIEW_REQUIRED",
            "REAL_GEOMETRIC_PENETRATION": "NO",
            "vtk_intersection_pair_count": int(len(pairs)),
            "candidate_pairs_checked": int(candidates),
        }, {"status": "NOT_APPLICABLE"}
    core_face_id = (
        target_pair[0] if entities[target_pair[0]] == 1 else target_pair[1]
    )
    extension_face_id = (
        target_pair[0] if entities[target_pair[0]] == 2 else target_pair[1]
    )
    diagnosis = {
        "status": "PASS",
        "REAL_GEOMETRIC_PENETRATION": "YES",
        "core_face_id": int(core_face_id),
        "extension_face_id": int(extension_face_id),
        "core_triangle_vertices_xyz": target_diagnosis[
            "first_triangle_vertices_xyz"
            if target_pair[0] == core_face_id
            else "second_triangle_vertices_xyz"
        ],
        "extension_triangle_vertices_xyz": target_diagnosis[
            "first_triangle_vertices_xyz"
            if target_pair[0] == extension_face_id
            else "second_triangle_vertices_xyz"
        ],
        "shared_vertex_count": target_diagnosis["shared_vertex_count"],
        "shared_edge_count": target_diagnosis["shared_edge_count"],
        "core_normal": target_diagnosis[
            "first_normal" if target_pair[0] == core_face_id else "second_normal"
        ],
        "extension_normal": target_diagnosis[
            "first_normal" if target_pair[0] == extension_face_id else "second_normal"
        ],
        "triangle_centroids": {
            "core": target_diagnosis[
                "first_triangle_centroid"
                if target_pair[0] == core_face_id
                else "second_triangle_centroid"
            ],
            "extension": target_diagnosis[
                "first_triangle_centroid"
                if target_pair[0] == extension_face_id
                else "second_triangle_centroid"
            ],
        },
        "minimum_triangle_distance_um": target_diagnosis[
            "minimum_triangle_distance_um"
        ],
        "true_triangle_triangle_intersection": True,
        "intersection_type": target_diagnosis["intersection_type"],
        "intersection_segment_start_xyz": target_diagnosis[
            "intersection_segment_start_xyz"
        ],
        "intersection_segment_end_xyz": target_diagnosis[
            "intersection_segment_end_xyz"
        ],
        "intersection_segment_length_um": target_diagnosis[
            "intersection_segment_length_um"
        ],
        "candidate_pairs_checked": int(candidates),
    }

    raw_faces = np.asarray(raw_mesh.faces, dtype=np.int64)
    raw_regions = np.asarray(raw_data.cell_data["SurfaceRegionId"], dtype=np.uint8)
    distances, seeds = extension_face_adjacency_layers(raw_faces, raw_regions)
    raw_extension_ids = np.flatnonzero(raw_regions == 1)
    failed_triangle = np.asarray(remeshed_mesh.vertices)[
        np.asarray(remeshed_mesh.faces)[extension_face_id]
    ]
    failed_center = failed_triangle.mean(axis=0)
    raw_extension = trimesh.Trimesh(
        vertices=np.asarray(raw_mesh.vertices),
        faces=raw_faces[raw_extension_ids],
        process=False,
    )
    _, nearest_distance, nearest_local = trimesh.proximity.closest_point(
        raw_extension, failed_center[None, :]
    )
    nearest_raw_face_id = int(raw_extension_ids[int(nearest_local[0])])
    edge_faces = _edge_face_map(raw_faces)
    interface_edges = np.asarray(
        [
            edge
            for edge, linked in edge_faces.items()
            if len(linked) == 2
            and {int(raw_regions[linked[0]]), int(raw_regions[linked[1]])} == {0, 1}
        ],
        dtype=np.int64,
    )
    interface_segments = np.asarray(raw_mesh.vertices)[interface_edges]
    spatial_distance = min(
        _point_to_segments_distance(point, interface_segments)
        for point in np.vstack((failed_triangle, failed_center))
    )
    boundary_list = list(boundaries)
    port_local = int(_boundary_assignment(failed_center[None, :], boundary_list)[0])
    boundary = boundary_list[port_local]
    axial = float((failed_center - boundary.center_um) @ boundary.outward_normal)
    layer_report = {
        "status": "PASS",
        "port_id": boundary.port_id,
        "boundary_index": boundary.index,
        "expected_port_id_suffix": "cut_002",
        "port_is_cut_002": boundary.port_id.endswith("cut_002"),
        "extension_face_id_previous_remesh": int(extension_face_id),
        "nearest_raw_extension_face_id": nearest_raw_face_id,
        "nearest_raw_triangle_distance_um": float(nearest_distance[0]),
        "axial_distance_from_original_cut_plane_um": axial,
        "axial_distance_in_source_radius": axial / boundary.source_radius_um,
        "axial_distance_in_diameter": axial / (2.0 * boundary.source_radius_um),
        "distance_to_CORE_EXTENSION_interface_um": spatial_distance,
        "face_adjacency_layer": int(distances[nearest_raw_face_id]),
        "interface_seed_face_count": int(len(seeds)),
        "raw_extension_face_count": int(len(raw_extension_ids)),
    }
    return diagnosis, layer_report


def guarded_region_mesh_quality(
    surface_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    *,
    entity_id: int,
    entity_label: str,
    local_target_edge_um: dict[int, float],
    quality: MeshQualityConfig,
    hard_body_gate: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Report GUARD and BODY separately and return spatial tail records."""

    data, mesh = polydata_mesh(surface_vtp)
    entities = np.asarray(data.cell_data["RemeshEntityId"], dtype=np.int32)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    selected_ids = np.flatnonzero(entities == entity_id)
    if len(selected_ids) == 0:
        raise SurfacePrepareError("VMTK_GUARD_ENTITY_ASSIGNMENT_FAILED:empty_entity")
    centers = np.asarray(mesh.triangles_center)[selected_ids]
    boundary_list = list(boundaries)
    assignment = _boundary_assignment(centers, boundary_list)
    raw_edge_faces = _edge_face_map(faces)
    interface_edges = np.asarray(
        [
            edge
            for edge, linked in raw_edge_faces.items()
            if len(linked) == 2
            and {int(entities[linked[0]]), int(entities[linked[1]])} == {1, 2}
        ],
        dtype=np.int64,
    )
    interface_segments = (
        np.asarray(mesh.vertices)[interface_edges]
        if len(interface_edges)
        else np.empty((0, 2, 3), dtype=float)
    )
    rows: list[dict[str, Any]] = []
    tail_records: list[dict[str, Any]] = []
    for local_index, boundary in enumerate(boundary_list):
        face_ids = selected_ids[assignment == local_index]
        selected_faces = faces[face_ids]
        metrics = triangle_metrics(np.asarray(mesh.vertices), selected_faces)
        target = float(local_target_edge_um[boundary.index])
        summary = summarize_extension_mesh(
            np.asarray(mesh.vertices),
            selected_faces,
            target_edge_length_um=target,
            local_original_median_edge_length_um=target,
            quality=quality,
        )
        size_ratio = float(summary["edge_length_median_um"] / target)
        mismatch = symmetric_mesh_size_mismatch(size_ratio)
        angle_tail = metrics.minimum_angles_deg < 5.0
        aspect_tail = metrics.aspect_ratios > 20.0
        tail = angle_tail | aspect_tail
        finite = all(
            np.isfinite(value)
            for value in summary.values()
            if isinstance(value, float)
        )
        checks = {"finite_metrics": bool(finite)}
        if hard_body_gate:
            checks.update(
                {
                    "bad_triangle_fraction": summary["bad_triangle_fraction"]
                    <= quality.maximum_bad_triangle_fraction,
                    "neighbor_area_ratio_p95": summary["neighbor_area_ratio_p95"]
                    <= quality.maximum_neighbor_area_ratio,
                    "symmetric_mesh_size_mismatch": mismatch <= 1.5,
                }
            )
        row = {
            "boundary_index": boundary.index,
            "port_id": boundary.port_id,
            "entity": entity_label,
            "triangle_count": int(len(selected_faces)),
            "minimum_angle_deg": summary["minimum_triangle_angle_deg"],
            "angle_P05_deg": summary["triangle_angle_p05_deg"],
            "aspect_ratio_median": summary["aspect_ratio_median"],
            "aspect_ratio_P95": summary["aspect_ratio_p95"],
            "aspect_ratio_max": summary["aspect_ratio_max"],
            "edge_length_median_um": summary["edge_length_median_um"],
            "edge_length_P95_um": summary["edge_length_p95_um"],
            "neighbor_area_ratio_P95": summary["neighbor_area_ratio_p95"],
            "bad_triangle_fraction": summary["bad_triangle_fraction"],
            "local_original_target_edge_um": target,
            "mesh_size_ratio": size_ratio,
            "symmetric_mesh_size_mismatch": mismatch,
            "triangle_count_angle_below_5deg": int(np.count_nonzero(angle_tail)),
            "triangle_count_aspect_above_20": int(np.count_nonzero(aspect_tail)),
            "checks": checks,
            "status": "PASS" if all(checks.values()) else "FAIL",
        }
        rows.append(row)
        for local_face_id in np.flatnonzero(tail):
            face_id = int(face_ids[local_face_id])
            center = np.asarray(mesh.triangles_center)[face_id]
            axial = float((center - boundary.center_um) @ boundary.outward_normal)
            distance = (
                _point_to_segments_distance(center, interface_segments)
                if len(interface_segments)
                else float("nan")
            )
            tail_records.append(
                {
                    "port_id": boundary.port_id,
                    "entity": entity_label,
                    "triangle_id": face_id,
                    "minimum_angle_deg": float(metrics.minimum_angles_deg[local_face_id]),
                    "aspect_ratio": float(metrics.aspect_ratios[local_face_id]),
                    "distance_to_original_interface_um": distance,
                    "axial_distance_um": axial,
                    "axial_distance_in_D": axial
                    / (2.0 * boundary.source_radius_um),
                }
            )
    report = {
        "status": (
            "PASS"
            if len(rows) == len(boundary_list)
            and all(row["status"] == "PASS" for row in rows)
            else "FAIL"
        ),
        "entity": entity_label,
        "hard_gate": hard_body_gate,
        "guard_sliver_policy": (
            "DIAGNOSTIC_ONLY" if not hard_body_gate else "BODY_HARD_BULK_QC"
        ),
        "maximum_symmetric_mesh_size_mismatch": 1.5,
        "boundaries": rows,
        "tail_triangle_count": len(tail_records),
    }
    return report, tail_records


def previous_tail_fraction_inside_guard(
    previous_remeshed_vtp: Path,
    guarded_raw_vtp: Path,
) -> dict[str, Any]:
    """Map previous sliver-tail centers to nearest RAW extension faces."""

    previous_data, previous_mesh = polydata_mesh(previous_remeshed_vtp)
    raw_data, raw_mesh = polydata_mesh(guarded_raw_vtp)
    previous_entities = np.asarray(
        previous_data.cell_data["RemeshEntityId"], dtype=np.int32
    )
    previous_ids = np.flatnonzero(previous_entities == 2)
    previous_faces = np.asarray(previous_mesh.faces)[previous_ids]
    metrics = triangle_metrics(np.asarray(previous_mesh.vertices), previous_faces)
    tail = (metrics.minimum_angles_deg < 5.0) | (metrics.aspect_ratios > 20.0)
    tail_ids = previous_ids[tail]
    raw_entities = np.asarray(raw_data.cell_data["RemeshEntityId"], dtype=np.int32)
    raw_extension_ids = np.flatnonzero(np.isin(raw_entities, (2, 3)))
    tree = cKDTree(np.asarray(raw_mesh.triangles_center)[raw_extension_ids])
    _, nearest = tree.query(np.asarray(previous_mesh.triangles_center)[tail_ids], k=1)
    nearest_raw = raw_extension_ids[np.asarray(nearest, dtype=np.int64)]
    inside = raw_entities[nearest_raw] == 2
    fraction = float(np.mean(inside)) if len(inside) else 0.0
    return {
        "previous_tail_triangle_count": int(len(tail_ids)),
        "previous_tail_inside_new_guard_count": int(np.count_nonzero(inside)),
        "previous_tail_fraction_within_new_guard": fraction,
        "PROXIMAL_TAIL_LOCALIZATION_CONFIRMED": "YES" if fraction > 0.5 else "NO",
        "role": "DIAGNOSTIC_ONLY",
    }
