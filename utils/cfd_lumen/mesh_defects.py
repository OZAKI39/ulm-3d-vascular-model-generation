"""Quantitative surface-defect and triangle-quality diagnostics for v2."""

from __future__ import annotations

from collections import defaultdict, deque
from typing import Any

import numpy as np
import trimesh
import vtk
from scipy.spatial import cKDTree


def _percentile(values: np.ndarray, q: float) -> float | None:
    return float(np.percentile(values, q)) if len(values) else None


def triangle_quality(mesh: trimesh.Trimesh, face_ids: np.ndarray | None = None) -> dict[str, Any]:
    """Return scale-independent triangle-shape metrics and per-face arrays."""

    ids = np.arange(len(mesh.faces), dtype=np.int64) if face_ids is None else np.asarray(face_ids)
    triangles = np.asarray(mesh.triangles[ids], dtype=float)
    if len(triangles) == 0:
        empty = np.empty(0, dtype=float)
        return {"face_ids": ids, "area": empty, "aspect_ratio": empty,
                "minimum_angle_deg": empty, "maximum_angle_deg": empty, "summary": {}}
    edges = np.stack(
        (
            np.linalg.norm(triangles[:, 1] - triangles[:, 0], axis=1),
            np.linalg.norm(triangles[:, 2] - triangles[:, 1], axis=1),
            np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1),
        ),
        axis=1,
    )
    area = np.asarray(mesh.area_faces[ids], dtype=float)
    longest = edges.max(axis=1)
    # 1.0 for an equilateral triangle; grows for slivers.
    aspect = np.divide(
        longest**2,
        4.0 * np.sqrt(3.0) * area,
        out=np.full_like(longest, np.inf),
        where=area > 0,
    )
    a, b, c = edges[:, 0], edges[:, 1], edges[:, 2]
    cosines = np.column_stack(
        (
            np.divide(b * b + c * c - a * a, 2.0 * b * c, out=np.ones_like(a), where=b * c > 0),
            np.divide(c * c + a * a - b * b, 2.0 * c * a, out=np.ones_like(a), where=c * a > 0),
            np.divide(a * a + b * b - c * c, 2.0 * a * b, out=np.ones_like(a), where=a * b > 0),
        )
    )
    angles = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
    minimum_angle = angles.min(axis=1)
    maximum_angle = angles.max(axis=1)
    finite_aspect = aspect[np.isfinite(aspect)]
    degenerate_count = int(np.count_nonzero(~np.isfinite(aspect) | (area <= 0.0)))
    summary = {
        "triangle_count": int(len(ids)),
        "triangle_area_min": float(area.min()),
        "triangle_area_p5": _percentile(area, 5),
        "triangle_area_median": _percentile(area, 50),
        "aspect_ratio_p95": _percentile(finite_aspect, 95),
        "aspect_ratio_max": float(finite_aspect.max()) if len(finite_aspect) else None,
        "degenerate_triangle_count": degenerate_count,
        "min_angle_p5_deg": _percentile(minimum_angle, 5),
        "minimum_angle_deg": float(minimum_angle.min()),
        "maximum_angle_deg": float(maximum_angle.max()),
    }
    return {
        "face_ids": ids,
        "area": area,
        "aspect_ratio": aspect,
        "minimum_angle_deg": minimum_angle,
        "maximum_angle_deg": maximum_angle,
        "summary": summary,
    }


def _edge_groups(edges: np.ndarray, vertices: np.ndarray) -> list[dict[str, Any]]:
    adjacency: dict[int, list[int]] = defaultdict(list)
    for edge_id, (first, second) in enumerate(edges):
        adjacency[int(first)].append(edge_id)
        adjacency[int(second)].append(edge_id)
    unseen = set(range(len(edges)))
    groups: list[dict[str, Any]] = []
    while unseen:
        seed = unseen.pop()
        queue = deque([seed])
        group = [seed]
        while queue:
            edge_id = queue.popleft()
            for vertex_id in edges[edge_id]:
                for neighbor in adjacency[int(vertex_id)]:
                    if neighbor in unseen:
                        unseen.remove(neighbor)
                        queue.append(neighbor)
                        group.append(neighbor)
        selected = edges[np.asarray(group, dtype=np.int64)]
        points = vertices[np.unique(selected)]
        length = np.linalg.norm(vertices[selected[:, 1]] - vertices[selected[:, 0]], axis=1).sum()
        degrees: dict[int, int] = defaultdict(int)
        for first, second in selected:
            degrees[int(first)] += 1
            degrees[int(second)] += 1
        groups.append(
            {
                "loop_id": len(groups),
                "edge_count": int(len(selected)),
                "loop_length_um": float(length),
                "centroid_um": np.mean(points, axis=0).tolist(),
                "closed_loop": bool(degrees and all(value == 2 for value in degrees.values())),
            }
        )
    return groups


def edge_defects(mesh: trimesh.Trimesh) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    inverse = np.asarray(mesh.edges_unique_inverse, dtype=np.int64)
    counts = np.bincount(inverse, minlength=len(mesh.edges_unique))
    boundary = np.asarray(mesh.edges_unique[counts == 1], dtype=np.int64)
    nonmanifold = np.asarray(mesh.edges_unique[counts > 2], dtype=np.int64)
    report = {
        "boundary_edge_count": int(len(boundary)),
        "boundary_loops": _edge_groups(boundary, np.asarray(mesh.vertices)) if len(boundary) else [],
        "non_manifold_edge_count": int(len(nonmanifold)),
    }
    return report, {"boundary_edges": boundary, "nonmanifold_edges": nonmanifold}


def _local_face_ids(
    mesh: trimesh.Trimesh,
    center: np.ndarray,
    radius: float,
) -> np.ndarray:
    triangles = np.asarray(mesh.triangles)
    minimum = center - radius
    maximum = center + radius
    overlaps = np.all(triangles.max(axis=1) >= minimum, axis=1) & np.all(
        triangles.min(axis=1) <= maximum, axis=1
    )
    return np.flatnonzero(overlaps)


def _triangle_intersections(
    mesh: trimesh.Trimesh,
    face_ids: np.ndarray,
) -> tuple[list[tuple[int, int]], int]:
    """Use an exact R-tree AABB index; never perform a global O(N^2) scan."""

    if len(face_ids) < 2:
        return [], 0
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces[face_ids], dtype=np.int64)
    pairs: list[tuple[int, int]] = []
    candidates_checked = 0
    triangle = vtk.vtkTriangle
    triangles = vertices[faces]
    triangle_minimum = triangles.min(axis=1)
    triangle_maximum = triangles.max(axis=1)
    bounds = np.column_stack((triangle_minimum, triangle_maximum))
    tree = trimesh.util.bounds_tree(bounds)
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
            other = triangles[other_id]
            if triangle.TrianglesIntersect(*current, *other):
                pairs.append((int(face_ids[local_id]), int(face_ids[other_id])))
    return pairs, candidates_checked


def _normal_diagnostics(mesh: trimesh.Trimesh) -> tuple[dict[str, Any], np.ndarray]:
    original = np.asarray(mesh.face_normals, dtype=float).copy()
    repaired = mesh.copy()
    trimesh.repair.fix_normals(repaired, multibody=True)
    recomputed = np.asarray(repaired.face_normals, dtype=float)
    flipped = np.flatnonzero(np.einsum("ij,ij->i", original, recomputed) < 0.0)
    return {
        "is_winding_consistent": bool(mesh.is_winding_consistent),
        "face_normal_consistency": float(1.0 - len(flipped) / max(1, len(mesh.faces))),
        "number_of_flipped_faces": int(len(flipped)),
    }, flipped


def _ray_cast_internal_faces(
    mesh: trimesh.Trimesh,
    face_ids: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Classify faces whose two normal-offset sides both lie inside the solid."""

    if len(face_ids) == 0:
        return np.empty(0, dtype=np.int64), {"ray_cast_faces_checked": 0, "ray_cast_error": None}
    triangles = np.asarray(mesh.triangles, dtype=float)
    local_edges = np.stack(
        (
            np.linalg.norm(triangles[:, 1] - triangles[:, 0], axis=1),
            np.linalg.norm(triangles[:, 2] - triangles[:, 1], axis=1),
            np.linalg.norm(triangles[:, 0] - triangles[:, 2], axis=1),
        ),
        axis=1,
    )
    # Scale each normal probe to its own triangle. A global median-edge offset
    # can jump across a nearby surface when testing marching-cubes slivers,
    # incorrectly classifying a valid exterior face as internal.
    face_epsilon = np.maximum(local_edges.min(axis=1) * 1.0e-2, 1.0e-9)
    internal: list[int] = []
    try:
        for start in range(0, len(face_ids), 512):
            selected = face_ids[start : start + 512]
            centers = np.asarray(mesh.triangles_center[selected])
            normals = np.asarray(mesh.face_normals[selected])
            epsilon = face_epsilon[selected, None]
            plus_inside = mesh.contains(centers + epsilon * normals)
            minus_inside = mesh.contains(centers - epsilon * normals)
            internal.extend(map(int, selected[plus_inside & minus_inside]))
        return np.asarray(internal, dtype=np.int64), {
            "ray_cast_faces_checked": int(len(face_ids)),
            "ray_cast_offset_method": "1% of each face minimum edge length",
            "ray_cast_offset_um_median": float(np.median(face_epsilon[face_ids])),
            "ray_cast_offset_um_min": float(np.min(face_epsilon[face_ids])),
            "ray_cast_offset_um_max": float(np.max(face_epsilon[face_ids])),
            "ray_cast_error": None,
        }
    except Exception as exc:
        return np.empty(0, dtype=np.int64), {
            "ray_cast_faces_checked": int(len(face_ids)),
            "ray_cast_offset_um": epsilon,
            "ray_cast_error": f"{type(exc).__name__}: {exc}",
        }


def diagnose_mesh_defects(
    mesh: trimesh.Trimesh,
    junctions: list[tuple[int, np.ndarray, float]],
    *,
    ray_face_ids: np.ndarray | None = None,
    ray_sample_limit: int | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Diagnose topology, internal/self intersections, normals, and sliver faces."""

    edge_report, artifacts = edge_defects(mesh)
    sorted_faces = np.sort(np.asarray(mesh.faces, dtype=np.int64), axis=1)
    _, duplicate_counts = np.unique(sorted_faces, axis=0, return_counts=True)
    duplicate_face_count = int(np.sum(np.maximum(duplicate_counts - 1, 0)))
    components = mesh.split(only_watertight=False)
    component_rows = [
        {
            "component_id": index,
            "area_um2": float(component.area),
            "volume_um3": float(abs(component.volume)) if component.is_watertight else None,
            "triangle_count": int(len(component.faces)),
        }
        for index, component in enumerate(components)
    ]
    global_quality = triangle_quality(mesh)
    junction_rows: list[dict[str, Any]] = []
    intersection_pairs: set[tuple[int, int]] = set()
    checked_candidates = 0
    local_face_union: set[int] = set()
    junction_bad_faces: set[int] = set()
    for node_id, center, radius in junctions:
        face_ids = _local_face_ids(mesh, np.asarray(center), max(3.0 * radius, 1.0e-6))
        local_face_union.update(map(int, face_ids))
        pairs, checked = _triangle_intersections(mesh, face_ids)
        intersection_pairs.update(tuple(sorted(pair)) for pair in pairs)
        checked_candidates += checked
        quality_result = triangle_quality(mesh, face_ids)
        bad_mask = (quality_result["aspect_ratio"] >= 20.0) | (
            quality_result["minimum_angle_deg"] <= 1.0
        )
        junction_bad_faces.update(map(int, face_ids[bad_mask]))
        quality = quality_result["summary"]
        junction_rows.append(
            {
                "junction_node_id": int(node_id),
                "local_face_count": int(len(face_ids)),
                "self_intersection_pair_count": int(len(pairs)),
                **quality,
            }
        )

    # Near-coincident opposing faces are an independent signal for retained internal seams.
    local_ids = np.asarray(sorted(local_face_union), dtype=np.int64)
    near_pairs: list[tuple[int, int]] = []
    if len(local_ids) > 1:
        centroids = np.asarray(mesh.triangles_center[local_ids])
        normals = np.asarray(mesh.face_normals[local_ids])
        edge_lengths = np.asarray(mesh.edges_unique_length)
        tolerance = max(float(np.median(edge_lengths)) * 1.0e-4, 1.0e-9)
        for first, second in cKDTree(centroids).query_pairs(tolerance):
            if float(np.dot(normals[first], normals[second])) < -0.95:
                near_pairs.append((int(local_ids[first]), int(local_ids[second])))

    ray_ids = local_ids
    if ray_sample_limit is not None and len(ray_ids) > ray_sample_limit:
        sampled = ray_ids[
            np.linspace(0, len(ray_ids) - 1, ray_sample_limit).round().astype(np.int64)
        ]
        mandatory = np.asarray(ray_face_ids if ray_face_ids is not None else [], dtype=np.int64)
        suspect = np.asarray(
            [face_id for pair in intersection_pairs.union(set(near_pairs)) for face_id in pair],
            dtype=np.int64,
        )
        ray_ids = np.unique(np.concatenate((sampled, mandatory, suspect)))
    elif ray_face_ids is not None:
        ray_ids = np.unique(np.concatenate((ray_ids, np.asarray(ray_face_ids, dtype=np.int64))))
    ray_internal, ray_report = _ray_cast_internal_faces(mesh, ray_ids)
    suspected_internal_ids = sorted(
        set(map(int, ray_internal))
        | {face_id for pair in near_pairs for face_id in pair}
    )
    normal_report, flipped_faces = _normal_diagnostics(mesh)
    intersection_rows: list[dict[str, Any]] = []
    for first, second in sorted(intersection_pairs):
        location = 0.5 * (
            np.asarray(mesh.triangles_center[first]) + np.asarray(mesh.triangles_center[second])
        )
        if junctions:
            nearest_node, nearest_center, _ = min(
                junctions, key=lambda item: float(np.linalg.norm(location - item[1]))
            )
            nearest_distance = float(np.linalg.norm(location - nearest_center))
        else:
            nearest_node, nearest_distance = None, None
        intersection_rows.append(
            {
                "face_id_a": int(first),
                "face_id_b": int(second),
                "location_um": location.tolist(),
                "nearest_junction_node_id": nearest_node,
                "distance_to_nearest_junction_um": nearest_distance,
            }
        )
    report = {
        **edge_report,
        "watertight": bool(mesh.is_watertight),
        "surface_connected_component_count": int(len(components)),
        "components": component_rows,
        "duplicate_face_count": duplicate_face_count,
        "suspected_internal_face_count": int(len(suspected_internal_ids)),
        "near_coincident_opposing_face_pair_count": int(len(near_pairs)),
        "self_intersection_count": int(len(intersection_pairs)),
        "self_intersection_pairs": intersection_rows,
        "self_intersection_candidate_pairs_checked": int(checked_candidates),
        "self_intersection_scope": "junction AABB via exact R-tree plus vtkTriangle",
        "internal_face_ray_scope": (
            "all junction AABB faces"
            if ray_sample_limit is None
            else f"mandatory suspects plus {ray_sample_limit} deterministic junction faces"
        ),
        **ray_report,
        **normal_report,
        "triangle_quality": global_quality["summary"],
        "junction_triangle_quality": junction_rows,
    }
    artifacts.update(
        {
            "self_intersection_pairs": np.asarray(sorted(intersection_pairs), dtype=np.int64).reshape((-1, 2)),
            "suspected_internal_faces": np.asarray(suspected_internal_ids, dtype=np.int64),
            "flipped_faces": flipped_faces,
            "aspect_ratio": global_quality["aspect_ratio"],
            "junction_bad_faces": np.asarray(sorted(junction_bad_faces), dtype=np.int64),
        }
    )
    return report, artifacts
