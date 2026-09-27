"""Exact triangle clipping restricted to one boundary-local surgery cylinder."""

from __future__ import annotations

from collections import defaultdict, deque
from typing import Any

import numpy as np
from shapely.geometry import LineString, Polygon

from .io import BoundaryInput, SurfacePrepareError
from .types import CutLoop, TaggedSurface


def orthogonal_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    normal = np.asarray(normal, dtype=float) / np.linalg.norm(normal)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(helper, normal))) > 0.85:
        helper = np.array([0.0, 1.0, 0.0])
    first = np.cross(normal, helper)
    first /= np.linalg.norm(first)
    second = np.cross(normal, first)
    return first, second


def polygon_metrics(
    points: np.ndarray, normal: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray]:
    first, second = orthogonal_basis(normal)
    origin = np.mean(points, axis=0)
    projected = np.column_stack(((points - origin) @ first, (points - origin) @ second))
    polygon = Polygon(projected)
    if not polygon.is_valid or not polygon.exterior.is_simple or polygon.area <= 0:
        raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
    center_2d = np.asarray(polygon.centroid.coords[0], dtype=float)
    center = origin + center_2d[0] * first + center_2d[1] * second
    return float(polygon.area), center, projected


def _candidate_components(faces: np.ndarray, candidate_ids: np.ndarray) -> list[set[int]]:
    candidate = set(map(int, candidate_ids))
    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for face_id in candidate:
        face = faces[face_id]
        for first, second in zip(face, np.roll(face, -1)):
            edge_faces[tuple(sorted((int(first), int(second))))].append(face_id)
    adjacency: dict[int, set[int]] = {face_id: set() for face_id in candidate}
    for linked in edge_faces.values():
        for face_id in linked:
            adjacency[face_id].update(other for other in linked if other != face_id)
    components: list[set[int]] = []
    remaining = set(candidate)
    while remaining:
        start = remaining.pop()
        component = {start}
        queue = deque((start,))
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        components.append(component)
    return components


def _intersection_vertex(
    vertices: list[np.ndarray],
    cache: dict[tuple[int, int], int],
    first_id: int,
    second_id: int,
    first_distance: float,
    second_distance: float,
) -> int:
    if abs(first_distance) <= 1.0e-12:
        return first_id
    if abs(second_distance) <= 1.0e-12:
        return second_id
    key = tuple(sorted((first_id, second_id)))
    if key in cache:
        return cache[key]
    fraction = first_distance / (first_distance - second_distance)
    point = vertices[first_id] + fraction * (vertices[second_id] - vertices[first_id])
    vertex_id = len(vertices)
    vertices.append(np.asarray(point, dtype=float))
    cache[key] = vertex_id
    return vertex_id


def _clip_face(
    face: np.ndarray,
    distances: np.ndarray,
    vertices: list[np.ndarray],
    cache: dict[tuple[int, int], int],
) -> tuple[list[int], tuple[int, int]]:
    output: list[int] = []
    previous = int(face[-1])
    previous_distance = float(distances[previous])
    previous_inside = previous_distance <= 1.0e-12
    for current_raw in face:
        current = int(current_raw)
        current_distance = float(distances[current])
        current_inside = current_distance <= 1.0e-12
        if current_inside:
            if not previous_inside:
                output.append(
                    _intersection_vertex(
                        vertices,
                        cache,
                        previous,
                        current,
                        previous_distance,
                        current_distance,
                    )
                )
            output.append(current)
        elif previous_inside:
            output.append(
                _intersection_vertex(
                    vertices,
                    cache,
                    previous,
                    current,
                    previous_distance,
                    current_distance,
                )
            )
        previous = current
        previous_distance = current_distance
        previous_inside = current_inside
    cleaned: list[int] = []
    for vertex_id in output:
        if not cleaned or cleaned[-1] != vertex_id:
            cleaned.append(vertex_id)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    new_ids = [vertex_id for vertex_id in cleaned if vertex_id >= len(distances)]
    on_plane = [
        vertex_id
        for vertex_id in cleaned
        if vertex_id < len(distances) and abs(distances[vertex_id]) <= 1.0e-12
    ]
    plane_ids = list(dict.fromkeys(new_ids + on_plane))
    if len(plane_ids) != 2:
        raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
    return cleaned, (plane_ids[0], plane_ids[1])


def _order_loop(edges: set[tuple[int, int]]) -> np.ndarray:
    graph: dict[int, set[int]] = defaultdict(set)
    for first, second in edges:
        graph[first].add(second)
        graph[second].add(first)
    if len(graph) < 3 or any(len(neighbors) != 2 for neighbors in graph.values()):
        raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
    start = min(graph)
    ordered = [start]
    previous = -1
    current = start
    while True:
        options = sorted(graph[current] - ({previous} if previous >= 0 else set()))
        if not options:
            raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
        following = options[0]
        if following == start:
            break
        if following in ordered:
            raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
        ordered.append(following)
        previous, current = current, following
        if len(ordered) > len(graph):
            raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
    if len(ordered) != len(graph):
        raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
    return np.asarray(ordered, dtype=np.int64)


def local_plane_cut(
    surface: TaggedSurface,
    boundary: BoundaryInput,
    *,
    radial_factor: float,
    axial_back_factor: float,
    axial_forward_factor: float,
) -> tuple[TaggedSurface, CutLoop, dict[str, Any]]:
    """Remove only the target rounded end and exactly clip crossing triangles."""

    vertices_array = np.asarray(surface.vertices, dtype=float)
    normal = boundary.outward_normal
    relative = vertices_array - boundary.center_um
    axial = relative @ normal
    face_points = vertices_array[surface.faces]
    centroids = np.mean(face_points, axis=1)
    face_relative = centroids - boundary.center_um
    face_axial = face_relative @ normal
    face_radial = np.linalg.norm(
        face_relative - np.outer(face_axial, normal), axis=1
    )
    radial_limit = radial_factor * boundary.source_radius_um
    back_limit = axial_back_factor * boundary.source_radius_um
    forward_limit = axial_forward_factor * boundary.source_radius_um
    candidate_ids = np.flatnonzero(
        (face_radial <= radial_limit)
        & (face_axial >= -back_limit)
        & (face_axial <= forward_limit)
    )
    components = _candidate_components(surface.faces, candidate_ids)
    if len(components) != 1:
        raise SurfacePrepareError("LOCAL_PORT_CUT_AMBIGUOUS")
    target_component = components[0]
    face_distances = axial[surface.faces]
    crossing_ids = {
        face_id
        for face_id in target_component
        if float(np.min(face_distances[face_id])) < -1.0e-12
        and float(np.max(face_distances[face_id])) > 1.0e-12
    }
    if not crossing_ids:
        raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
    outward_ids = {
        face_id
        for face_id in target_component
        if face_id not in crossing_ids
        and float(np.min(face_distances[face_id])) >= -1.0e-12
    }
    vertices = [point.copy() for point in vertices_array]
    intersection_cache: dict[tuple[int, int], int] = {}
    loop_edges: set[tuple[int, int]] = set()
    new_faces: list[tuple[int, int, int]] = []
    new_type: list[int] = []
    new_index: list[int] = []
    new_origin: list[int] = []
    new_kind: list[int] = []
    new_extension_index: list[int] = []
    new_extension_band: list[int] = []
    for face_id, face in enumerate(surface.faces):
        if face_id in outward_ids:
            continue
        if face_id not in crossing_ids:
            new_faces.append(tuple(map(int, face)))
            new_type.append(int(surface.boundary_type[face_id]))
            new_index.append(int(surface.boundary_index[face_id]))
            new_origin.append(int(surface.boundary_origin[face_id]))
            new_kind.append(int(surface.face_kind[face_id]))
            new_extension_index.append(int(surface.extension_index[face_id]))
            new_extension_band.append(int(surface.extension_band[face_id]))
            continue
        polygon, segment = _clip_face(
            face, axial, vertices, intersection_cache
        )
        loop_edges.add(tuple(sorted(segment)))
        for index in range(1, len(polygon) - 1):
            new_faces.append((polygon[0], polygon[index], polygon[index + 1]))
            new_type.append(0)
            new_index.append(-1)
            new_origin.append(0)
            new_kind.append(0)
            new_extension_index.append(-1)
            new_extension_band.append(-1)
    loop_ids = _order_loop(loop_edges)
    loop_points = np.asarray([vertices[index] for index in loop_ids], dtype=float)
    area, loop_center, projected = polygon_metrics(loop_points, normal)
    if not LineString(np.vstack((projected, projected[0]))).is_simple:
        raise SurfacePrepareError("PORT_CUT_LOOP_INVALID")
    area_vector = 0.5 * np.sum(
        np.cross(loop_points, np.roll(loop_points, -1, axis=0)), axis=0
    )
    if float(np.dot(area_vector, normal)) < 0:
        loop_ids = loop_ids[::-1].copy()
        loop_points = loop_points[::-1].copy()
    residual = float(np.max(np.abs((loop_points - boundary.center_um) @ normal)))
    equivalent_radius = float(np.sqrt(area / np.pi))
    output = TaggedSurface(
        vertices=np.asarray(vertices, dtype=float),
        faces=np.asarray(new_faces, dtype=np.int64),
        boundary_type=np.asarray(new_type, dtype=np.uint8),
        boundary_index=np.asarray(new_index, dtype=np.int32),
        boundary_origin=np.asarray(new_origin, dtype=np.uint8),
        face_kind=np.asarray(new_kind, dtype=np.uint8),
        extension_index=np.asarray(new_extension_index, dtype=np.int32),
        extension_band=np.asarray(new_extension_band, dtype=np.int32),
        source_vertex_index=np.concatenate(
            (
                surface.source_vertex_index,
                np.full(len(vertices) - len(surface.vertices), -1, dtype=np.int64),
            )
        ),
    )
    qc = {
        "port_id": boundary.port_id,
        "candidate_surface_component_count": len(components),
        "crossing_triangle_count": len(crossing_ids),
        "removed_outward_triangle_count": len(outward_ids),
        "cut_loop_count": 1,
        "cut_loop_vertex_count": len(loop_ids),
        "cut_loop_area_um2": area,
        "equivalent_radius_um": equivalent_radius,
        "equivalent_radius_relative_error": abs(
            equivalent_radius - boundary.source_radius_um
        )
        / boundary.source_radius_um,
        "loop_center_um": loop_center.tolist(),
        "expected_center_um": boundary.center_um.tolist(),
        "loop_center_distance_um": float(np.linalg.norm(loop_center - boundary.center_um)),
        "plane_residual_um": residual,
        "local_radial_limit_um": radial_limit,
        "local_axial_back_limit_um": back_limit,
        "local_axial_forward_limit_um": forward_limit,
    }
    return output, CutLoop(
        boundary_index=boundary.index,
        vertex_ids=loop_ids,
        center_um=loop_center,
        outward_normal=normal.copy(),
        area_um2=area,
        equivalent_radius_um=equivalent_radius,
        plane_residual_um=residual,
    ), qc
