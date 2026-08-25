"""Straight extrusion of actual cut loops and robust planar cap triangulation."""

from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import triangulate

from .io import BoundaryInput, SurfacePrepareError
from .local_cut import orthogonal_basis
from .types import BoundarySurfaceResult, CutLoop, TaggedSurface


def triangulate_cap(points: np.ndarray, normal: np.ndarray) -> list[tuple[int, int, int]]:
    """Triangulate a simple, possibly non-convex planar polygon with GEOS Delaunay filtering."""

    first, second = orthogonal_basis(normal)
    origin = np.mean(points, axis=0)
    points_2d = np.column_stack(((points - origin) @ first, (points - origin) @ second))
    polygon = Polygon(points_2d)
    if not polygon.is_valid or not polygon.exterior.is_simple or polygon.area <= 0:
        raise SurfacePrepareError("PORT_CAP_TRIANGULATION_FAILED")
    triangles = [
        triangle
        for triangle in triangulate(polygon)
        if polygon.covers(triangle.representative_point())
        and polygon.intersection(triangle).area >= triangle.area * (1.0 - 1.0e-10)
    ]
    result: list[tuple[int, int, int]] = []
    accumulated_area = 0.0
    for triangle in triangles:
        coordinates = np.asarray(triangle.exterior.coords[:-1], dtype=float)
        indices: list[int] = []
        for coordinate in coordinates:
            distances = np.linalg.norm(points_2d - coordinate, axis=1)
            index = int(np.argmin(distances))
            if distances[index] > 1.0e-8:
                raise SurfacePrepareError("PORT_CAP_TRIANGULATION_FAILED")
            indices.append(index)
        if len(set(indices)) != 3:
            raise SurfacePrepareError("PORT_CAP_TRIANGULATION_FAILED")
        face = tuple(indices)
        face_normal = np.cross(points[face[1]] - points[face[0]], points[face[2]] - points[face[0]])
        if float(np.dot(face_normal, normal)) < 0:
            face = (face[0], face[2], face[1])
        result.append(face)
        accumulated_area += float(triangle.area)
    if not result or abs(accumulated_area - polygon.area) > max(1.0e-10, polygon.area * 1.0e-8):
        raise SurfacePrepareError("PORT_CAP_TRIANGULATION_FAILED")
    return result


def extrude_and_cap(
    surface: TaggedSurface,
    loop: CutLoop,
    boundary: BoundaryInput,
) -> tuple[TaggedSurface, BoundarySurfaceResult]:
    """Rigidly translate the actual cut loop, add a wall, then add one tagged flat cap."""

    proximal_ids = np.asarray(loop.vertex_ids, dtype=np.int64)
    proximal = surface.vertices[proximal_ids]
    distal = proximal + boundary.outward_normal * boundary.extension_length_um
    first_distal_id = len(surface.vertices)
    distal_ids = np.arange(first_distal_id, first_distal_id + len(distal), dtype=np.int64)
    vertices = np.vstack((surface.vertices, distal))
    side_faces: list[tuple[int, int, int]] = []
    count = len(proximal_ids)
    for index in range(count):
        following = (index + 1) % count
        first = int(proximal_ids[index])
        next_first = int(proximal_ids[following])
        second = int(distal_ids[index])
        next_second = int(distal_ids[following])
        side_faces.extend(((first, next_first, next_second), (first, next_second, second)))
    local_cap_faces = triangulate_cap(distal, boundary.outward_normal)
    cap_faces = [
        tuple(int(distal_ids[index]) for index in face) for face in local_cap_faces
    ]
    side_start = len(surface.faces)
    cap_start = side_start + len(side_faces)
    faces = np.vstack(
        (
            surface.faces,
            np.asarray(side_faces, dtype=np.int64),
            np.asarray(cap_faces, dtype=np.int64),
        )
    )
    side_count = len(side_faces)
    cap_count = len(cap_faces)
    boundary_type_code = 1 if boundary.role == "ASSUMED_INLET" else 2
    origin_code = 1 if boundary.boundary_origin == "CUT_PORT" else 2
    output = TaggedSurface(
        vertices=vertices,
        faces=faces,
        boundary_type=np.concatenate(
            (
                surface.boundary_type,
                np.zeros(side_count, dtype=np.uint8),
                np.full(cap_count, boundary_type_code, dtype=np.uint8),
            )
        ),
        boundary_index=np.concatenate(
            (
                surface.boundary_index,
                np.full(side_count, -1, dtype=np.int32),
                np.full(cap_count, boundary.index, dtype=np.int32),
            )
        ),
        boundary_origin=np.concatenate(
            (
                surface.boundary_origin,
                np.zeros(side_count, dtype=np.uint8),
                np.full(cap_count, origin_code, dtype=np.uint8),
            )
        ),
        face_kind=np.concatenate(
            (
                surface.face_kind,
                np.ones(side_count, dtype=np.uint8),
                np.full(cap_count, 2, dtype=np.uint8),
            )
        ),
    )
    cap_indices = np.arange(cap_start, cap_start + cap_count, dtype=np.int64)
    side_indices = np.arange(side_start, side_start + side_count, dtype=np.int64)
    cap_triangles = vertices[faces[cap_indices]]
    normals = np.cross(
        cap_triangles[:, 1] - cap_triangles[:, 0],
        cap_triangles[:, 2] - cap_triangles[:, 0],
    )
    normal_lengths = np.linalg.norm(normals, axis=1)
    if np.any(normal_lengths <= np.finfo(float).eps):
        raise SurfacePrepareError("PORT_CAP_TRIANGULATION_FAILED")
    normals /= normal_lengths[:, None]
    normal_dots = normals @ boundary.outward_normal
    cap_residual = np.abs((distal - boundary.extension_end_um) @ boundary.outward_normal)
    measured_lengths = (distal - proximal) @ boundary.outward_normal
    displacement = distal - proximal
    displacement_norms = np.linalg.norm(displacement, axis=1)
    axis_dots = (
        displacement @ boundary.outward_normal / displacement_norms
    )
    return output, BoundarySurfaceResult(
        boundary_index=boundary.index,
        port_id=boundary.port_id,
        boundary_origin=boundary.boundary_origin,
        role=boundary.role,
        source_radius_um=boundary.source_radius_um,
        extension_length_um=boundary.extension_length_um,
        actual_cap_area_um2=loop.area_um2,
        equivalent_radius_um=loop.equivalent_radius_um,
        cap_planarity_error_um=float(np.max(cap_residual)),
        minimum_cap_normal_dot=float(np.min(normal_dots)),
        extension_length_error_um=float(
            np.max(np.abs(measured_lengths - boundary.extension_length_um))
        ),
        extension_axis_dot=float(np.min(axis_dots)),
        cut_loop_vertex_count=len(proximal_ids),
        cap_face_indices=cap_indices,
        side_face_indices=side_indices,
    )
