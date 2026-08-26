"""Straight extrusion of actual cut loops and robust planar cap triangulation."""

from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import triangulate

from .config import ExtensionMeshConfig
from .io import BoundaryInput, SurfacePrepareError
from .local_cut import orthogonal_basis, polygon_metrics
from .mesh_refinement import build_refined_rings, choose_quad_diagonal
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
    *,
    local_original_median_edge_length_um: float,
    mesh_config: ExtensionMeshConfig,
) -> tuple[TaggedSurface, BoundarySurfaceResult]:
    """Create a structured multi-ring extension without moving Ring 0 or the core."""

    proximal_ids = np.asarray(loop.vertex_ids, dtype=np.int64)
    proximal = surface.vertices[proximal_ids]
    refined = build_refined_rings(
        proximal,
        boundary.outward_normal,
        boundary.extension_length_um,
        local_original_median_edge_length_um,
        boundary.source_radius_um,
        mesh_config,
    )
    new_points = refined.points[1:].reshape((-1, 3))
    vertices = np.vstack((surface.vertices, new_points))
    ring_ids = [proximal_ids]
    next_id = len(surface.vertices)
    for _ in range(1, refined.ring_count):
        ids = np.arange(next_id, next_id + len(proximal_ids), dtype=np.int64)
        ring_ids.append(ids)
        next_id += len(proximal_ids)
    distal_ids = ring_ids[-1]
    distal = vertices[distal_ids]
    side_faces: list[tuple[int, int, int]] = []
    side_bands: list[int] = []
    for band, (first_ring, second_ring) in enumerate(
        zip(ring_ids[:-1], ring_ids[1:], strict=True)
    ):
        for index in range(len(proximal_ids)):
            faces = choose_quad_diagonal(
                vertices, first_ring, second_ring, index
            )
            side_faces.extend(faces)
            side_bands.extend((band, band))
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
        extension_index=np.concatenate(
            (
                surface.extension_index,
                np.full(side_count, boundary.index, dtype=np.int32),
                np.full(cap_count, boundary.index, dtype=np.int32),
            )
        ),
        extension_band=np.concatenate(
            (
                surface.extension_band,
                np.asarray(side_bands, dtype=np.int32),
                np.full(cap_count, -1, dtype=np.int32),
            )
        ),
        source_vertex_index=np.concatenate(
            (
                surface.source_vertex_index,
                np.full(len(new_points), -1, dtype=np.int64),
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
    actual_cap_area, distal_center, _ = polygon_metrics(
        distal, boundary.outward_normal
    )
    _, proximal_center, proximal_projected = polygon_metrics(
        proximal, boundary.outward_normal
    )
    first, second = orthogonal_basis(boundary.outward_normal)
    intermediate_centers: list[np.ndarray] = []
    intermediate_areas: list[float] = []
    intermediate_validity: list[bool] = []
    intermediate_signed_areas: list[float] = []
    intermediate_axial_errors: list[float] = []
    for ring, station in zip(
        refined.points[1:-1], refined.stations_um[1:-1], strict=True
    ):
        origin = np.mean(ring, axis=0)
        projected = np.column_stack(
            ((ring - origin) @ first, (ring - origin) @ second)
        )
        polygon = Polygon(projected)
        center_2d = np.asarray(polygon.centroid.coords[0], dtype=float)
        intermediate_centers.append(
            origin + center_2d[0] * first + center_2d[1] * second
        )
        intermediate_areas.append(float(polygon.area))
        intermediate_validity.append(
            bool(
                polygon.is_valid
                and polygon.exterior.is_simple
                and np.isfinite(polygon.area)
                and polygon.area > 0.0
            )
        )
        following = np.roll(projected, -1, axis=0)
        intermediate_signed_areas.append(
            0.5
            * float(
                np.sum(
                    projected[:, 0] * following[:, 1]
                    - following[:, 0] * projected[:, 1]
                )
            )
        )
        intermediate_axial_errors.append(
            float(
                np.max(
                    np.abs(
                        (ring - proximal_center) @ boundary.outward_normal
                        - float(station)
                    )
                )
            )
        )
    intermediate_centers_array = np.asarray(intermediate_centers)
    centerline_displacements = intermediate_centers_array - proximal_center
    centerline_radial_offsets = centerline_displacements - np.outer(
        centerline_displacements @ boundary.outward_normal,
        boundary.outward_normal,
    )
    intermediate_drifts = np.linalg.norm(centerline_radial_offsets, axis=1)
    worst_intermediate_offset = int(np.argmax(intermediate_drifts))
    intermediate_areas_array = np.asarray(intermediate_areas)
    target_intermediate_areas = refined.target_areas_um2[1:-1]
    area_relative_errors = np.abs(
        intermediate_areas_array - target_intermediate_areas
    ) / target_intermediate_areas
    proximal_following = np.roll(proximal_projected, -1, axis=0)
    proximal_signed_area = 0.5 * float(
        np.sum(
            proximal_projected[:, 0] * proximal_following[:, 1]
            - proximal_following[:, 0] * proximal_projected[:, 1]
        )
    )
    measured_lengths = (distal - proximal) @ boundary.outward_normal
    center_displacement = distal_center - proximal_center
    center_displacement_norm = float(np.linalg.norm(center_displacement))
    if center_displacement_norm <= np.finfo(float).eps:
        raise SurfacePrepareError("PORT_EXTENSION_AXIS_INVALID")
    extension_axis_dot = float(
        np.dot(center_displacement, boundary.outward_normal)
        / center_displacement_norm
    )
    return output, BoundarySurfaceResult(
        boundary_index=boundary.index,
        port_id=boundary.port_id,
        boundary_origin=boundary.boundary_origin,
        role=boundary.role,
        source_radius_um=boundary.source_radius_um,
        extension_length_um=boundary.extension_length_um,
        actual_cap_area_um2=actual_cap_area,
        equivalent_radius_um=float(np.sqrt(actual_cap_area / np.pi)),
        cap_planarity_error_um=float(np.max(cap_residual)),
        minimum_cap_normal_dot=float(np.min(normal_dots)),
        extension_length_error_um=float(
            np.max(np.abs(measured_lengths - boundary.extension_length_um))
        ),
        extension_axis_dot=extension_axis_dot,
        intermediate_ring_centerline_max_deviation_um=(
            float(np.max(intermediate_drifts))
        ),
        intermediate_ring_centerline_p95_deviation_um=float(
            np.percentile(intermediate_drifts, 95)
        ),
        intermediate_ring_centerline_mean_deviation_um=float(
            np.mean(intermediate_drifts)
        ),
        intermediate_ring_centerline_worst_ring_index=(
            worst_intermediate_offset + 1
        ),
        intermediate_ring_centerline_worst_ring_axial_station_um=float(
            refined.stations_um[worst_intermediate_offset + 1]
        ),
        intermediate_ring_axial_station_max_error_um=float(
            np.max(intermediate_axial_errors)
        ),
        intermediate_ring_area_relative_error_max=float(
            np.max(area_relative_errors)
        ),
        intermediate_ring_all_areas_finite_positive=bool(
            np.all(np.isfinite(intermediate_areas_array))
            and np.all(intermediate_areas_array > 0.0)
        ),
        intermediate_ring_all_polygons_simple_valid=bool(
            all(intermediate_validity)
        ),
        intermediate_ring_all_orientations_consistent=bool(
            np.all(
                np.sign(np.asarray(intermediate_signed_areas))
                == np.sign(proximal_signed_area)
            )
        ),
        local_original_median_edge_length_um=local_original_median_edge_length_um,
        target_edge_length_um=refined.target_edge_length_um,
        ring_count=refined.ring_count,
        transition_length_um=refined.transition_length_um,
        proximal_ring_max_motion_um=float(
            np.max(np.linalg.norm(vertices[proximal_ids] - proximal, axis=1))
        ),
        distal_ring_max_motion_um=float(
            np.max(
                np.linalg.norm(
                    distal
                    - (
                        refined.regularized_loop
                        + boundary.extension_length_um * boundary.outward_normal
                    ),
                    axis=1,
                )
            )
        ),
        cut_loop_vertex_count=len(proximal_ids),
        cap_face_indices=cap_indices,
        side_face_indices=side_indices,
    )
