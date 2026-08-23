"""v5 continuous port centerlines and collar boundary-loop stitching."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
import trimesh
from shapely import LinearRing

from .config import CFDLumenConfig
from .local_implicit_junction import sample_from_junction
from .lumen_builder import build_variable_radius_tube
from .surface_qc import _orthogonal_basis
from .types import BranchGeometry, GeometryValidationError, JunctionCollar, PortGeometry


SOURCE_CENTERLINE = np.uint8(0)
CORE_BOUNDARY_EXTENSION_ORIGIN = np.uint8(1)
CFD_DERIVED_EXTENSION = np.uint8(2)

JUNCTION_CORE_FACE = np.uint8(1)
TRANSITION_COLLAR_FACE = np.uint8(2)
EXPLICIT_BRANCH_FACE = np.uint8(3)


def extend_branches_to_cfd_ports(
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> tuple[list[BranchGeometry], list[dict[str, Any]]]:
    """Append/prepend each straight constant-radius port before one tube generation."""

    port_by_node = {int(port.local_node_id): port for port in ports}
    output: list[BranchGeometry] = []
    rows: list[dict[str, Any]] = []
    for branch in branches:
        extended = replace(branch)
        points = np.asarray(branch.points_um, dtype=float).copy()
        radii = np.asarray(branch.radius_um, dtype=float).copy()
        point_type = np.full(len(points), SOURCE_CENTERLINE, dtype=np.uint8)
        cut_ids = [""] * len(points)
        edge_ids = np.full(len(points), -1, dtype=np.int64)
        core_distance = np.full(len(points), np.nan, dtype=float)
        core_positions = np.full((len(points), 3), np.nan, dtype=float)
        for endpoint in (-1, 0):
            local_node_id = int(branch.local_node_ids[endpoint])
            port = port_by_node.get(local_node_id)
            if port is None:
                continue
            spacing = float(
                np.clip(
                    config.geometry.resample_radius_fraction * port.radius_um,
                    config.geometry.min_resample_spacing_um,
                    config.geometry.max_resample_spacing_um,
                )
            )
            count = max(1, int(np.ceil(port.extension_length_um / spacing)))
            distances = np.linspace(0.0, port.extension_length_um, count + 1)
            extension_points = (
                port.original_position_um[None, :]
                + distances[:, None] * port.outward_tangent[None, :]
            )
            original_core = np.asarray(
                port.original_core_cut_position_um
                if port.original_core_cut_position_um is not None
                else port.original_position_um,
                dtype=float,
            )
            source_cut_id = port.source_core_cut_port_id or port.cut_port_id
            for index, (distance, point) in enumerate(zip(distances, extension_points)):
                rows.append(
                    {
                        "branch_id": branch.branch_id,
                        "extension_point_index": index,
                        "point_type": (
                            "CFD_EXTENSION_ORIGIN" if index == 0 else "CFD_EXTENSION"
                        ),
                        "source_cut_port_id": source_cut_id,
                        "active_cfd_port_id": port.cut_port_id,
                        "source_global_edge_id": port.global_edge_id,
                        "distance_from_cfd_extension_origin_um": float(distance),
                        "distance_from_core_boundary_um": float(
                            port.source_core_to_cfd_cut_length_um + distance
                        ),
                        "x_um": float(point[0]),
                        "y_um": float(point[1]),
                        "z_um": float(point[2]),
                        "radius_um": port.radius_um,
                        "original_core_cut_x_um": float(original_core[0]),
                        "original_core_cut_y_um": float(original_core[1]),
                        "original_core_cut_z_um": float(original_core[2]),
                        "original_source_boundary_type": "CORE_ROI_BOUNDARY",
                        "extension_origin_type": "CFD_DOMAIN_SOURCE_CUT",
                        "boundary_type_at_end": "CFD_BOUNDARY_PORT",
                        "point_boundary_type": (
                            "CFD_BOUNDARY_PORT"
                            if index == len(distances) - 1
                            else "CFD_DOMAIN_SOURCE_CUT"
                            if index == 0
                            else "NONE"
                        ),
                    }
                )
            derived_type = np.full(count + 1, CFD_DERIVED_EXTENSION, dtype=np.uint8)
            derived_type[0] = CORE_BOUNDARY_EXTENSION_ORIGIN
            derived_ids = [source_cut_id] * (count + 1)
            derived_edges = np.full(count + 1, port.global_edge_id, dtype=np.int64)
            derived_distances = port.source_core_to_cfd_cut_length_um + distances
            derived_core = np.repeat(original_core[None, :], count + 1, axis=0)
            if endpoint == -1:
                points = np.vstack((points[:-1], extension_points))
                radii = np.concatenate(
                    (radii[:-1], np.full(count + 1, port.radius_um, dtype=float))
                )
                point_type = np.concatenate((point_type[:-1], derived_type))
                cut_ids = [*cut_ids[:-1], *derived_ids]
                edge_ids = np.concatenate((edge_ids[:-1], derived_edges))
                core_distance = np.concatenate((core_distance[:-1], derived_distances))
                core_positions = np.vstack((core_positions[:-1], derived_core))
            else:
                points = np.vstack((extension_points[::-1], points[1:]))
                radii = np.concatenate(
                    (np.full(count + 1, port.radius_um, dtype=float), radii[1:])
                )
                point_type = np.concatenate((derived_type[::-1], point_type[1:]))
                cut_ids = [*derived_ids[::-1], *cut_ids[1:]]
                edge_ids = np.concatenate((derived_edges[::-1], edge_ids[1:]))
                core_distance = np.concatenate((derived_distances[::-1], core_distance[1:]))
                core_positions = np.vstack((derived_core[::-1], core_positions[1:]))
        extended.points_um = points
        extended.radius_um = radii
        extended.arc_length_um = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
        )
        extended.construction_point_type = point_type
        extended.construction_source_cut_port_id = tuple(cut_ids)
        extended.construction_source_global_edge_id = edge_ids
        extended.construction_distance_from_core_boundary_um = core_distance
        extended.construction_original_core_cut_position_um = core_positions
        output.append(extended)
    return output, rows


def _boundary_loops(mesh: trimesh.Trimesh) -> list[np.ndarray]:
    edges = np.sort(np.asarray(mesh.edges, dtype=np.int64), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    boundary = unique[counts == 1]
    if not len(boundary):
        return []
    adjacency: dict[int, list[int]] = {}
    for first, second in boundary:
        adjacency.setdefault(int(first), []).append(int(second))
        adjacency.setdefault(int(second), []).append(int(first))
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise GeometryValidationError(
            "COLLAR_LOOP_EXTRACTION_FAILED: boundary loop contains an open/branched vertex"
        )
    remaining = {tuple(map(int, edge)) for edge in boundary}
    loops: list[np.ndarray] = []
    while remaining:
        first_edge = min(remaining)
        start, current = first_edge
        ordered = [start]
        previous = start
        remaining.discard(tuple(sorted((start, current))))
        for _ in range(len(boundary) + 1):
            ordered.append(current)
            candidates = [node for node in adjacency[current] if node != previous]
            if len(candidates) != 1:
                raise GeometryValidationError(
                    "COLLAR_LOOP_EXTRACTION_FAILED: boundary ordering is ambiguous"
                )
            following = candidates[0]
            if following == start:
                remaining.discard(tuple(sorted((current, following))))
                break
            remaining.discard(tuple(sorted((current, following))))
            previous, current = current, following
        else:
            raise GeometryValidationError(
                "COLLAR_LOOP_EXTRACTION_FAILED: boundary loop did not close"
            )
        if len(ordered) < 3:
            raise GeometryValidationError(
                "COLLAR_LOOP_EXTRACTION_FAILED: boundary loop has fewer than three vertices"
            )
        loops.append(np.asarray(ordered, dtype=np.int64))
    return loops


def _validate_loop_geometry(
    mesh: trimesh.Trimesh,
    loop: np.ndarray,
    center: np.ndarray,
    normal: np.ndarray,
) -> None:
    first, second = _orthogonal_basis(normal)
    relative = np.asarray(mesh.vertices[loop], dtype=float) - center
    projected = np.column_stack((relative @ first, relative @ second))
    ring = LinearRing(projected)
    if not ring.is_closed or not ring.is_simple:
        raise GeometryValidationError(
            "COLLAR_LOOP_EXTRACTION_FAILED: loop is open or self-intersecting"
        )


def _loop_at_plane(
    mesh: trimesh.Trimesh,
    loops: list[np.ndarray],
    center: np.ndarray,
    normal: np.ndarray,
    tolerance_um: float,
) -> np.ndarray:
    candidates = [
        loop
        for loop in loops
        if float(
            np.max(
                np.abs((np.asarray(mesh.vertices[loop]) - center[None, :]) @ normal)
            )
        )
        <= tolerance_um
    ]
    if len(candidates) != 1:
        raise GeometryValidationError(
            "COLLAR_LOOP_EXTRACTION_FAILED: expected exactly one closed loop at collar "
            f"plane, found {len(candidates)}"
        )
    _validate_loop_geometry(mesh, candidates[0], center, normal)
    return candidates[0]


def _collapse_boundary_loops_to_count(
    mesh: trimesh.Trimesh,
    loops: list[np.ndarray],
    target_count: int,
) -> trimesh.Trimesh:
    """Locally collapse only clipped boundary edges to uniform transition rings."""

    vertices = np.asarray(mesh.vertices, dtype=float).copy()
    remap = np.arange(len(vertices), dtype=np.int64)
    for loop in loops:
        if len(loop) < target_count:
            raise GeometryValidationError(
                "COLLAR_LOOP_EXTRACTION_FAILED: implicit loop has fewer vertices than N_loop"
            )
        points = vertices[loop]
        closed = np.vstack((points, points[0]))
        edge_lengths = np.linalg.norm(np.diff(closed, axis=0), axis=1)
        perimeter = float(edge_lengths.sum())
        vertex_arc = np.concatenate(([0.0], np.cumsum(edge_lengths[:-1])))
        targets = np.arange(target_count, dtype=float) * perimeter / target_count
        representatives: list[int] = []
        for target in targets:
            circular = np.minimum(
                np.abs(vertex_arc - target), perimeter - np.abs(vertex_arc - target)
            )
            candidate = int(loop[int(np.argmin(circular))])
            if candidate in representatives:
                raise GeometryValidationError(
                    "COLLAR_LOOP_EXTRACTION_FAILED: loop collapse representatives are not unique"
                )
            representatives.append(candidate)
        # Keep each surviving boundary vertex at its original marching-cubes
        # position. Moving a representative onto an interpolated perimeter
        # sample can fold its incident core triangles across a nearby face.
        # The following transition rings perform the geometric resampling;
        # this operation is only the local, order-preserving boundary remesh.
        group = np.floor(vertex_arc / perimeter * target_count + 0.5).astype(int)
        group %= target_count
        for vertex_id, group_id in zip(loop, group):
            remap[int(vertex_id)] = representatives[int(group_id)]
    faces = remap[np.asarray(mesh.faces, dtype=np.int64)]
    keep = (
        (faces[:, 0] != faces[:, 1])
        & (faces[:, 1] != faces[:, 2])
        & (faces[:, 2] != faces[:, 0])
    )
    collapsed = trimesh.Trimesh(vertices=vertices, faces=faces[keep], process=False)
    collapsed.update_faces(collapsed.unique_faces())
    collapsed.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(collapsed, multibody=True)
    return collapsed


def clip_implicit_patch_to_transition_inner(
    mesh: trimesh.Trimesh,
    collars: list[JunctionCollar],
    branches_by_id: dict[int, BranchGeometry],
    tolerance_um: float,
    loop_vertex_count: int,
) -> tuple[trimesh.Trimesh, dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    clipped = mesh.copy()
    planes: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for collar in collars:
        branch = branches_by_id[collar.branch_id]
        center, _, tangent = sample_from_junction(
            branch, collar.endpoint_index, collar.explicit_cap_distance_um
        )
        candidate = trimesh.intersections.slice_mesh_plane(
            clipped,
            plane_normal=-tangent,
            plane_origin=center,
            cap=False,
        )
        if candidate is None or not len(candidate.faces):
            raise GeometryValidationError(
                f"COLLAR_LOOP_EXTRACTION_FAILED: implicit clip empty for branch {collar.branch_id}"
            )
        clipped = candidate
        clipped.merge_vertices(digits_vertex=8)
        clipped.remove_unreferenced_vertices()
        planes[collar.branch_id] = (center, tangent, np.asarray(collar.collar_position_um))
    loops = _boundary_loops(clipped)
    for collar in collars:
        center, tangent, _ = planes[collar.branch_id]
        _loop_at_plane(clipped, loops, center, tangent, tolerance_um)
    # Keep the exact clipped implicit boundary. It is resampled to N_loop by
    # the first transition-ring zipper in ``build_transition_strip``. Directly
    # collapsing marching-cubes boundary vertices can fold adjacent core faces,
    # violating v5's zero-self-intersection requirement.
    output: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for collar in collars:
        center, tangent, _ = planes[collar.branch_id]
        loop = _loop_at_plane(clipped, loops, center, tangent, tolerance_um)
        output[collar.branch_id] = (loop, center, tangent)
    return clipped, output


def trim_branches_to_transition_outer(
    branches: list[BranchGeometry],
    collars_by_junction: dict[int, list[JunctionCollar]],
) -> tuple[list[BranchGeometry], dict[int, list[JunctionCollar]]]:
    by_branch: dict[int, list[JunctionCollar]] = {}
    for collars in collars_by_junction.values():
        for collar in collars:
            by_branch.setdefault(collar.branch_id, []).append(collar)
    output: list[BranchGeometry] = []
    for branch in branches:
        lower = 0.0
        upper = float(branch.arc_length_um[-1])
        for collar in by_branch.get(branch.branch_id, []):
            if collar.endpoint_index == 0:
                lower = max(lower, collar.implicit_extent_um)
            else:
                upper = min(
                    upper, float(branch.arc_length_um[-1]) - collar.implicit_extent_um
                )
        if upper - lower <= 1.0e-8:
            raise GeometryValidationError(
                f"COLLAR_LOOP_EXTRACTION_FAILED: branch {branch.branch_id} has no explicit region"
            )
        inside = branch.arc_length_um[
            (branch.arc_length_um > lower) & (branch.arc_length_um < upper)
        ]
        targets = np.concatenate(([lower], inside, [upper]))
        points = np.column_stack(
            [
                np.interp(targets, branch.arc_length_um, branch.points_um[:, axis])
                for axis in range(3)
            ]
        )
        radii = np.interp(targets, branch.arc_length_um, branch.radius_um)
        trimmed = replace(branch)
        trimmed.points_um = points
        trimmed.radius_um = radii
        trimmed.arc_length_um = targets - lower
        output.append(trimmed)
    return output, by_branch


def _remove_collar_caps(
    mesh: trimesh.Trimesh,
    branch: BranchGeometry,
    collars: list[JunctionCollar],
) -> tuple[trimesh.Trimesh, dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    keep = np.ones(len(mesh.faces), dtype=bool)
    planes: dict[int, tuple[np.ndarray, np.ndarray, float]] = {}
    for collar in collars:
        endpoint = 0 if collar.endpoint_index == 0 else -1
        center = np.asarray(branch.points_um[endpoint], dtype=float)
        tangent = (
            np.asarray(branch.points_um[1] - branch.points_um[0], dtype=float)
            if endpoint == 0
            else np.asarray(branch.points_um[-2] - branch.points_um[-1], dtype=float)
        )
        tangent /= np.linalg.norm(tangent)
        radius = float(branch.radius_um[endpoint])
        relative = np.asarray(mesh.triangles_center) - center[None, :]
        axial = np.abs(relative @ tangent)
        radial = np.linalg.norm(
            relative - (relative @ tangent)[:, None] * tangent[None, :], axis=1
        )
        tolerance = max(radius * 1.0e-4, 1.0e-8)
        cap_faces = (
            (axial <= tolerance)
            & (radial <= 1.01 * radius)
            & (np.abs(np.asarray(mesh.face_normals) @ tangent) >= 0.95)
        )
        if not np.any(cap_faces):
            raise GeometryValidationError(
                f"COLLAR_LOOP_EXTRACTION_FAILED: no explicit cap at branch {branch.branch_id}"
            )
        keep &= ~cap_faces
        planes[collar.junction_node_id] = (center, tangent, radius)
    opened = mesh.copy()
    opened.update_faces(keep)
    opened.remove_unreferenced_vertices()
    loops = _boundary_loops(opened)
    if len(loops) != len(planes):
        raise GeometryValidationError(
            "COLLAR_LOOP_EXTRACTION_FAILED: explicit tube boundary loop count "
            f"{len(loops)} != collar count {len(planes)}"
        )
    output: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    available = set(range(len(loops)))
    for junction_id, (center, tangent, radius) in sorted(planes.items()):
        loop_index = min(
            available,
            key=lambda index: float(
                np.linalg.norm(np.asarray(opened.vertices[loops[index]]).mean(axis=0) - center)
            ),
        )
        loop = loops[loop_index]
        center_error = float(
            np.linalg.norm(np.asarray(opened.vertices[loop]).mean(axis=0) - center)
        )
        if center_error > 0.25 * radius:
            raise GeometryValidationError(
                "COLLAR_LOOP_EXTRACTION_FAILED: explicit loop center differs from collar "
                f"by {center_error:.6g} um"
            )
        available.remove(loop_index)
        _validate_loop_geometry(opened, loop, center, tangent)
        output[junction_id] = (loop, center, tangent)
    return opened, output


def build_open_explicit_tubes(
    branches: list[BranchGeometry],
    collars_by_branch: dict[int, list[JunctionCollar]],
    tube_sides: int,
) -> tuple[
    dict[int, trimesh.Trimesh],
    dict[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]],
]:
    meshes: dict[int, trimesh.Trimesh] = {}
    loops: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for branch in branches:
        mesh = build_variable_radius_tube(branch, tube_sides)
        collars = collars_by_branch.get(branch.branch_id, [])
        if collars:
            mesh, branch_loops = _remove_collar_caps(mesh, branch, collars)
            for junction_id, payload in branch_loops.items():
                loops[(junction_id, branch.branch_id)] = payload
        meshes[branch.branch_id] = mesh
    return meshes, loops


def _resample_loop(points: np.ndarray, count: int) -> np.ndarray:
    closed = np.vstack((points, points[0]))
    lengths = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    perimeter = float(lengths.sum())
    if perimeter <= 0:
        raise GeometryValidationError("COLLAR_LOOP_EXTRACTION_FAILED: zero loop perimeter")
    arc = np.concatenate(([0.0], np.cumsum(lengths)))
    targets = np.arange(count, dtype=float) * perimeter / count
    output = np.column_stack(
        [np.interp(targets, arc, closed[:, axis]) for axis in range(3)]
    )
    return output


def _phase_align(
    implicit: np.ndarray,
    explicit: np.ndarray,
) -> tuple[np.ndarray, int, bool, float]:
    best: tuple[float, np.ndarray, int, bool] | None = None
    for reverse in (False, True):
        candidate = explicit[::-1] if reverse else explicit
        for shift in range(len(candidate)):
            shifted = np.roll(candidate, -shift, axis=0)
            error = float(np.sum((implicit - shifted) ** 2))
            if best is None or error < best[0]:
                best = (error, shifted, shift, reverse)
    assert best is not None
    return best[1], best[2], best[3], best[0]


def _smoothstep(value: float, method: str) -> float:
    if method == "cubic":
        return 3.0 * value**2 - 2.0 * value**3
    return 6.0 * value**5 - 15.0 * value**4 + 10.0 * value**3


def _plane_polygon_area(
    points: np.ndarray,
    center: np.ndarray,
    normal: np.ndarray,
) -> float:
    first, second = _orthogonal_basis(normal)
    relative = points - center[None, :]
    planar = np.column_stack((relative @ first, relative @ second))
    return float(
        0.5
        * abs(
            np.dot(planar[:, 0], np.roll(planar[:, 1], -1))
            - np.dot(planar[:, 1], np.roll(planar[:, 0], -1))
        )
    )


def _radius_constrained_ring(
    points: np.ndarray,
    center: np.ndarray,
    tangent: np.ndarray,
    source_radius_um: float,
) -> tuple[np.ndarray, float]:
    """Project one transition ring and preserve the source-radius polygon area."""

    relative = points - center[None, :]
    radial = relative - (relative @ tangent)[:, None] * tangent[None, :]
    projected = center[None, :] + radial
    area_before = _plane_polygon_area(projected, center, tangent)
    target_area = (
        0.5
        * len(points)
        * np.sin(2.0 * np.pi / len(points))
        * source_radius_um**2
    )
    if area_before <= 0.0:
        raise GeometryValidationError(
            "COLLAR_LOOP_EXTRACTION_FAILED: transition ring has zero projected area"
        )
    constrained = center[None, :] + radial * np.sqrt(target_area / area_before)
    area_after = _plane_polygon_area(constrained, center, tangent)
    return constrained, float((area_after - target_area) / target_area)


def _zipper_faces(first_count: int, second_count: int) -> np.ndarray:
    faces: list[tuple[int, int, int]] = []
    first = 0
    second = 0
    tolerance = 1.0e-12
    while first < first_count or second < second_count:
        first_fraction = (first + 1) / first_count if first < first_count else np.inf
        second_fraction = (
            (second + 1) / second_count if second < second_count else np.inf
        )
        first_now = first % first_count
        second_now = first_count + second % second_count
        if abs(first_fraction - second_fraction) <= tolerance:
            first_next = (first + 1) % first_count
            second_next = first_count + (second + 1) % second_count
            faces.extend(
                (
                    (first_now, first_next, second_next),
                    (first_now, second_next, second_now),
                )
            )
            first += 1
            second += 1
        elif first_fraction < second_fraction:
            faces.append((first_now, (first + 1) % first_count, second_now))
            first += 1
        else:
            faces.append(
                (first_now, first_count + (second + 1) % second_count, second_now)
            )
            second += 1
    return np.asarray(faces, dtype=np.int64)


def _dimensionless_laplacian_roughness(mesh: trimesh.Trimesh) -> float:
    vertices = np.asarray(mesh.vertices, dtype=float)
    values: list[float] = []
    for vertex_id, neighbors in enumerate(mesh.vertex_neighbors):
        if not len(neighbors):
            continue
        neighbor_points = vertices[np.asarray(neighbors, dtype=np.int64)]
        scale = float(np.mean(np.linalg.norm(neighbor_points - vertices[vertex_id], axis=1)))
        if scale > 0:
            values.append(
                float(np.linalg.norm(vertices[vertex_id] - neighbor_points.mean(axis=0)))
                / scale
            )
    return float(np.mean(values)) if values else 0.0


def _constrained_taubin_transition(
    strip: trimesh.Trimesh,
    implicit_count: int,
    tube_sides: int,
    ring_count: int,
    ring_centers: list[np.ndarray],
    config: CFDLumenConfig,
) -> trimesh.Trimesh:
    vertices = np.asarray(strip.vertices, dtype=float).copy()
    neighbors = [np.asarray(row, dtype=np.int64) for row in strip.vertex_neighbors]
    movable = np.zeros(len(vertices), dtype=bool)
    ring_slices: list[slice] = []
    for ring_index in range(ring_count):
        start = implicit_count + ring_index * tube_sides
        ring_slice = slice(start, start + tube_sides)
        ring_slices.append(ring_slice)
        movable[ring_slice] = True
    target_rms: list[float] = []
    axis = np.asarray(ring_centers[-1] - ring_centers[0], dtype=float)
    axis /= np.linalg.norm(axis)
    for ring_slice, center in zip(ring_slices, ring_centers[1:-1]):
        relative = vertices[ring_slice] - center[None, :]
        relative -= (relative @ axis)[:, None] * axis[None, :]
        target_rms.append(float(np.sqrt(np.mean(np.sum(relative**2, axis=1)))))

    def step(weight: float) -> None:
        previous = vertices.copy()
        for vertex_id in np.flatnonzero(movable):
            if len(neighbors[vertex_id]):
                vertices[vertex_id] = previous[vertex_id] + weight * (
                    previous[neighbors[vertex_id]].mean(axis=0) - previous[vertex_id]
                )
        for ring_slice, center, radius in zip(
            ring_slices, ring_centers[1:-1], target_rms
        ):
            relative = vertices[ring_slice] - vertices[ring_slice].mean(axis=0)[None, :]
            relative -= (relative @ axis)[:, None] * axis[None, :]
            current = float(np.sqrt(np.mean(np.sum(relative**2, axis=1))))
            if current > 0:
                relative *= radius / current
            vertices[ring_slice] = center[None, :] + relative

    for _ in range(config.hybrid_transition.smoothing_iterations):
        step(config.hybrid_transition.smoothing_lambda)
        step(config.hybrid_transition.smoothing_mu)
    return trimesh.Trimesh(vertices=vertices, faces=np.asarray(strip.faces), process=False)


def build_transition_strip(
    implicit_mesh: trimesh.Trimesh,
    implicit_loop: np.ndarray,
    implicit_center: np.ndarray,
    explicit_mesh: trimesh.Trimesh,
    explicit_loop: np.ndarray,
    explicit_center: np.ndarray,
    branch: BranchGeometry,
    collar: JunctionCollar,
    tube_sides: int,
    config: CFDLumenConfig,
) -> tuple[trimesh.Trimesh, dict[str, Any]]:
    implicit_points = np.asarray(implicit_mesh.vertices[implicit_loop], dtype=float)
    explicit_points = np.asarray(explicit_mesh.vertices[explicit_loop], dtype=float)
    implicit_resampled = _resample_loop(implicit_points, tube_sides)
    explicit_resampled = _resample_loop(explicit_points, tube_sides)
    explicit_aligned, shift, reversed_order, phase_error = _phase_align(
        implicit_resampled, explicit_resampled
    )
    # Preserve the exact explicit boundary while allowing a VTK cap loop to
    # contain a few extra vertices on curved tubes. The transition uses N_loop
    # samples; a final order-preserving zipper welds those samples to the exact
    # tube boundary without moving explicit-region vertices.
    explicit_exact = explicit_points[::-1] if reversed_order else explicit_points
    explicit_start = int(
        np.argmin(np.linalg.norm(explicit_exact - explicit_aligned[0], axis=1))
    )
    explicit_exact = np.roll(explicit_exact, -explicit_start, axis=0)
    intermediate: list[np.ndarray] = []
    ring_centers = [np.asarray(implicit_center, dtype=float)]
    ring_source_radii: list[float] = []
    ring_area_relative_errors: list[float] = []
    for ring_index in range(1, config.hybrid_transition.transition_rings + 1):
        fraction = ring_index / (config.hybrid_transition.transition_rings + 1)
        weight = _smoothstep(fraction, config.hybrid_transition.smoothstep)
        distance = (
            (1.0 - fraction) * collar.explicit_cap_distance_um
            + fraction * collar.implicit_extent_um
        )
        center, source_radius, tangent = sample_from_junction(
            branch, collar.endpoint_index, distance
        )
        implicit_radial = implicit_resampled - implicit_center[None, :]
        explicit_radial = explicit_aligned - explicit_center[None, :]
        candidate = (
            center[None, :]
            + (1.0 - weight) * implicit_radial
            + weight * explicit_radial
        )
        constrained, area_error = _radius_constrained_ring(
            candidate, center, tangent, source_radius
        )
        intermediate.append(constrained)
        ring_source_radii.append(float(source_radius))
        ring_area_relative_errors.append(area_error)
        ring_centers.append(np.asarray(center, dtype=float))
    ring_centers.append(np.asarray(explicit_center, dtype=float))
    point_blocks = [implicit_points, *intermediate, explicit_exact]
    points = np.vstack(point_blocks)
    faces = _zipper_faces(len(implicit_points), tube_sides).tolist()
    previous_offset = len(implicit_points)
    for ring_index in range(max(0, len(intermediate) - 1)):
        next_offset = previous_offset + tube_sides
        for index in range(tube_sides):
            following = (index + 1) % tube_sides
            faces.extend(
                (
                    (
                        previous_offset + index,
                        previous_offset + following,
                        next_offset + following,
                    ),
                    (
                        previous_offset + index,
                        next_offset + following,
                        next_offset + index,
                    ),
                )
            )
        previous_offset = next_offset
    final_offset = len(implicit_points) + len(intermediate) * tube_sides
    final_zipper = _zipper_faces(tube_sides, len(explicit_exact))
    final_zipper = np.where(
        final_zipper < tube_sides,
        final_zipper + previous_offset,
        final_zipper - tube_sides + final_offset,
    )
    faces.extend(final_zipper.tolist())
    strip = trimesh.Trimesh(
        vertices=points,
        faces=np.asarray(faces, dtype=np.int64),
        process=False,
    )
    roughness_before = _dimensionless_laplacian_roughness(strip)
    if config.hybrid_transition.constrained_smoothing:
        strip = _constrained_taubin_transition(
            strip,
            len(implicit_points),
            tube_sides,
            len(intermediate),
            ring_centers,
            config,
        )
    roughness_after = _dimensionless_laplacian_roughness(strip)
    return strip, {
        "implicit_loop_vertex_count": int(len(implicit_points)),
        "resampled_loop_vertex_count": tube_sides,
        "explicit_loop_vertex_count": int(len(explicit_points)),
        "transition_ring_count": config.hybrid_transition.transition_rings,
        "transition_triangle_count": int(len(strip.faces)),
        "phase_shift": shift,
        "explicit_order_reversed": reversed_order,
        "phase_alignment_squared_error_um2": phase_error,
        "implicit_center_um": implicit_center.tolist(),
        "explicit_center_um": explicit_center.tolist(),
        "transition_length_um": float(np.linalg.norm(explicit_center - implicit_center)),
        "smoothstep": config.hybrid_transition.smoothstep,
        "ring_source_radius_um": ring_source_radii,
        "maximum_ring_area_absolute_relative_error": max(
            map(abs, ring_area_relative_errors), default=None
        ),
        "constrained_smoothing_applied": config.hybrid_transition.constrained_smoothing,
        "transition_roughness_before_smoothing": roughness_before,
        "transition_roughness_after_smoothing": roughness_after,
        "transition_roughness_smoothing_ratio": (
            roughness_after / roughness_before if roughness_before > 0 else None
        ),
    }


def combine_stitched_surfaces(
    core_meshes: list[trimesh.Trimesh],
    transition_meshes: list[trimesh.Trimesh],
    explicit_meshes: list[trimesh.Trimesh],
) -> tuple[trimesh.Trimesh, np.ndarray]:
    meshes = [*core_meshes, *transition_meshes, *explicit_meshes]
    labels = np.concatenate(
        [
            *[
                np.full(len(mesh.faces), JUNCTION_CORE_FACE, dtype=np.uint8)
                for mesh in core_meshes
            ],
            *[
                np.full(len(mesh.faces), TRANSITION_COLLAR_FACE, dtype=np.uint8)
                for mesh in transition_meshes
            ],
            *[
                np.full(len(mesh.faces), EXPLICIT_BRANCH_FACE, dtype=np.uint8)
                for mesh in explicit_meshes
            ],
        ]
    )
    combined = trimesh.util.concatenate(meshes)
    combined.merge_vertices(digits_vertex=8)
    keep = combined.nondegenerate_faces(height=1.0e-12) & combined.unique_faces()
    labels = labels[keep]
    combined.update_faces(keep)
    combined.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(combined, multibody=True)
    if combined.is_watertight and combined.volume < 0:
        combined.invert()
    return combined, labels
