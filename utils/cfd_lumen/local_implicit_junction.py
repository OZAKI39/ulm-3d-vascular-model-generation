"""Local polyball junction patches and source-preserving collar construction."""

from __future__ import annotations

import time
from dataclasses import replace
from typing import Any

import numpy as np
import trimesh
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .types import (
    BranchGeometry,
    GeometryValidationError,
    JunctionCollar,
    LocalJunctionPatch,
)


def _geometry_conflict(
    code: str,
    message: str,
    *,
    junction_node_id: int | None = None,
    branch_id: int | None = None,
) -> GeometryValidationError:
    error = GeometryValidationError(f"{code}: {message}")
    error.code = code  # type: ignore[attr-defined]
    error.junction_node_id = junction_node_id  # type: ignore[attr-defined]
    error.branch_id = branch_id  # type: ignore[attr-defined]
    return error


def _endpoint_view(branch: BranchGeometry, endpoint_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if endpoint_index == 0:
        return branch.points_um, branch.radius_um, branch.arc_length_um
    length = float(branch.arc_length_um[-1])
    return branch.points_um[::-1], branch.radius_um[::-1], length - branch.arc_length_um[::-1]


def sample_from_junction(
    branch: BranchGeometry,
    endpoint_index: int,
    distance_um: float,
) -> tuple[np.ndarray, float, np.ndarray]:
    points, radii, arc = _endpoint_view(branch, endpoint_index)
    distance = float(np.clip(distance_um, 0.0, arc[-1]))
    index = min(max(int(np.searchsorted(arc, distance)), 1), len(points) - 1)
    span = float(arc[index] - arc[index - 1])
    fraction = 0.0 if span <= 0 else (distance - float(arc[index - 1])) / span
    point = (1.0 - fraction) * points[index - 1] + fraction * points[index]
    radius = float((1.0 - fraction) * radii[index - 1] + fraction * radii[index])
    tangent = points[index] - points[index - 1]
    tangent /= np.linalg.norm(tangent)
    return np.asarray(point), radius, np.asarray(tangent)


def _junction_incident_branches(
    branches: list[BranchGeometry],
    junction_node_id: int,
) -> list[tuple[BranchGeometry, int]]:
    incident: list[tuple[BranchGeometry, int]] = []
    for branch in branches:
        if branch.local_node_ids[0] == junction_node_id:
            incident.append((branch, 0))
        elif branch.local_node_ids[-1] == junction_node_id:
            incident.append((branch, -1))
    return incident


def define_junction_collars(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    config: CFDLumenConfig,
    *,
    junction_node_ids: set[int] | None = None,
) -> dict[int, list[JunctionCollar]]:
    """Place collars outside the core without crossing another topology node or CUT_PORT."""

    degree = np.bincount(np.asarray(roi.local_edges).ravel(), minlength=roi.node_count)
    cut_ids = {int(port.local_node_id): port for port in roi.cut_ports}
    output: dict[int, list[JunctionCollar]] = {}
    for node_id in map(int, np.flatnonzero(degree >= 3)):
        if junction_node_ids is not None and node_id not in junction_node_ids:
            continue
        incident = _junction_incident_branches(branches, node_id)
        distances: dict[int, float] = {}
        maxima: dict[int, float] = {}
        for branch, endpoint_index in incident:
            endpoint_radius = float(branch.radius_um[endpoint_index])
            diameter = 2.0 * endpoint_radius
            overlap = config.junction.overlap_diameters * diameter
            other_node = branch.local_node_ids[-1] if endpoint_index == 0 else branch.local_node_ids[0]
            reserved = max(0.5 * float(branch.radius_um[-1 if endpoint_index == 0 else 0]), 1.0e-6)
            if int(other_node) in cut_ids:
                reserved += config.ports.overlap_diameters * 2.0 * float(
                    cut_ids[int(other_node)].radius_at_cut_um
                )
            maximum = branch.length_um - overlap - reserved
            base = config.junction.collar_diameters * diameter
            if maximum <= overlap or base > maximum:
                conflict = "JUNCTION_PORT_REGION_CONFLICT" if int(other_node) in cut_ids else "JUNCTION_COLLAR_TOPOLOGY_CONFLICT"
                raise _geometry_conflict(
                    conflict,
                    f"junction {node_id}, branch {branch.branch_id}, "
                    f"required collar {base:.6g} um, available {maximum:.6g} um",
                    junction_node_id=node_id,
                    branch_id=branch.branch_id,
                )
            distances[branch.branch_id] = base
            maxima[branch.branch_id] = maximum

        # Acute branches may still overlap at 2D. Move only the involved collars
        # outward in 0.25D steps, never beyond the next topology endpoint.  The
        # complete Boolean overlap interval must be tube-like, so separation is
        # checked both at the collar and at its inner explicit-cap plane.
        for _ in range(64):
            samples: dict[int, tuple[tuple[np.ndarray, float, np.ndarray], tuple[np.ndarray, float, np.ndarray]]] = {}
            for branch, endpoint_index in incident:
                distance = distances[branch.branch_id]
                collar_sample = sample_from_junction(branch, endpoint_index, distance)
                collar_radius = collar_sample[1]
                overlap = config.junction.overlap_diameters * 2.0 * collar_radius
                cap_distance = max(distance - overlap, 0.25 * collar_radius)
                samples[branch.branch_id] = (
                    collar_sample,
                    sample_from_junction(branch, endpoint_index, cap_distance),
                )
            failing: set[int] = set()
            for first in range(len(incident)):
                for second in range(first + 1, len(incident)):
                    branch_a = incident[first][0]
                    branch_b = incident[second][0]
                    for sample_index in (0, 1):
                        point_a, radius_a, _ = samples[branch_a.branch_id][sample_index]
                        point_b, radius_b, _ = samples[branch_b.branch_id][sample_index]
                        radius_sum = radius_a + radius_b
                        tolerance = config.junction.separation_tolerance_fraction * radius_sum
                        if float(np.linalg.norm(point_a - point_b)) <= radius_sum + tolerance:
                            failing.update((branch_a.branch_id, branch_b.branch_id))
                            break
            if not failing:
                break
            for branch, endpoint_index in incident:
                if branch.branch_id not in failing:
                    continue
                diameter = 2.0 * float(branch.radius_um[endpoint_index])
                following = distances[branch.branch_id] + 0.25 * diameter
                if following > maxima[branch.branch_id]:
                    raise _geometry_conflict(
                        "JUNCTION_COLLAR_SEPARATION_FAILED",
                        f"junction {node_id}, branch {branch.branch_id} cannot separate "
                        "before next topology node",
                        junction_node_id=node_id,
                        branch_id=branch.branch_id,
                    )
                distances[branch.branch_id] = following
        else:
            raise _geometry_conflict(
                "JUNCTION_COLLAR_SEPARATION_FAILED",
                f"junction {node_id}",
                junction_node_id=node_id,
            )

        collars: list[JunctionCollar] = []
        for branch, endpoint_index in incident:
            distance = distances[branch.branch_id]
            point, radius, tangent = sample_from_junction(branch, endpoint_index, distance)
            overlap = config.junction.overlap_diameters * 2.0 * radius
            collars.append(
                JunctionCollar(
                    junction_node_id=node_id,
                    branch_id=branch.branch_id,
                    endpoint_index=endpoint_index,
                    collar_distance_um=distance,
                    overlap_length_um=overlap,
                    implicit_extent_um=distance + overlap,
                    explicit_cap_distance_um=max(distance - overlap, 0.25 * radius),
                    collar_position_um=point,
                    collar_radius_um=radius,
                    tangent_away=tangent,
                )
            )
        output[node_id] = collars

    # Patches from opposite ends of an inter-junction branch must remain local.
    collars_by_branch: dict[int, list[JunctionCollar]] = {}
    for collars in output.values():
        for collar in collars:
            collars_by_branch.setdefault(collar.branch_id, []).append(collar)
    for branch in branches:
        pair = collars_by_branch.get(branch.branch_id, [])
        if len(pair) == 2 and sum(item.implicit_extent_um for item in pair) >= branch.length_um:
            raise _geometry_conflict(
                "JUNCTION_PATCH_OVERLAP_CONFLICT",
                f"branch {branch.branch_id} is too short for two local patches",
                branch_id=branch.branch_id,
            )
    return output


def _dense_local_samples(
    branch: BranchGeometry,
    collar: JunctionCollar,
    spacing: float,
) -> tuple[np.ndarray, np.ndarray]:
    count = max(2, int(np.ceil(collar.implicit_extent_um / spacing)) + 1)
    distances = np.linspace(0.0, collar.implicit_extent_um, count)
    points: list[np.ndarray] = []
    radii: list[float] = []
    for distance in distances:
        point, radius, _ = sample_from_junction(branch, collar.endpoint_index, float(distance))
        points.append(point)
        radii.append(radius)
    return np.asarray(points, dtype=float), np.asarray(radii, dtype=float)


def _clean_local_mesh(mesh: trimesh.Trimesh) -> trimesh.Trimesh:
    cleaned = mesh.copy()
    cleaned.update_faces(cleaned.nondegenerate_faces())
    cleaned.update_faces(cleaned.unique_faces())
    cleaned.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(cleaned, multibody=True)
    if cleaned.is_watertight and cleaned.volume < 0:
        cleaned.invert()
    return cleaned


def build_local_junction_patch(
    roi: ROIRecord,
    branches_by_id: dict[int, BranchGeometry],
    junction_node_id: int,
    collars: list[JunctionCollar],
    config: CFDLumenConfig,
    *,
    cells_across_min_diameter: int | None = None,
    controlled: bool = False,
) -> LocalJunctionPatch:
    """Build a local float32 swept-sphere/polyball zero-level surface."""

    started = time.perf_counter()
    cells = int(cells_across_min_diameter or config.junction.implicit.cells_across_min_diameter)
    minimum_diameter = 2.0 * min(collar.collar_radius_um for collar in collars)
    grid_spacing = minimum_diameter / cells
    sample_spacing = 0.5 * grid_spacing
    point_blocks: list[np.ndarray] = []
    radius_blocks: list[np.ndarray] = []
    for collar in collars:
        points, radii = _dense_local_samples(
            branches_by_id[collar.branch_id], collar, sample_spacing
        )
        point_blocks.append(points)
        radius_blocks.append(radii)
    sample_points = np.vstack(point_blocks)
    sample_radii = np.concatenate(radius_blocks)
    junction_center = np.asarray(roi.local_node_positions_um[junction_node_id], dtype=float)
    junction_radius = float(roi.local_node_radius_um[junction_node_id])
    max_radius = float(sample_radii.max())
    padding = config.junction.bbox_padding_radius * max_radius
    minimum = np.min(sample_points - sample_radii[:, None], axis=0) - padding
    maximum = np.max(sample_points + sample_radii[:, None], axis=0) + padding
    minimum = np.minimum(minimum, junction_center - junction_radius - padding)
    maximum = np.maximum(maximum, junction_center + junction_radius + padding)
    axes = [
        np.arange(minimum[axis], maximum[axis] + 0.5 * grid_spacing, grid_spacing)
        for axis in range(3)
    ]
    nx, ny, nz = map(len, axes)
    grid_cells = int(nx * ny * nz)
    if grid_cells > config.junction.implicit.max_grid_cells:
        raise GeometryValidationError(
            f"Local junction {junction_node_id} grid has {grid_cells:,} cells, exceeding "
            f"{config.junction.implicit.max_grid_cells:,}"
        )
    field = np.empty((nz, ny, nx), dtype=np.float32)
    tree = cKDTree(sample_points)
    k = min(config.junction.implicit.k_nearest, len(sample_points))
    controlled_trees = [cKDTree(points) for points in point_blocks]
    plane_size = nx * ny
    slab = max(1, config.junction.implicit.chunk_size // max(1, plane_size))
    for z_start in range(0, nz, slab):
        z_end = min(nz, z_start + slab)
        zz, yy, xx = np.meshgrid(
            axes[2][z_start:z_end], axes[1], axes[0], indexing="ij"
        )
        query = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        if controlled:
            phi = np.linalg.norm(query - junction_center[None, :], axis=1) - junction_radius
            for collar, branch_tree, radii in zip(
                collars, controlled_trees, radius_blocks
            ):
                branch_k = min(config.junction.implicit.k_nearest, len(radii))
                distances, indices = branch_tree.query(query, k=branch_k, workers=-1)
                if branch_k == 1:
                    branch_phi = distances - radii[indices]
                else:
                    branch_phi = np.min(distances - radii[indices], axis=1)
                epsilon = (
                    config.junction.implicit.direction_clip_epsilon_radius_fraction
                    * junction_radius
                )
                valid = (query - junction_center[None, :]) @ collar.tangent_away >= -epsilon
                phi = np.minimum(phi, np.where(valid, branch_phi, np.inf))
        else:
            distances, indices = tree.query(query, k=k, workers=-1)
            if k == 1:
                phi = distances - sample_radii[indices]
            else:
                phi = np.min(distances - sample_radii[indices], axis=1)
        field[z_start:z_end] = phi.reshape((z_end - z_start, ny, nx)).astype(np.float32)
    if not (float(field.min()) < 0.0 < float(field.max())):
        raise GeometryValidationError(f"Local junction {junction_node_id} field does not bracket zero")
    vertices_zyx, faces, _, _ = marching_cubes(
        field, level=0.0, spacing=(grid_spacing, grid_spacing, grid_spacing)
    )
    vertices = vertices_zyx[:, ::-1] + minimum
    raw = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    clean = _clean_local_mesh(raw)
    if not clean.is_watertight:
        raise GeometryValidationError(f"Local junction {junction_node_id} patch is not watertight")
    metadata: dict[str, Any] = {
        "junction_node_id": junction_node_id,
        "method": (
            "controlled_direction_clipped_variable_radius_polyball"
            if controlled
            else "variable_radius_polyball"
        ),
        "controlled_local_implicit": controlled,
        "junction_core_field": "norm(x - x_J) - r_J" if controlled else None,
        "direction_clip_epsilon_um": (
            config.junction.implicit.direction_clip_epsilon_radius_fraction
            * junction_radius
            if controlled
            else None
        ),
        "smooth_union": False,
        "dtype": "float32",
        "cells_across_min_diameter": cells,
        "grid_spacing_um": grid_spacing,
        "grid_dimensions_xyz": [nx, ny, nz],
        "grid_cell_count": grid_cells,
        "bbox_min_um": minimum.tolist(),
        "bbox_max_um": maximum.tolist(),
        "centerline_sample_count": int(len(sample_points)),
        "raw_triangle_count": int(len(raw.faces)),
        "clean_triangle_count": int(len(clean.faces)),
        "degenerate_triangles_removed": int(len(raw.faces) - len(clean.faces)),
        "local_volume_um3": float(abs(clean.volume)),
        "runtime_s": time.perf_counter() - started,
        "collars": [
            {
                "branch_id": collar.branch_id,
                "endpoint_index": collar.endpoint_index,
                "collar_distance_um": collar.collar_distance_um,
                "overlap_length_um": collar.overlap_length_um,
                "implicit_extent_um": collar.implicit_extent_um,
                "explicit_cap_distance_um": collar.explicit_cap_distance_um,
                "collar_position_um": collar.collar_position_um.tolist(),
                "collar_radius_um": collar.collar_radius_um,
            }
            for collar in collars
        ],
    }
    return LocalJunctionPatch(
        junction_node_id=junction_node_id,
        raw_mesh=raw,
        clean_mesh=clean,
        collars=collars,
        centerline_points_um=sample_points,
        centerline_radius_um=sample_radii,
        metadata=metadata,
    )


def trim_explicit_branches(
    branches: list[BranchGeometry],
    collars_by_junction: dict[int, list[JunctionCollar]],
) -> tuple[list[BranchGeometry], list[dict[str, Any]]]:
    by_branch: dict[int, list[JunctionCollar]] = {}
    for collars in collars_by_junction.values():
        for collar in collars:
            by_branch.setdefault(collar.branch_id, []).append(collar)
    output: list[BranchGeometry] = []
    cap_planes: list[dict[str, Any]] = []
    for branch in branches:
        lower = 0.0
        upper = float(branch.arc_length_um[-1])
        for collar in by_branch.get(branch.branch_id, []):
            if collar.endpoint_index == 0:
                lower = max(lower, collar.explicit_cap_distance_um)
            else:
                upper = min(upper, float(branch.arc_length_um[-1]) - collar.explicit_cap_distance_um)
        if upper - lower <= 1.0e-8:
            raise _geometry_conflict(
                "JUNCTION_COLLAR_TOPOLOGY_CONFLICT",
                f"branch {branch.branch_id} has no explicit segment",
                branch_id=branch.branch_id,
            )
        inside = branch.arc_length_um[(branch.arc_length_um > lower) & (branch.arc_length_um < upper)]
        targets = np.concatenate(([lower], inside, [upper]))
        points = np.column_stack(
            [np.interp(targets, branch.arc_length_um, branch.points_um[:, axis]) for axis in range(3)]
        )
        radii = np.interp(targets, branch.arc_length_um, branch.radius_um)
        trimmed = replace(branch)
        trimmed.points_um = points
        trimmed.radius_um = radii
        trimmed.arc_length_um = targets - lower
        output.append(trimmed)
        for collar in by_branch.get(branch.branch_id, []):
            if collar.endpoint_index == 0:
                center = points[0]
                tangent = points[1] - points[0]
                radius = radii[0]
            else:
                center = points[-1]
                tangent = points[-2] - points[-1]
                radius = radii[-1]
            tangent /= np.linalg.norm(tangent)
            cap_planes.append(
                {
                    "junction_node_id": collar.junction_node_id,
                    "branch_id": branch.branch_id,
                    "center_um": np.asarray(center).tolist(),
                    "normal_away_from_junction": np.asarray(tangent).tolist(),
                    "radius_um": float(radius),
                }
            )
    return output, cap_planes
