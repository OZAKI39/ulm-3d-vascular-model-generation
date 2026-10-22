"""Display-only interpolating tubes; saved ROI positions/radii stay untouched."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

from ..sampling.sampling_types import ROIRecord


TUBE_SIDES = 32
MIN_EDGE_SUBDIVISIONS = 8
MAX_EDGE_SUBDIVISIONS = 64


def _branch_paths(edges: np.ndarray) -> tuple[list[list[int]], dict[int, int]]:
    """Trace every edge once; do not smooth across a bifurcation."""

    adjacency: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for edge_id, (start, end) in enumerate(edges):
        adjacency[int(start)].append((edge_id, int(end)))
        adjacency[int(end)].append((edge_id, int(start)))
    degree = {node: len(links) for node, links in adjacency.items()}
    visited: set[int] = set()
    paths: list[list[int]] = []

    def trace(start: int, edge_id: int, neighbor: int) -> None:
        path = [start]
        while edge_id not in visited:
            visited.add(edge_id)
            path.append(neighbor)
            if degree[neighbor] != 2:
                break
            next_links = [link for link in adjacency[neighbor] if link[0] not in visited]
            if not next_links:
                break
            edge_id, neighbor = next_links[0]
        paths.append(path)

    for node in sorted(adjacency):
        if degree[node] != 2:
            for edge_id, neighbor in adjacency[node]:
                if edge_id not in visited:
                    trace(node, edge_id, neighbor)
    # Also handle closed degree-two loops without dropping edges or hanging.
    for edge_id, (start, end) in enumerate(edges):
        if edge_id not in visited:
            trace(int(start), edge_id, int(end))
    return paths, degree


def _interpolate_branch(
    points: np.ndarray, radii: np.ndarray, *, closed: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """C1 interpolation with exact source knots and no radius overshoot.

    A local Hermite curve passes through every centerline point. Unit bisector
    tangents round bends without moving source points. PCHIP interpolates radii
    independently: it does not average, filter, or fit away their original values.
    """

    deltas = np.diff(points, axis=0)
    lengths = np.linalg.norm(deltas, axis=1)
    if np.any(lengths <= 0):
        raise ValueError("Cannot interpolate a zero-length ROI edge while preserving its radius samples")
    arc = np.concatenate(([0.0], np.cumsum(lengths)))
    directions = deltas / lengths[:, None]
    tangents = np.vstack((directions[0], directions[:-1] + directions[1:], directions[-1]))
    if closed:
        tangents[0] = tangents[-1] = directions[-1] + directions[0]
    norms = np.linalg.norm(tangents, axis=1)
    tangents = np.divide(tangents, norms[:, None], out=np.zeros_like(tangents), where=norms[:, None] > 0)
    centerline = CubicHermiteSpline(arc, points, tangents, axis=0, extrapolate=False)
    if closed:
        radius_curve = PchipInterpolator(
            np.concatenate(([arc[-2] - arc[-1]], arc, [arc[1] + arc[-1]])),
            np.concatenate(([radii[-2]], radii, [radii[1]])),
            extrapolate=False,
        )
    else:
        radius_curve = PchipInterpolator(arc, radii, extrapolate=False)
    subdivisions = np.clip(
        np.ceil(lengths / (0.5 * np.minimum(radii[:-1], radii[1:]))),
        MIN_EDGE_SUBDIVISIONS,
        MAX_EDGE_SUBDIVISIONS,
    ).astype(np.int64)
    sample_arc = np.concatenate(
        [np.linspace(arc[i], arc[i + 1], int(count), endpoint=False) for i, count in enumerate(subdivisions)]
        + [arc[-1:]]
    )
    knots = np.concatenate(([0], np.cumsum(subdivisions)))
    sampled_points = centerline(sample_arc)
    sampled_radii = radius_curve(sample_arc)
    # Include each original sample explicitly, without even interpolation roundoff.
    sampled_points[knots] = points
    sampled_radii[knots] = radii
    return sampled_points, sampled_radii, knots, sample_arc


def build_roi_display_tubes(roi: ROIRecord) -> tuple[Any, Any] | None:
    """Return interpolated display centerlines and tubes with true-end caps only.

    This is a visualization, not a boolean-unioned lumen or simulation mesh.
    At bifurcations the original branches meet without interior end disks.
    """

    import pyvista as pv

    edges = np.asarray(roi.local_edges, dtype=np.int64)
    if not len(edges):
        return None
    node_ids = np.asarray(roi.local_node_ids, dtype=np.int64)
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    if not np.all(np.isfinite(positions)) or not np.all(np.isfinite(radii) & (radii > 0)):
        raise ValueError("ROI display requires finite coordinates and positive source radii")
    index_by_id = {int(node_id): index for index, node_id in enumerate(node_ids)}
    edge_indices = np.asarray([[index_by_id[int(node)] for node in edge] for edge in edges])
    if not np.allclose(positions[edge_indices], roi.local_edge_points_um, rtol=1e-12, atol=1e-10):
        raise ValueError("ROI node and edge coordinates disagree")
    if not np.allclose(radii[edge_indices], roi.local_edge_radius_um, rtol=1e-12, atol=1e-10):
        raise ValueError("ROI node and edge radii disagree; refusing to average source values")

    paths, degree = _branch_paths(edges)
    point_chunks: list[np.ndarray] = []
    radius_chunks: list[np.ndarray] = []
    source_chunks: list[np.ndarray] = []
    arc_chunks: list[np.ndarray] = []
    cells: list[np.ndarray] = []
    end_caps: list[tuple[int, np.ndarray]] = []
    offset = 0
    for path in paths:
        indices = np.asarray([index_by_id[node] for node in path])
        closed = path[0] == path[-1]
        points, radius, knots, arc = _interpolate_branch(positions[indices], radii[indices], closed=closed)
        source_nodes = np.full(len(points), -1, dtype=np.int64)
        source_nodes[knots] = path
        if closed:
            # A closed polyline reuses its first point ID, avoiding duplicate rims.
            points, radius, source_nodes, arc = points[:-1], radius[:-1], source_nodes[:-1], arc[:-1]
        else:
            if degree[path[0]] == 1:
                end_caps.append((offset, points[0] - points[1]))
            if degree[path[-1]] == 1:
                end_caps.append((offset + len(points) - 1, points[-1] - points[-2]))
        point_chunks.append(points)
        radius_chunks.append(radius)
        source_chunks.append(source_nodes)
        arc_chunks.append(arc)
        cell = np.arange(offset, offset + len(points), dtype=np.int64)
        if closed:
            cell = np.append(cell, offset)
        cells.append(np.concatenate(([len(cell)], cell)))
        offset += len(points)

    centerlines = pv.PolyData(np.concatenate(point_chunks))
    centerlines.lines = np.concatenate(cells)
    centerlines.point_data["radius_um"] = np.concatenate(radius_chunks)
    centerlines.point_data["diameter_um"] = 2.0 * centerlines.point_data["radius_um"]
    centerlines.point_data["source_local_node_id"] = np.concatenate(source_chunks)
    centerlines.point_data["source_arc_um"] = np.concatenate(arc_chunks)
    centerlines.point_data["display_sample_id"] = np.arange(centerlines.n_points, dtype=np.int64)
    tubes = centerlines.tube(scalars="radius_um", absolute=True, n_sides=TUBE_SIDES, capping=False)

    # VTK shares rings along each complete polyline. Cap only degree-one ends;
    # ordinary source nodes and bifurcations never receive an internal disk.
    faces: list[np.ndarray] = []
    sample_ids = np.asarray(tubes.point_data["display_sample_id"])
    for sample_id, outward in end_caps:
        ring = np.flatnonzero(sample_ids == sample_id)
        if len(ring) != TUBE_SIDES:
            raise ValueError("Tube filter did not preserve an endpoint's radius ring")
        rim = np.asarray(tubes.points[ring])
        if np.dot(np.cross(rim[1] - rim[0], rim[2] - rim[0]), outward) < 0:
            ring = ring[::-1]
        faces.append(np.concatenate(([len(ring)], ring)))
    if faces:
        tubes.faces = np.concatenate(faces)
    return centerlines, tubes
