"""Source-preserving branch extraction, validation, and arc-length resampling."""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
from scipy.interpolate import PchipInterpolator

from utils.sampling.sampling_types import ROIRecord
from utils.sampling.structural_features import _branch_paths

from .config import CFDLumenConfig
from .types import BranchGeometry, GeometryValidationError


def _roi_graph(roi: ROIRecord) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(int(value) for value in roi.local_node_ids)
    graph.add_edges_from((int(first), int(second)) for first, second in roi.local_edges)
    return graph


def validate_and_extract_branches(
    roi: ROIRecord,
    config: CFDLumenConfig,
) -> tuple[list[BranchGeometry], dict[str, Any]]:
    """Validate all reconstruction inputs without mutating the saved source ROI."""

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
    invalid_coordinate = np.flatnonzero(~np.all(np.isfinite(positions), axis=1))
    for index in invalid_coordinate.tolist():
        failures.append({"reason": "nonfinite_coordinate", "node_id": int(node_ids[index])})
    invalid_radius = np.flatnonzero(~np.isfinite(radii) | (radii <= 0.0))
    for index in invalid_radius.tolist():
        failures.append(
            {
                "reason": "invalid_radius",
                "node_id": int(node_ids[index]),
                "radius_um": float(radii[index]),
            }
        )
    if len(edges):
        if int(edges.min()) < 0 or int(edges.max()) >= len(node_ids):
            failures.append({"reason": "edge_node_index_out_of_range"})
        lengths = np.linalg.norm(positions[edges[:, 1]] - positions[edges[:, 0]], axis=1)
    else:
        lengths = np.empty(0, dtype=float)
        failures.append({"reason": "roi_has_no_edges"})
    invalid_edges = np.flatnonzero(
        ~np.isfinite(lengths) | (lengths <= config.geometry.minimum_edge_length_um)
    )
    for edge_index in invalid_edges.tolist():
        failures.append(
            {
                "reason": "zero_or_near_zero_edge",
                "edge_id": int(roi.local_edge_ids[edge_index]),
                "global_edge_id": int(roi.local_edge_global_ids[edge_index]),
                "length_um": float(lengths[edge_index]),
            }
        )

    graph = _roi_graph(roi)
    component_count = nx.number_connected_components(graph) if graph else 0
    duplicate_edges = len(edges) - len({tuple(sorted(map(int, edge))) for edge in edges})
    if duplicate_edges:
        failures.append({"reason": "duplicate_local_edges", "count": duplicate_edges})
    if component_count != 1:
        failures.append({"reason": "roi_not_connected", "component_count": component_count})
    if graph and component_count == 1:
        branch_membership: dict[int, list[int]] = {int(node): [] for node in graph}
        for branch_id, path in enumerate(_branch_paths(graph)):
            for node in path:
                branch_membership[int(node)].append(branch_id)
        for failure in failures:
            if failure.get("reason") == "invalid_radius" and "node_id" in failure:
                failure["branch_id"] = branch_membership.get(int(failure["node_id"]), [])
    if failures:
        error = GeometryValidationError(
            f"ROI {roi.roi_id} failed pre-geometry validation: {failures[0]['reason']}"
        )
        error.failures = failures  # type: ignore[attr-defined]
        raise error

    edge_index_by_nodes: dict[tuple[int, int], int] = {}
    for edge_index, (first, second) in enumerate(edges):
        edge_index_by_nodes[tuple(sorted((int(first), int(second))))] = edge_index
    paths = _branch_paths(graph)
    branches: list[BranchGeometry] = []
    for branch_id, path in enumerate(paths):
        source_edge_ids: list[int] = []
        for first, second in zip(path[:-1], path[1:]):
            edge_index = edge_index_by_nodes[tuple(sorted((int(first), int(second))))]
            source_edge_ids.append(int(roi.local_edge_global_ids[edge_index]))
        global_nodes = tuple(int(roi.local_node_global_ids[node]) for node in path)
        branches.append(
            BranchGeometry(
                branch_id=branch_id,
                local_node_ids=tuple(map(int, path)),
                source_global_nodes=global_nodes,
                source_global_edges=tuple(source_edge_ids),
                raw_points_um=positions[path].copy(),
                raw_radius_um=radii[path].copy(),
            )
        )

    covered = [edge for branch in branches for edge in branch.source_global_edges]
    expected = list(map(int, roi.local_edge_global_ids))
    if sorted(covered) != sorted(expected) or len(covered) != len(expected):
        raise GeometryValidationError(
            f"ROI {roi.roi_id} branch extraction did not preserve every source edge exactly once"
        )
    report = {
        "roi_id": roi.roi_id,
        "geometry_unit": "um",
        "node_count": int(len(node_ids)),
        "edge_count": int(len(edges)),
        "branch_count": len(branches),
        "bifurcation_count": int(sum(graph.degree(node) >= 3 for node in graph)),
        "connected_component_count": component_count,
        "cycle_rank": int(graph.number_of_edges() - graph.number_of_nodes() + component_count),
        "cut_port_count": len(roi.cut_ports),
        "true_terminal_count": len(roi.true_terminal_local_ids),
        "radius_min_um": float(radii.min()),
        "radius_median_um": float(np.median(radii)),
        "radius_max_um": float(radii.max()),
        "centerline_total_length_um": float(lengths.sum()),
        "source_topology_preserved": True,
        "source_edge_coverage_exactly_once": True,
        "validation_failures": [],
    }
    return branches, report


def _light_smooth(points: np.ndarray, window: int) -> np.ndarray:
    if len(points) <= 2 or window <= 1:
        return points.copy()
    window = max(3, int(window) | 1)
    radius = window // 2
    padded = np.pad(points, ((radius, radius), (0, 0)), mode="edge")
    kernel = np.full(window, 1.0 / window)
    smoothed = np.column_stack(
        [np.convolve(padded[:, axis], kernel, mode="valid") for axis in range(3)]
    )
    smoothed[0] = points[0]
    smoothed[-1] = points[-1]
    return smoothed


def resample_branch(branch: BranchGeometry, config: CFDLumenConfig) -> BranchGeometry:
    points = branch.raw_points_um
    if config.geometry.centerline_smoothing:
        points = _light_smooth(points, config.geometry.smoothing_window_points)
    raw_step = np.linalg.norm(np.diff(points, axis=0), axis=1)
    source_s = np.concatenate(([0.0], np.cumsum(raw_step)))
    if source_s[-1] <= config.geometry.minimum_edge_length_um:
        raise GeometryValidationError(f"Branch {branch.branch_id} has no positive arc length")
    radius_interpolator = PchipInterpolator(source_s, branch.raw_radius_um, extrapolate=False)
    sample_s = [0.0]
    while sample_s[-1] < source_s[-1]:
        current_radius = float(radius_interpolator(sample_s[-1]))
        spacing = float(
            np.clip(
                config.geometry.resample_radius_fraction * current_radius,
                config.geometry.min_resample_spacing_um,
                config.geometry.max_resample_spacing_um,
            )
        )
        following = min(sample_s[-1] + spacing, float(source_s[-1]))
        if following - sample_s[-1] <= config.geometry.minimum_edge_length_um:
            break
        sample_s.append(following)
    if sample_s[-1] < source_s[-1]:
        sample_s.append(float(source_s[-1]))
    target_s = np.asarray(sample_s, dtype=float)
    target_points = np.column_stack(
        [np.interp(target_s, source_s, points[:, axis]) for axis in range(3)]
    )
    target_radius = np.asarray(radius_interpolator(target_s), dtype=float)
    if np.any(~np.isfinite(target_radius)) or np.any(target_radius <= 0):
        raise GeometryValidationError(
            f"Branch {branch.branch_id} acquired a non-positive/non-finite radius during resampling"
        )
    branch.points_um = target_points
    branch.radius_um = target_radius
    branch.arc_length_um = target_s
    return branch


def resample_branches(
    branches: list[BranchGeometry],
    config: CFDLumenConfig,
) -> list[BranchGeometry]:
    return [resample_branch(branch, config) for branch in branches]
