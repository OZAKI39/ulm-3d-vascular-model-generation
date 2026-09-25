"""Loss-aware aggregation of a voxel skeleton into vascular branches."""

from __future__ import annotations

import logging
from collections import Counter
from itertools import combinations

import networkx as nx
import numpy as np
from scipy import ndimage

from ..config import HierarchicalGraphConfig
from .cycles import identify_cycles
from .geometry import cumulative_length, indices_to_lps, smooth_and_measure, tortuosity
from .model import (
    BranchRecord,
    HierarchicalGraphResult,
    JunctionGeometry,
    NodeRecord,
)
from .representations import build_branch_as_node_graph, build_junction_graph


def _neighbor_offsets(connectivity: int) -> list[tuple[int, int, int]]:
    limit = {6: 1, 18: 2, 26: 3}[connectivity]
    offsets: list[tuple[int, int, int]] = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == dy == dz == 0:
                    continue
                if abs(dx) + abs(dy) + abs(dz) <= limit:
                    offsets.append((dx, dy, dz))
    return offsets


def _ndimage_structure(connectivity: int) -> np.ndarray:
    return ndimage.generate_binary_structure(3, {6: 1, 18: 2, 26: 3}[connectivity])


def _voxel_adjacency(
    skeleton: np.ndarray, connectivity: int
) -> tuple[np.ndarray, list[list[int]]]:
    indices = np.argwhere(skeleton).astype(np.int32, copy=False)
    lookup = {tuple(int(value) for value in index): voxel_id for voxel_id, index in enumerate(indices)}
    offsets = _neighbor_offsets(connectivity)
    adjacency: list[list[int]] = [[] for _ in range(len(indices))]
    for voxel_id, index in enumerate(indices):
        x, y, z = (int(value) for value in index)
        for dx, dy, dz in offsets:
            neighbor_id = lookup.get((x + dx, y + dy, z + dz))
            if neighbor_id is not None:
                adjacency[voxel_id].append(neighbor_id)
        adjacency[voxel_id].sort()
    return indices, adjacency


def _connected_groups(voxel_ids: set[int], adjacency: list[list[int]]) -> list[list[int]]:
    remaining = set(voxel_ids)
    groups: list[list[int]] = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        stack = [start]
        group: list[int] = []
        while stack:
            current = stack.pop()
            group.append(current)
            for neighbor in adjacency[current]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        groups.append(sorted(group))
    return groups


def _representative_voxel(region_indices: np.ndarray) -> tuple[int, int, int]:
    centroid = np.mean(region_indices, axis=0)
    selected = int(np.argmin(np.linalg.norm(region_indices - centroid, axis=1)))
    return tuple(int(value) for value in region_indices[selected])


def _edge_key(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _trace_branches(
    indices: np.ndarray,
    adjacency: list[list[int]],
    node_regions: list[list[int]],
) -> tuple[list[tuple[int, int, list[int]]], dict[int, int], set[tuple[int, int]]]:
    voxel_to_node: dict[int, int] = {}
    for node_id, region in enumerate(node_regions):
        for voxel_id in region:
            voxel_to_node[voxel_id] = node_id

    visited_edges: set[tuple[int, int]] = set()
    internal_node_edges: set[tuple[int, int]] = set()
    for voxel_id, node_id in voxel_to_node.items():
        for neighbor in adjacency[voxel_id]:
            if voxel_to_node.get(neighbor) == node_id:
                internal_node_edges.add(_edge_key(voxel_id, neighbor))

    traced: list[tuple[int, int, list[int]]] = []
    for node_u, region in enumerate(node_regions):
        for start in region:
            for neighbor in adjacency[start]:
                if voxel_to_node.get(neighbor) == node_u:
                    continue
                first_edge = _edge_key(start, neighbor)
                if first_edge in visited_edges:
                    continue
                path = [start]
                previous = start
                current = neighbor
                visited_edges.add(first_edge)
                while current not in voxel_to_node:
                    path.append(current)
                    candidates = [item for item in adjacency[current] if item != previous]
                    if len(candidates) != 1:
                        coordinate = tuple(int(value) for value in indices[current])
                        raise RuntimeError(
                            "A non-node skeleton voxel does not have exactly two neighbors: "
                            f"index={coordinate}, remaining_neighbors={len(candidates)}"
                        )
                    following = candidates[0]
                    edge = _edge_key(current, following)
                    if edge in visited_edges:
                        raise RuntimeError("Skeleton tracing encountered an already visited path edge")
                    visited_edges.add(edge)
                    previous, current = current, following
                path.append(current)
                traced.append((node_u, voxel_to_node[current], path))

    all_edges = {
        _edge_key(voxel_id, neighbor)
        for voxel_id, neighbors in enumerate(adjacency)
        for neighbor in neighbors
        if voxel_id < neighbor
    }
    unexplained_edges = all_edges - visited_edges - internal_node_edges
    if unexplained_edges:
        raise RuntimeError(
            f"Branch tracing left {len(unexplained_edges)} non-node skeleton edges unexplained"
        )
    return traced, voxel_to_node, internal_node_edges


def _node_type(degree: int, provisional_type: str) -> str:
    if provisional_type == "cycle_anchor":
        return "cycle_anchor"
    if degree == 0:
        return "isolated"
    if degree == 1:
        return "terminal"
    if degree == 2:
        return "connector"
    if degree == 3:
        return "junction"
    return "complex_junction"


def _outward_direction(
    branch: BranchRecord,
    endpoint: str,
    distance_um: float,
) -> tuple[np.ndarray, float]:
    if endpoint == "start":
        points = branch.points_smoothed_lps_um
        radii = branch.coarse_radius_smoothed_um
    else:
        points = branch.points_smoothed_lps_um[::-1]
        radii = branch.coarse_radius_smoothed_um[::-1]
    distances = cumulative_length(points)
    target_index = min(int(np.searchsorted(distances, distance_um)), len(points) - 1)
    if target_index == 0 and len(points) > 1:
        target_index = 1
    vector = points[target_index] - points[0]
    norm = float(np.linalg.norm(vector))
    direction = vector / norm if norm > np.finfo(float).eps else np.zeros(3, dtype=float)
    return direction, float(radii[0])


def _junction_geometry(
    nodes: list[NodeRecord],
    branches: list[BranchRecord],
    direction_distance_um: float,
) -> list[JunctionGeometry]:
    output: list[JunctionGeometry] = []
    for node in nodes:
        if node.graph_degree < 3:
            continue
        incident: list[dict[str, object]] = []
        for branch in branches:
            endpoints: list[str] = []
            if branch.node_u == node.node_id:
                endpoints.append("start")
            if branch.node_v == node.node_id:
                endpoints.append("end")
            for endpoint in endpoints:
                direction, radius = _outward_direction(
                    branch, endpoint, direction_distance_um
                )
                incident.append(
                    {
                        "branch_id": branch.branch_id,
                        "endpoint": endpoint,
                        "outward_direction_lps": direction.tolist(),
                        "coarse_attachment_radius_um": radius,
                        "direction_sampling_distance_um": direction_distance_um,
                    }
                )
        angles: list[dict[str, object]] = []
        for left, right in combinations(incident, 2):
            left_direction = np.asarray(left["outward_direction_lps"], dtype=float)
            right_direction = np.asarray(right["outward_direction_lps"], dtype=float)
            denominator = float(np.linalg.norm(left_direction) * np.linalg.norm(right_direction))
            angle = (
                float(
                    np.degrees(
                        np.arccos(
                            np.clip(np.dot(left_direction, right_direction) / denominator, -1.0, 1.0)
                        )
                    )
                )
                if denominator > np.finfo(float).eps
                else None
            )
            angles.append(
                {
                    "branch_a": left["branch_id"],
                    "branch_a_endpoint": left["endpoint"],
                    "branch_b": right["branch_id"],
                    "branch_b_endpoint": right["endpoint"],
                    "angle_deg": angle,
                }
            )
        output.append(JunctionGeometry(node.node_id, incident, angles))
    return output


def build_hierarchical_graph(
    skeleton: np.ndarray,
    mask: np.ndarray,
    origin_lps_um: tuple[float, float, float],
    spacing_um: tuple[float, float, float],
    config: HierarchicalGraphConfig,
    logger: logging.Logger | None = None,
) -> HierarchicalGraphResult:
    """Aggregate degree-2 chains while retaining raw voxel geometry."""

    logger = logger or logging.getLogger("ulm_3d_vascular")
    config.validate()
    skeleton = np.asarray(skeleton, dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    if skeleton.ndim != 3 or mask.shape != skeleton.shape:
        raise ValueError("Skeleton and mask must be 3-D arrays with the same shape")
    if not np.any(skeleton):
        raise ValueError("A non-empty skeleton is required")
    if np.any(skeleton & ~mask):
        raise ValueError("Skeleton contains voxels outside the vascular mask")

    logger.info(
        "Building voxel adjacency with %d-neighbor connectivity", config.neighbor_connectivity
    )
    indices, adjacency = _voxel_adjacency(skeleton, config.neighbor_connectivity)
    degrees = np.asarray([len(items) for items in adjacency], dtype=int)
    junction_ids = set(np.flatnonzero(degrees >= 3).tolist())
    terminal_ids = np.flatnonzero(degrees <= 1).tolist()
    junction_groups = _connected_groups(junction_ids, adjacency)

    provisional_regions: list[tuple[str, list[int]]] = [
        ("junction", group) for group in junction_groups
    ]
    provisional_regions.extend(("terminal", [voxel_id]) for voxel_id in terminal_ids)
    if not provisional_regions:
        anchor = int(np.lexsort((indices[:, 2], indices[:, 1], indices[:, 0]))[0])
        provisional_regions.append(("cycle_anchor", [anchor]))
        logger.info("The skeleton is a pure degree-2 cycle; inserted one cycle anchor")
    provisional_regions.sort(key=lambda item: tuple(int(value) for value in indices[min(item[1])]))
    node_regions = [item[1] for item in provisional_regions]
    provisional_types = [item[0] for item in provisional_regions]

    traced, voxel_to_node, _ = _trace_branches(indices, adjacency, node_regions)
    radius_field = ndimage.distance_transform_edt(mask, sampling=spacing_um)
    resample_step_um = config.resample_step_um or float(min(spacing_um))

    branches: list[BranchRecord] = []
    for branch_id, (node_u, node_v, path_ids) in enumerate(traced):
        path_indices = indices[np.asarray(path_ids, dtype=int)]
        raw_points = indices_to_lps(path_indices, origin_lps_um, spacing_um)
        raw_length = cumulative_length(raw_points)
        radius_raw = radius_field[tuple(path_indices.T)].astype(float, copy=False)
        (
            smoothed_points,
            smoothed_length,
            radius_smoothed,
            directions,
            curvature,
            max_deviation,
        ) = smooth_and_measure(
            raw_points,
            radius_raw,
            resample_step_um=resample_step_um,
            smoothing_enabled=config.smoothing_enabled,
            smoothing_window_points=config.smoothing_window_points,
        )
        straight_distance = float(np.linalg.norm(raw_points[-1] - raw_points[0]))
        branches.append(
            BranchRecord(
                branch_id=branch_id,
                node_u=node_u,
                node_v=node_v,
                voxel_indices_xyz=path_indices,
                points_raw_lps_um=raw_points,
                arc_length_raw_um=raw_length,
                coarse_radius_raw_um=radius_raw,
                points_smoothed_lps_um=smoothed_points,
                arc_length_smoothed_um=smoothed_length,
                coarse_radius_smoothed_um=radius_smoothed,
                local_direction_smoothed=directions,
                curvature_smoothed_per_um=curvature,
                straight_distance_um=straight_distance,
                tortuosity_raw=tortuosity(float(raw_length[-1]), straight_distance),
                tortuosity_smoothed=tortuosity(
                    float(smoothed_length[-1]), straight_distance
                ),
                smoothing_max_deviation_um=max_deviation,
            )
        )

    incident_ids: dict[int, list[int]] = {node_id: [] for node_id in range(len(node_regions))}
    degree_counts: Counter[int] = Counter()
    for branch in branches:
        incident_ids[branch.node_u].append(branch.branch_id)
        incident_ids[branch.node_v].append(branch.branch_id)
        if branch.node_u == branch.node_v:
            degree_counts[branch.node_u] += 2
        else:
            degree_counts[branch.node_u] += 1
            degree_counts[branch.node_v] += 1

    nodes: list[NodeRecord] = []
    for node_id, region_ids in enumerate(node_regions):
        region_indices = indices[np.asarray(region_ids, dtype=int)]
        representative_index = _representative_voxel(region_indices)
        representative_lps = tuple(
            float(value)
            for value in indices_to_lps(
                np.asarray([representative_index]), origin_lps_um, spacing_um
            )[0]
        )
        graph_degree = int(degree_counts[node_id])
        nodes.append(
            NodeRecord(
                node_id=node_id,
                node_type=_node_type(graph_degree, provisional_types[node_id]),
                graph_degree=graph_degree,
                voxel_indices_xyz=region_indices,
                representative_index_xyz=representative_index,
                representative_lps_um=representative_lps,
                incident_branch_ids=sorted(set(incident_ids[node_id])),
            )
        )

    junction_graph = build_junction_graph(nodes, branches)
    cycles, cycle_rank = identify_cycles(junction_graph, nodes, branches)
    junction_graph = build_junction_graph(nodes, branches)
    branch_as_node_graph = build_branch_as_node_graph(nodes, branches)
    junctions = _junction_geometry(
        nodes, branches, config.junction_direction_distance_um
    )

    reconstructed = np.zeros_like(skeleton, dtype=bool)
    for node in nodes:
        reconstructed[tuple(node.voxel_indices_xyz.T)] = True
    regular_usage: Counter[tuple[int, int, int]] = Counter()
    node_voxels = {
        tuple(int(value) for value in indices[voxel_id]) for voxel_id in voxel_to_node
    }
    for branch in branches:
        reconstructed[tuple(branch.voxel_indices_xyz.T)] = True
        for index in branch.voxel_indices_xyz:
            key = tuple(int(value) for value in index)
            if key not in node_voxels:
                regular_usage[key] += 1
    duplicate_interior = int(sum(count > 1 for count in regular_usage.values()))
    missing = int(np.count_nonzero(skeleton & ~reconstructed))
    extra = int(np.count_nonzero(reconstructed & ~skeleton))
    structure = _ndimage_structure(config.neighbor_connectivity)
    skeleton_components = int(ndimage.label(skeleton, structure=structure)[1])
    graph_components = (
        int(nx.number_connected_components(junction_graph))
        if junction_graph.number_of_nodes()
        else 0
    )
    logger.info(
        "Hierarchical graph: %d nodes, %d branches, %d cycle-rank; "
        "round-trip missing=%d, extra=%d",
        len(nodes),
        len(branches),
        cycle_rank,
        missing,
        extra,
    )
    return HierarchicalGraphResult(
        nodes=nodes,
        branches=branches,
        junctions=junctions,
        cycles=cycles,
        junction_graph=junction_graph,
        branch_as_node_graph=branch_as_node_graph,
        source_skeleton=skeleton,
        reconstructed_skeleton=reconstructed,
        skeleton_component_count=skeleton_components,
        graph_component_count=graph_components,
        skeleton_voxel_count=int(np.count_nonzero(skeleton)),
        represented_voxel_count=int(np.count_nonzero(reconstructed & skeleton)),
        missing_voxel_count=missing,
        extra_voxel_count=extra,
        duplicate_interior_voxel_count=duplicate_interior,
        cycle_rank=cycle_rank,
        coordinate_system="LPS",
        physical_unit="micrometer",
        origin_lps_um=tuple(float(value) for value in origin_lps_um),
        spacing_um=tuple(float(value) for value in spacing_um),
        radius_source="voxel distance transform; coarse navigation estimate",
    )
