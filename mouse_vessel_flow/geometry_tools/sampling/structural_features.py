"""Topology and material-size features for connected ROI graphs."""

from __future__ import annotations

import networkx as nx
import numpy as np

from .sampling_types import ROIRecord


def _graph(roi: ROIRecord) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(int(value) for value in roi.local_node_ids)
    graph.add_edges_from((int(a), int(b)) for a, b in roi.local_edges)
    return graph


def _branch_paths(graph: nx.Graph) -> list[list[int]]:
    critical = {int(node) for node in graph if graph.degree(node) != 2}
    visited: set[frozenset[int]] = set()
    paths: list[list[int]] = []
    for start in sorted(critical):
        for neighbor in sorted(graph.neighbors(start)):
            marker = frozenset((start, int(neighbor)))
            if marker in visited:
                continue
            path = [start, int(neighbor)]
            visited.add(marker)
            previous, current = start, int(neighbor)
            while current not in critical:
                following = [int(value) for value in graph.neighbors(current) if int(value) != previous]
                if not following:
                    break
                next_node = following[0]
                visited.add(frozenset((current, next_node)))
                path.append(next_node)
                previous, current = current, next_node
            paths.append(path)
    for first, second in sorted((int(a), int(b)) for a, b in graph.edges):
        marker = frozenset((first, second))
        if marker in visited:
            continue
        path = [first, second]
        visited.add(marker)
        previous, current = first, second
        while True:
            following = [
                int(value)
                for value in graph.neighbors(current)
                if int(value) != previous and frozenset((current, int(value))) not in visited
            ]
            if not following:
                break
            next_node = following[0]
            visited.add(frozenset((current, next_node)))
            path.append(next_node)
            previous, current = current, next_node
        paths.append(path)
    return paths


def compute_structural_features(roi: ROIRecord) -> dict[str, float]:
    graph = _graph(roi)
    component_count = nx.number_connected_components(graph) if graph else 0
    cycle_rank = graph.number_of_edges() - graph.number_of_nodes() + component_count
    branch_paths = _branch_paths(graph)
    total_length = float(
        np.sum(np.linalg.norm(roi.local_edge_points_um[:, 1] - roi.local_edge_points_um[:, 0], axis=1))
    )
    volume = float(np.prod(np.asarray(roi.bbox_size_um, dtype=float)))
    tortuosities: list[float] = []
    for path in branch_paths:
        points = roi.local_node_positions_um[path]
        length = float(np.sum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
        chord = float(np.linalg.norm(points[-1] - points[0])) if len(points) >= 2 else 0.0
        if chord > 1.0e-12:
            tortuosities.append(length / chord)
    branch_count = len(branch_paths)
    bifurcation_count = sum(graph.degree(node) >= 3 for node in graph)
    return {
        "branch_count": float(branch_count),
        "bifurcation_count": float(bifurcation_count),
        "total_vessel_length_um": total_length,
        "vessel_length_density_um_per_um3": total_length / volume if volume > 0 else float("nan"),
        "branch_density_per_um3": branch_count / volume if volume > 0 else float("nan"),
        "bifurcation_density_per_um3": bifurcation_count / volume if volume > 0 else float("nan"),
        "cycle_rank": float(cycle_rank),
        "mean_branch_tortuosity": float(np.mean(tortuosities)) if tortuosities else 1.0,
    }

