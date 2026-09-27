"""Independent distribution, topology, and traceability validation."""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
from scipy.stats import wasserstein_distance

from .overlap_control import overlap_statistics
from .radius_features import aggregate_radius_distribution, weighted_quantile
from .sampling_types import ClusteringResult, ROIRecord, SelectionResult


def _empirical_ks(first: np.ndarray, second: np.ndarray) -> float | None:
    if not len(first) or not len(second):
        return None
    support = np.unique(np.concatenate((first, second)))
    first_sorted = np.sort(first)
    second_sorted = np.sort(second)
    first_cdf = np.searchsorted(first_sorted, support, side="right") / len(first_sorted)
    second_cdf = np.searchsorted(second_sorted, support, side="right") / len(second_sorted)
    return float(np.max(np.abs(first_cdf - second_cdf)))


def _weighted_ks(
    first: np.ndarray,
    first_weights: np.ndarray,
    second: np.ndarray,
    second_weights: np.ndarray,
) -> float | None:
    if not len(first) or not len(second):
        return None
    support = np.unique(np.concatenate((first, second)))

    def cdf(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
        order = np.argsort(values, kind="stable")
        sorted_values = values[order]
        cumulative = np.cumsum(weights[order])
        indices = np.searchsorted(sorted_values, support, side="right") - 1
        result = np.zeros(len(support), dtype=float)
        valid = indices >= 0
        result[valid] = cumulative[indices[valid]] / cumulative[-1]
        return result

    return float(np.max(np.abs(cdf(first, first_weights) - cdf(second, second_weights))))


def _quantile_error(first: np.ndarray, second: np.ndarray) -> float | None:
    if not len(first) or not len(second):
        return None
    quantiles = np.asarray((0.10, 0.25, 0.50, 0.75, 0.90))
    reference = np.quantile(first, quantiles)
    selected = np.quantile(second, quantiles)
    return float(np.mean(np.abs(reference - selected) / (np.abs(reference) + 1.0e-12)))


def _scalar_metrics(first: np.ndarray, second: np.ndarray) -> dict[str, float | None]:
    if not len(first) or not len(second):
        return {"wasserstein": None, "ks": None, "quantile_error": None}
    return {
        "wasserstein": float(wasserstein_distance(first, second)),
        "ks": _empirical_ks(first, second),
        "quantile_error": _quantile_error(first, second),
    }


def validate_roi_integrity(rois: list[ROIRecord]) -> dict[str, Any]:
    connected = 0
    traceable = 0
    boundary_valid = 0
    for roi in rois:
        graph = nx.Graph()
        graph.add_nodes_from(int(value) for value in roi.local_node_ids)
        graph.add_edges_from((int(a), int(b)) for a, b in roi.local_edges)
        connected += int(bool(graph) and nx.is_connected(graph))
        traceable += int(
            len(roi.global_edge_ids) == len(roi.local_edge_global_ids)
            and set(map(int, roi.local_edge_global_ids)) == set(roi.global_edge_ids)
        )
        lower = np.asarray(roi.bbox_min_um)
        upper = np.asarray(roi.bbox_max_um)
        boundary_valid += int(
            all(
                np.any(np.isclose(port.intersection_position_um, lower, atol=1.0e-6))
                or np.any(np.isclose(port.intersection_position_um, upper, atol=1.0e-6))
                for port in roi.cut_ports
            )
        )
    count = len(rois)
    return {
        "candidate_count": count,
        "connected_count": connected,
        "globally_traceable_count": traceable,
        "boundary_semantics_valid_count": boundary_valid,
        "all_connected": connected == count,
        "all_globally_traceable": traceable == count,
        "all_boundary_semantics_valid": boundary_valid == count,
    }


def validate_sampling(
    rois: list[ROIRecord],
    clustering: ClusteringResult,
    selection: SelectionResult,
) -> dict[str, object]:
    selected = [rois[index] for index in selection.selected_indices]
    candidate_radius, candidate_weight = aggregate_radius_distribution(rois)
    selected_radius, selected_weight = aggregate_radius_distribution(selected)
    radius_quantiles = np.asarray((0.10, 0.25, 0.50, 0.75, 0.90))
    candidate_q = weighted_quantile(candidate_radius, candidate_weight, radius_quantiles)
    selected_q = weighted_quantile(selected_radius, selected_weight, radius_quantiles)
    radius_quantile_error = (
        float(np.mean(np.abs(candidate_q - selected_q) / (np.abs(candidate_q) + 1.0e-12)))
        if len(candidate_radius) and len(selected_radius)
        else None
    )
    radius_metrics = {
        "wasserstein": (
            float(
                wasserstein_distance(
                    candidate_radius,
                    selected_radius,
                    u_weights=candidate_weight,
                    v_weights=selected_weight,
                )
            )
            if len(candidate_radius) and len(selected_radius)
            else None
        ),
        "ks": _weighted_ks(candidate_radius, candidate_weight, selected_radius, selected_weight),
        "quantile_error": radius_quantile_error,
    }
    feature_metrics: dict[str, dict[str, float | None]] = {}
    for name in (
        "branch_count",
        "bifurcation_count",
        "total_vessel_length_um",
        "cycle_rank",
    ):
        candidate_values = np.asarray([roi.structural_features[name] for roi in rois], dtype=float)
        selected_values = np.asarray([roi.structural_features[name] for roi in selected], dtype=float)
        feature_metrics[name] = _scalar_metrics(candidate_values, selected_values)
    represented = set(int(clustering.assignments[index]) for index in selection.selected_indices)
    return {
        "radius": radius_metrics,
        "structure": feature_metrics,
        "cluster_coverage": len(represented) / clustering.n_clusters if clustering.n_clusters else 0.0,
        "represented_cluster_count": len(represented),
        "cluster_count": clustering.n_clusters,
        "selected_count": len(selected),
        "representatives_are_real_candidates": all(0 <= index < len(rois) for index in selection.selected_indices),
        "spatial_redundancy": overlap_statistics(selected),
        "roi_integrity": validate_roi_integrity(rois),
    }
