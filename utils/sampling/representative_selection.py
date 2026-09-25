"""Select actual candidate ROIs while enforcing spatial redundancy limits."""

from __future__ import annotations

import numpy as np

from .overlap_control import bbox_overlap_fraction, source_space_distance
from .sampling_config import SamplingConfig
from .sampling_types import ClusteringResult, ROIRecord, SelectionResult


def _cluster_quotas(clustering: ClusteringResult, config: SamplingConfig) -> np.ndarray:
    sizes = np.asarray(clustering.cluster_sizes, dtype=int)
    quotas = np.zeros(clustering.n_clusters, dtype=int)
    nonempty = np.flatnonzero(sizes > 0)
    target = min(
        config.target_selected_count,
        len(clustering.assignments),
        config.representatives_per_cluster * len(nonempty),
    )
    if config.selection_mode == "coverage_balanced":
        base, remainder = divmod(target, len(nonempty))
        quotas[nonempty] = base
        for cluster_id in nonempty[:remainder]:
            quotas[cluster_id] += 1
    else:
        exact = target * sizes / max(1, np.sum(sizes))
        quotas = np.floor(exact).astype(int)
        remaining = target - int(np.sum(quotas))
        order = np.argsort(-(exact - quotas), kind="stable")
        for cluster_id in order[:remaining]:
            if sizes[cluster_id] > 0:
                quotas[cluster_id] += 1
    quotas = np.minimum(quotas, sizes)
    quotas = np.minimum(quotas, config.representatives_per_cluster)
    return quotas


def select_representatives(
    rois: list[ROIRecord],
    clustering: ClusteringResult,
    config: SamplingConfig,
) -> SelectionResult:
    quotas = _cluster_quotas(clustering, config)
    selected: list[int] = []
    overlap_rejections = 0
    cluster_order = list(range(clustering.n_clusters))
    if config.selection_mode == "coverage_balanced":
        cluster_order.sort(key=lambda cluster_id: (clustering.cluster_sizes[cluster_id], cluster_id))
    for cluster_id in cluster_order:
        quota = int(quotas[cluster_id])
        members = np.flatnonzero(clustering.assignments == cluster_id)
        ordered = sorted(
            members.tolist(),
            key=lambda index: (
                float(clustering.distances_to_center[index]),
                rois[index].roi_id,
            ),
        )
        accepted = 0
        for index in ordered:
            candidate = rois[index]
            excessive_overlap = any(
                bbox_overlap_fraction(candidate, rois[chosen]) > config.max_selected_overlap
                for chosen in selected
            )
            too_close = any(
                source_space_distance(candidate, rois[chosen])
                < config.min_representative_distance_um
                for chosen in selected
            )
            if excessive_overlap or too_close:
                overlap_rejections += 1
                continue
            selected.append(index)
            accepted += 1
            if accepted >= quota:
                break
    scientific_label = (
        "population-representative"
        if config.selection_mode == "distribution_preserving"
        else "benchmark-balanced"
    )
    return SelectionResult(
        config.selection_mode,
        scientific_label,
        tuple(selected),
        overlap_rejections,
        int(np.sum(quotas)),
    )
