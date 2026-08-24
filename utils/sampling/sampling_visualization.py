"""Acceptance figures for candidate coverage, selection, and distributions."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .radius_features import aggregate_radius_distribution
from .sampling_types import GlobalVascularModel, ROIRecord, SamplingExperiment


def _box_edges(minimum: np.ndarray, maximum: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    corners = np.asarray(
        [
            (x, y, z)
            for x in (minimum[0], maximum[0])
            for y in (minimum[1], maximum[1])
            for z in (minimum[2], maximum[2])
        ],
        dtype=float,
    )
    return [
        (corners[first], corners[second])
        for first in range(8)
        for second in range(first + 1, 8)
        if np.count_nonzero(corners[first] != corners[second]) == 1
    ]


def _style_3d(axis: object) -> None:
    axis.set_xlabel("X (um)")
    axis.set_ylabel("Y (um)")
    axis.set_zlabel("Z (um)")
    axis.set_box_aspect((1, 1, 1))


def global_roi_overview(
    models: list[GlobalVascularModel],
    rois: list[ROIRecord],
    path: Path,
    *,
    selected_only: bool,
) -> Path:
    figure = plt.figure(figsize=(12, 9))
    axis = figure.add_subplot(111, projection="3d")
    for model in models:
        for edge in model.edges:
            points = np.asarray((edge.upstream_position_um, edge.downstream_position_um))
            axis.plot(points[:, 0], points[:, 1], points[:, 2], color="#78DCE8", alpha=0.22, linewidth=0.45)
    palette = plt.get_cmap("tab20")
    for index, roi in enumerate(rois):
        color = palette(roi.cluster_id % 20) if selected_only else (1.0, 0.35, 0.2, 0.38)
        minimum = np.asarray(roi.bbox_min_um)
        maximum = np.asarray(roi.bbox_max_um)
        for start, end in _box_edges(minimum, maximum):
            axis.plot(
                (start[0], end[0]),
                (start[1], end[1]),
                (start[2], end[2]),
                color=color,
                linewidth=1.8 if selected_only else 0.75,
                alpha=0.9 if selected_only else 0.35,
            )
        center = np.asarray(roi.bbox_center_um)
        axis.scatter(*center, color=color, s=28 if selected_only else 8)
        if selected_only:
            axis.text(*center, f"C{roi.cluster_id}", color=color, fontsize=7)
    axis.set_title(
        "Selected representative real ROIs" if selected_only else "All connected candidate ROIs"
    )
    _style_3d(axis)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def radius_distribution_figure(
    candidates: list[ROIRecord],
    selected: list[ROIRecord],
    path: Path,
) -> Path:
    candidate_values, candidate_weights = aggregate_radius_distribution(candidates)
    selected_values, selected_weights = aggregate_radius_distribution(selected)
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    if len(candidate_values) and len(selected_values):
        limits = (min(candidate_values.min(), selected_values.min()), max(candidate_values.max(), selected_values.max()))
        bins = np.linspace(limits[0], limits[1], 28)
        axes[0].hist(candidate_values, bins=bins, weights=candidate_weights / candidate_weights.sum(), alpha=0.55, label="candidates")
        axes[0].hist(selected_values, bins=bins, weights=selected_weights / selected_weights.sum(), alpha=0.55, label="selected")
        for values, weights, label in (
            (candidate_values, candidate_weights, "candidates"),
            (selected_values, selected_weights, "selected"),
        ):
            order = np.argsort(values)
            axes[1].step(values[order], np.cumsum(weights[order]) / np.sum(weights), where="post", label=label)
    axes[0].set_title("Arc-length weighted radius histogram")
    axes[1].set_title("Arc-length weighted radius CDF")
    for axis in axes:
        axis.set_xlabel("Radius (um)")
        axis.grid(alpha=0.2)
        axis.legend()
    axes[0].set_ylabel("Weighted probability")
    axes[1].set_ylabel("Weighted CDF")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def scalar_distribution_figure(
    candidates: list[ROIRecord],
    selected: list[ROIRecord],
    feature_name: str,
    title: str,
    path: Path,
) -> Path:
    candidate = np.asarray([roi.structural_features[feature_name] for roi in candidates], dtype=float)
    chosen = np.asarray([roi.structural_features[feature_name] for roi in selected], dtype=float)
    figure, axis = plt.subplots(figsize=(7.5, 4.5))
    if len(candidate):
        bins = min(18, max(5, int(np.sqrt(len(candidate)) * 2)))
        axis.hist(candidate, bins=bins, density=True, alpha=0.55, label="candidates")
        if len(chosen):
            axis.hist(chosen, bins=bins, density=True, alpha=0.55, label="selected")
    axis.set_title(title)
    axis.set_xlabel(feature_name)
    axis.set_ylabel("Density")
    axis.grid(alpha=0.2)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def cluster_size_figure(experiment: SamplingExperiment, path: Path) -> Path:
    figure, axis = plt.subplots(figsize=(8, 4.5))
    x = np.arange(experiment.clustering.n_clusters)
    axis.bar(x, experiment.clustering.cluster_sizes, color="#4C78A8")
    axis.set_xlabel("Cluster ID")
    axis.set_ylabel("Candidate ROI count")
    axis.set_title("Cluster size distribution")
    axis.set_xticks(x)
    axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def pca_feature_space_figure(
    experiment: SamplingExperiment,
    path: Path,
) -> Path:
    values = experiment.clustering.scaled_features
    centered = values - np.mean(values, axis=0)
    _, singular, right = np.linalg.svd(centered, full_matrices=False)
    components = centered @ right[: min(2, len(right))].T
    if components.shape[1] == 1:
        components = np.column_stack((components[:, 0], np.zeros(len(components))))
    total_variance = float(np.sum(singular ** 2))
    explained = (singular[:2] ** 2 / total_variance) if total_variance > 0 else np.zeros(2)
    figure, axis = plt.subplots(figsize=(8, 6))
    scatter = axis.scatter(
        components[:, 0],
        components[:, 1],
        c=experiment.clustering.assignments,
        cmap="tab20",
        s=28,
        alpha=0.7,
        label="candidate ROI",
    )
    selected = np.asarray(experiment.selection.selected_indices, dtype=int)
    if len(selected):
        axis.scatter(
            components[selected, 0],
            components[selected, 1],
            facecolors="none",
            edgecolors="black",
            linewidths=1.8,
            s=110,
            label="selected real ROI",
        )
    axis.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    axis.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    axis.set_title(f"Scaled feature space: {experiment.feature_mode}")
    axis.legend()
    figure.colorbar(scatter, ax=axis, label="Cluster ID")
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def silhouette_scan_figure(scan: list[dict[str, float | int | None]], path: Path) -> Path:
    figure, first = plt.subplots(figsize=(8, 4.5))
    counts = np.asarray([row["n_clusters"] for row in scan], dtype=int)
    silhouette = np.asarray([row["silhouette_score"] if row["silhouette_score"] is not None else np.nan for row in scan], dtype=float)
    inertia = np.asarray([row["inertia"] for row in scan], dtype=float)
    first.plot(counts, silhouette, marker="o", color="#E45756", label="silhouette")
    first.set_xlabel("K")
    first.set_ylabel("Silhouette score", color="#E45756")
    second = first.twinx()
    second.plot(counts, inertia, marker="s", color="#4C78A8", label="inertia")
    second.set_ylabel("Inertia", color="#4C78A8")
    first.set_title("Exploratory K scan (no automatic K choice)")
    first.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def roi_preview(roi: ROIRecord, path: Path) -> Path:
    figure = plt.figure(figsize=(7, 6))
    axis = figure.add_subplot(111, projection="3d")
    for edge in roi.local_edge_points_um:
        axis.plot(edge[:, 0], edge[:, 1], edge[:, 2], color="#45D6E8", linewidth=2.0)
    if roi.true_terminal_local_ids:
        points = roi.local_node_positions_um[list(roi.true_terminal_local_ids)]
        axis.scatter(points[:, 0], points[:, 1], points[:, 2], color="#F87171", s=34, label="TRUE_TERMINAL")
    if roi.cut_ports:
        points = np.asarray([port.intersection_position_um for port in roi.cut_ports])
        axis.scatter(points[:, 0], points[:, 1], points[:, 2], color="#F59E0B", marker="^", s=48, label="CUT_PORT")
    for start, end in _box_edges(np.asarray(roi.bbox_min_um), np.asarray(roi.bbox_max_um)):
        axis.plot((start[0], end[0]), (start[1], end[1]), (start[2], end[2]), color="#FF5C5C", linewidth=1.0)
    axis.set_title(f"{roi.roi_id}\ncluster={roi.cluster_id}, branches={roi.branch_count}, cuts={roi.cut_port_count}")
    _style_3d(axis)
    if roi.true_terminal_local_ids or roi.cut_ports:
        axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def feature_mode_comparison_figure(
    experiments: dict[str, SamplingExperiment],
    path: Path,
) -> Path:
    metric_names = ("branch_count", "bifurcation_count", "total_vessel_length_um", "cycle_rank")
    modes = list(experiments)
    matrix = np.asarray(
        [
            [
                float(experiments[mode].validation["structure"][name]["ks"] or 0.0)
                for name in metric_names
            ]
            for mode in modes
        ]
    )
    figure, axis = plt.subplots(figsize=(10, 5))
    x = np.arange(len(metric_names))
    width = 0.8 / max(1, len(modes))
    for index, mode in enumerate(modes):
        axis.bar(x + (index - (len(modes) - 1) / 2) * width, matrix[index], width, label=mode)
    axis.set_xticks(x, ("branches", "bifurcations", "length", "cycle rank"))
    axis.set_ylabel("Candidate vs selected KS statistic")
    axis.set_title("Radius-only vs radius-plus-structure sampling")
    axis.legend()
    axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path

