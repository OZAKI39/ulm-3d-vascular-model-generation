"""Source-space ROI overlap and redundancy calculations."""

from __future__ import annotations

import numpy as np

from .sampling_types import ROIRecord


def bbox_overlap_fraction(first: ROIRecord, second: ROIRecord) -> float:
    first_min = np.asarray(first.bbox_min_um, dtype=float)
    first_max = np.asarray(first.bbox_max_um, dtype=float)
    second_min = np.asarray(second.bbox_min_um, dtype=float)
    second_max = np.asarray(second.bbox_max_um, dtype=float)
    intersection = np.maximum(0.0, np.minimum(first_max, second_max) - np.maximum(first_min, second_min))
    intersection_volume = float(np.prod(intersection))
    first_volume = float(np.prod(first_max - first_min))
    second_volume = float(np.prod(second_max - second_min))
    denominator = min(first_volume, second_volume)
    return intersection_volume / denominator if denominator > 0 else 0.0


def source_space_distance(first: ROIRecord, second: ROIRecord) -> float:
    if first.source_model_id != second.source_model_id:
        return float("inf")
    return float(
        np.linalg.norm(np.asarray(first.bbox_center_um) - np.asarray(second.bbox_center_um))
    )


def overlap_statistics(rois: list[ROIRecord]) -> dict[str, float | list[float]]:
    overlaps: list[float] = []
    nearest: list[float] = []
    for index, first in enumerate(rois):
        local_distances: list[float] = []
        for second_index, second in enumerate(rois):
            if index == second_index:
                continue
            if second_index > index:
                overlaps.append(bbox_overlap_fraction(first, second))
            distance = source_space_distance(first, second)
            if np.isfinite(distance):
                local_distances.append(distance)
        if local_distances:
            nearest.append(min(local_distances))
    return {
        "mean_spatial_overlap": float(np.mean(overlaps)) if overlaps else 0.0,
        "median_spatial_overlap": float(np.median(overlaps)) if overlaps else 0.0,
        "maximum_spatial_overlap": float(np.max(overlaps)) if overlaps else 0.0,
        "nearest_selected_roi_distance_um": nearest,
        "mean_nearest_selected_roi_distance_um": float(np.mean(nearest)) if nearest else 0.0,
    }

