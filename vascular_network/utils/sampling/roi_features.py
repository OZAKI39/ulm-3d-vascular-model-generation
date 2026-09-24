"""Unified, configurable sampling descriptors."""

from __future__ import annotations

import numpy as np

from .radius_features import compute_radius_features
from .sampling_config import SamplingConfig
from .sampling_types import ROIRecord
from .structural_features import compute_structural_features


RADIUS_FEATURE_NAMES = ("r10", "r25", "r50", "r75", "r90")
STRUCTURE_COUNT_NAMES = (
    "branch_count",
    "bifurcation_count",
    "total_vessel_length_um",
    "cycle_rank",
)
STRUCTURE_DENSITY_NAMES = (
    "branch_density_per_um3",
    "bifurcation_density_per_um3",
    "vessel_length_density_um_per_um3",
    "cycle_rank",
)


def populate_roi_features(roi: ROIRecord) -> None:
    roi.radius_features = compute_radius_features(roi)
    roi.structural_features = compute_structural_features(roi)


def feature_names_for_mode(
    mode: str,
    *,
    variable_roi_size: bool,
) -> tuple[str, ...]:
    if mode == "radius_only":
        return RADIUS_FEATURE_NAMES
    structure = STRUCTURE_DENSITY_NAMES if variable_roi_size else STRUCTURE_COUNT_NAMES
    if mode == "radius_plus_structure":
        return RADIUS_FEATURE_NAMES + structure
    if mode == "extended_morphology":
        return RADIUS_FEATURE_NAMES + structure + ("mean_branch_tortuosity",)
    raise ValueError(f"Unsupported feature mode: {mode}")


def build_feature_matrix(
    rois: list[ROIRecord],
    config: SamplingConfig,
    *,
    feature_mode: str | None = None,
) -> tuple[np.ndarray, tuple[str, ...], tuple[str, ...]]:
    if not rois:
        raise ValueError("Cannot build features for an empty ROI collection")
    sizes = np.asarray([roi.bbox_size_um for roi in rois], dtype=float)
    variable_size = bool(np.any(np.ptp(sizes, axis=0) > 1.0e-8))
    mode = feature_mode or config.feature_mode
    names = feature_names_for_mode(mode, variable_roi_size=variable_size)
    matrix = np.asarray(
        [
            [
                roi.radius_features.get(name, roi.structural_features.get(name, float("nan")))
                for name in names
            ]
            for roi in rois
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(matrix)):
        invalid = np.argwhere(~np.isfinite(matrix))
        raise ValueError(f"Non-finite ROI features at indices {invalid[:10].tolist()}")
    groups = tuple("radius" if name in RADIUS_FEATURE_NAMES else "structure" for name in names)
    for roi, vector in zip(rois, matrix):
        if mode == config.feature_mode:
            roi.feature_vector = vector.copy()
    return matrix, names, groups
