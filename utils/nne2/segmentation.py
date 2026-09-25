"""Fast image-stack normalization and vessel segmentation for NNE2."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import ndimage
from skimage.filters import threshold_otsu
from .config import NNE2Config
from .components import clean_components


@dataclass(slots=True)
class SegmentationResult:
    normalized_zyx: np.ndarray
    vessel_score_zyx: np.ndarray
    candidate_mask_zyx: np.ndarray
    mask_zyx: np.ndarray
    removed_mask_zyx: np.ndarray
    threshold: float
    component_count_before: int
    component_count_after: int
    removed_voxel_count: int
    component_decisions: list[dict[str, Any]]

    def report(self) -> dict[str, Any]:
        total = int(self.mask_zyx.size)
        foreground = int(np.count_nonzero(self.mask_zyx))
        return {
            "shape_zyx": self.mask_zyx.shape,
            "threshold": self.threshold,
            "foreground_voxel_count": foreground,
            "foreground_fraction": foreground / total if total else 0.0,
            "component_count_before_cleanup": self.component_count_before,
            "component_count_after_cleanup": self.component_count_after,
            "removed_voxel_count": self.removed_voxel_count,
            "removed_fraction_of_candidate": (
                self.removed_voxel_count / max(1, int(np.count_nonzero(self.candidate_mask_zyx)))
            ),
            "component_decisions": self.component_decisions,
            "method": (
                "per-plane robust normalization, background suppression, global Otsu/quantile "
                "threshold, 3-D closing and small-component removal"
            ),
        }


def _normalize_per_plane(volume_zyx: np.ndarray) -> np.ndarray:
    volume = np.asarray(volume_zyx, dtype=np.float32)
    low = np.percentile(volume, 2.0, axis=(1, 2)).astype(np.float32)
    high = np.percentile(volume, 99.7, axis=(1, 2)).astype(np.float32)
    scale = np.maximum(high - low, 1.0)
    normalized = (volume - low[:, None, None]) / scale[:, None, None]
    return np.clip(normalized, 0.0, 1.0, out=normalized)


def segment_vessels(
    volume_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    config: NNE2Config,
    logger: logging.Logger | None = None,
) -> SegmentationResult:
    logger = logger or logging.getLogger("ulm_3d_vascular")
    if volume_zyx.ndim != 3 or not volume_zyx.size:
        raise ValueError("A non-empty 3-D NNE2 stack is required")
    normalized = _normalize_per_plane(volume_zyx)
    sx, sy, sz = spacing_xyz_um
    smooth_sigma = (
        max(0.25, config.gaussian_sigma_um / sz),
        max(0.25, config.gaussian_sigma_um / sy),
        max(0.25, config.gaussian_sigma_um / sx),
    )
    smooth = ndimage.gaussian_filter(normalized, sigma=smooth_sigma, mode="nearest")
    background_sigma = (
        0.0,
        max(1.0, config.background_sigma_um / sy),
        max(1.0, config.background_sigma_um / sx),
    )
    background = ndimage.gaussian_filter(smooth, sigma=background_sigma, mode="nearest")
    local_contrast = np.maximum(smooth - background, 0.0)
    score = np.asarray(0.65 * smooth + 0.35 * local_contrast, dtype=np.float32)

    sample = score[:: max(1, score.shape[0] // 64), ::2, ::2].reshape(-1)
    sample = sample[np.isfinite(sample)]
    if not len(sample):
        raise ValueError("Vessel score contains no finite values")
    quantile_threshold = float(np.quantile(sample, config.foreground_quantile))
    positive = sample[sample > 0]
    otsu_threshold = float(threshold_otsu(positive)) if len(positive) > 1 else quantile_threshold
    threshold = max(quantile_threshold, otsu_threshold)
    raw_mask = score >= threshold
    closing_structure = ndimage.generate_binary_structure(3, 1)
    component_structure = ndimage.generate_binary_structure(3, 3)
    before_components = int(ndimage.label(raw_mask, structure=component_structure)[1])
    candidate_mask = raw_mask
    if config.closing_iterations:
        candidate_mask = ndimage.binary_closing(
            candidate_mask, structure=closing_structure, iterations=config.closing_iterations
        )
    mask, removed_mask, component_decisions = clean_components(
        np.asarray(candidate_mask, dtype=bool),
        spacing_xyz_um,
        min_component_voxels=config.min_component_voxels,
    )
    after_components = int(ndimage.label(mask, structure=component_structure)[1])
    removed = int(np.count_nonzero(removed_mask))
    foreground_fraction = float(np.mean(mask))
    logger.info(
        "Segmented vessels: threshold=%.4f, foreground=%.3f%%, components %d -> %d",
        threshold,
        100.0 * foreground_fraction,
        before_components,
        after_components,
    )
    if foreground_fraction <= 0 or foreground_fraction > 0.35:
        raise ValueError(
            f"Implausible segmented foreground fraction: {foreground_fraction:.4f}"
        )
    return SegmentationResult(
        normalized_zyx=normalized,
        vessel_score_zyx=score,
        candidate_mask_zyx=np.asarray(candidate_mask, dtype=bool),
        mask_zyx=np.asarray(mask, dtype=bool),
        removed_mask_zyx=np.asarray(removed_mask, dtype=bool),
        threshold=threshold,
        component_count_before=before_components,
        component_count_after=after_components,
        removed_voxel_count=removed,
        component_decisions=component_decisions,
    )
