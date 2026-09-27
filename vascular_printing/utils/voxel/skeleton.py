"""Extract a navigation-quality coarse skeleton from a voxel mask."""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy import ndimage
from skimage.morphology import skeletonize


@dataclass(slots=True)
class SkeletonResult:
    skeleton: np.ndarray
    skeleton_voxel_count: int
    connected_component_count: int
    voxels_outside_mask: int
    approximate_length_um: float

    def report(self) -> dict[str, Any]:
        payload = asdict(self)
        payload.pop("skeleton")
        return payload


def extract_coarse_skeleton(
    mask: np.ndarray,
    spacing_um: tuple[float, float, float],
    logger: logging.Logger | None = None,
) -> SkeletonResult:
    logger = logger or logging.getLogger("ulm_3d_vascular")
    if mask.ndim != 3 or not np.any(mask):
        raise ValueError("A non-empty 3-D mask is required for skeletonization")
    logger.info("Extracting coarse 3-D skeleton with the Lee method")
    skeleton = np.asarray(skeletonize(mask, method="lee"), dtype=bool)
    voxel_count = int(np.count_nonzero(skeleton))
    if voxel_count == 0:
        raise ValueError("Skeletonization produced an empty result")
    outside = int(np.count_nonzero(skeleton & ~mask))
    _, component_count = ndimage.label(
        skeleton, structure=ndimage.generate_binary_structure(3, 3)
    )
    return SkeletonResult(
        skeleton=skeleton,
        skeleton_voxel_count=voxel_count,
        connected_component_count=int(component_count),
        voxels_outside_mask=outside,
        approximate_length_um=voxel_count * float(np.mean(spacing_um)),
    )

