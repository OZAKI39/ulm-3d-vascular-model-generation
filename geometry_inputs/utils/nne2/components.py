"""Explainable connected-component decisions for NNE2 vessel masks."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import ndimage


def clean_components(
    candidate_mask_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    *,
    min_component_voxels: int,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Remove only components below the configured size and record every decision."""
    candidate = np.asarray(candidate_mask_zyx, dtype=bool)
    structure = ndimage.generate_binary_structure(3, 3)
    labels, count = ndimage.label(candidate, structure=structure)
    objects = ndimage.find_objects(labels)
    spacing_zyx = np.asarray(spacing_xyz_um[::-1], dtype=float)
    voxel_volume = float(np.prod(spacing_xyz_um))
    cleaned = np.zeros_like(candidate)
    decisions: list[dict[str, Any]] = []
    shape = np.asarray(candidate.shape, dtype=int)
    for label_id in range(1, count + 1):
        component = labels == label_id
        voxel_count = int(np.count_nonzero(component))
        slices = objects[label_id - 1]
        if slices is None:
            continue
        start = np.asarray([item.start for item in slices], dtype=int)
        stop = np.asarray([item.stop for item in slices], dtype=int)
        extent_um_zyx = (stop - start) * spacing_zyx
        touches_boundary = bool(np.any(start == 0) or np.any(stop == shape))
        keep = voxel_count >= min_component_voxels
        if keep:
            cleaned |= component
        decisions.append(
            {
                "component_id": label_id,
                "voxel_count": voxel_count,
                "volume_um3": voxel_count * voxel_volume,
                "bbox_z_min": int(start[0]),
                "bbox_z_max_exclusive": int(stop[0]),
                "bbox_y_min": int(start[1]),
                "bbox_y_max_exclusive": int(stop[1]),
                "bbox_x_min": int(start[2]),
                "bbox_x_max_exclusive": int(stop[2]),
                "extent_x_um": float(extent_um_zyx[2]),
                "extent_y_um": float(extent_um_zyx[1]),
                "extent_z_um": float(extent_um_zyx[0]),
                "touches_volume_boundary": touches_boundary,
                "decision": "keep" if keep else "remove",
                "reason": (
                    "meets_min_component_voxels"
                    if keep
                    else "below_min_component_voxels"
                ),
            }
        )
    return cleaned, candidate & ~cleaned, decisions
