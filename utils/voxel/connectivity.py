"""Connected-component filtering for 3-D binary voxel masks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage


@dataclass(slots=True)
class VoxelConnectivityResult:
    main_mask: np.ndarray
    removed_islands_mask: np.ndarray
    initial_component_count: int
    final_component_count: int
    initial_foreground_voxel_count: int
    final_foreground_voxel_count: int
    removed_island_voxel_count: int
    removed_island_fraction: float
    component_voxel_counts_top20: list[int]


def keep_largest_voxel_component(
    mask: np.ndarray,
    *,
    enabled: bool = True,
) -> VoxelConnectivityResult:
    if mask.ndim != 3 or not np.any(mask):
        raise ValueError("A non-empty 3-D mask is required")
    structure = ndimage.generate_binary_structure(3, 3)
    labels, component_count = ndimage.label(mask, structure=structure)
    counts = np.bincount(labels.ravel())
    foreground_counts = counts[1:]
    order = np.argsort(foreground_counts)[::-1]
    top20 = [int(foreground_counts[index]) for index in order[:20]]
    initial_foreground = int(np.count_nonzero(mask))

    if not enabled or component_count <= 1:
        main_mask = np.asarray(mask, dtype=bool)
        removed = np.zeros_like(main_mask, dtype=bool)
        final_component_count = int(component_count)
    else:
        largest_label = int(np.argmax(foreground_counts)) + 1
        main_mask = labels == largest_label
        removed = np.asarray(mask & ~main_mask, dtype=bool)
        final_component_count = 1

    final_foreground = int(np.count_nonzero(main_mask))
    removed_count = initial_foreground - final_foreground
    return VoxelConnectivityResult(
        main_mask=main_mask,
        removed_islands_mask=removed,
        initial_component_count=int(component_count),
        final_component_count=final_component_count,
        initial_foreground_voxel_count=initial_foreground,
        final_foreground_voxel_count=final_foreground,
        removed_island_voxel_count=removed_count,
        removed_island_fraction=removed_count / initial_foreground,
        component_voxel_counts_top20=top20,
    )

