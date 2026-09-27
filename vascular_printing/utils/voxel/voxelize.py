"""Convert a closed vascular surface into a uniform binary voxel mask."""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from ..config import VoxelizationConfig
from .connectivity import keep_largest_voxel_component


@dataclass(slots=True)
class VoxelizationResult:
    mask: np.ndarray
    removed_islands_mask: np.ndarray
    origin_lps_um: tuple[float, float, float]
    spacing_um: tuple[float, float, float]
    dimensions_xyz: tuple[int, int, int]
    initial_foreground_voxel_count: int
    foreground_voxel_count: int
    foreground_fraction: float
    initial_connected_component_count: int
    connected_component_count: int
    removed_island_voxel_count: int
    removed_island_fraction: float
    component_voxel_counts_top20: list[int]
    mask_volume_um3: float

    def report(self) -> dict[str, Any]:
        payload = asdict(self)
        payload.pop("mask")
        payload.pop("removed_islands_mask")
        return payload


def _grid_geometry(
    bounds: tuple[float, float, float, float, float, float],
    config: VoxelizationConfig,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[int, int, int]]:
    spacing = (config.voxel_size_um,) * 3
    padding = config.padding_voxels * config.voxel_size_um
    origin = (bounds[0] - padding, bounds[2] - padding, bounds[4] - padding)
    maxima = (bounds[1] + padding, bounds[3] + padding, bounds[5] + padding)
    dimensions = tuple(
        int(np.ceil((maximum - start) / config.voxel_size_um)) + 1
        for start, maximum in zip(origin, maxima)
    )
    voxel_count = int(np.prod(dimensions, dtype=np.int64))
    if voxel_count > config.max_voxel_count:
        gib = voxel_count / (1024**3)
        raise MemoryError(
            "Requested voxel grid is too large: "
            f"{dimensions} = {voxel_count:,} voxels (~{gib:.2f} GiB per uint8 array). "
            "Increase --voxel-size or --max-voxel-count deliberately."
        )
    return origin, spacing, dimensions


def voxelize_surface(
    mesh: vtk.vtkPolyData,
    config: VoxelizationConfig,
    logger: logging.Logger | None = None,
) -> VoxelizationResult:
    logger = logger or logging.getLogger("ulm_3d_vascular")
    bounds = tuple(float(value) for value in mesh.GetBounds())
    origin, spacing, dimensions = _grid_geometry(bounds, config)
    extent = (0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1)
    logger.info(
        "Voxel grid: dimensions=%s, spacing=%.3f um, total=%s",
        dimensions,
        config.voxel_size_um,
        f"{int(np.prod(dimensions, dtype=np.int64)):,}",
    )

    white_image = vtk.vtkImageData()
    white_image.SetSpacing(spacing)
    white_image.SetOrigin(origin)
    white_image.SetExtent(extent)
    white_image.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
    vtk_to_numpy(white_image.GetPointData().GetScalars()).fill(1)

    surface_to_stencil = vtk.vtkPolyDataToImageStencil()
    surface_to_stencil.SetInputData(mesh)
    surface_to_stencil.SetOutputOrigin(origin)
    surface_to_stencil.SetOutputSpacing(spacing)
    surface_to_stencil.SetOutputWholeExtent(extent)
    surface_to_stencil.Update()

    stencil = vtk.vtkImageStencil()
    stencil.SetInputData(white_image)
    stencil.SetStencilConnection(surface_to_stencil.GetOutputPort())
    stencil.ReverseStencilOff()
    stencil.SetBackgroundValue(0)
    stencil.Update()

    flat = vtk_to_numpy(stencil.GetOutput().GetPointData().GetScalars())
    raw_mask = np.asarray(flat.reshape(dimensions, order="F"), dtype=bool)
    if not np.any(raw_mask):
        raise ValueError("Voxelization produced an empty mask")
    connectivity = keep_largest_voxel_component(
        raw_mask,
        enabled=config.keep_largest_connected_component,
    )
    logger.info(
        "Voxel connectivity: %d component(s) before filtering, %d after; "
        "removed %d voxels (%.4f%%)",
        connectivity.initial_component_count,
        connectivity.final_component_count,
        connectivity.removed_island_voxel_count,
        connectivity.removed_island_fraction * 100,
    )
    total = int(raw_mask.size)
    voxel_volume = float(np.prod(spacing))
    return VoxelizationResult(
        mask=connectivity.main_mask,
        removed_islands_mask=connectivity.removed_islands_mask,
        origin_lps_um=origin,
        spacing_um=spacing,
        dimensions_xyz=dimensions,
        initial_foreground_voxel_count=connectivity.initial_foreground_voxel_count,
        foreground_voxel_count=connectivity.final_foreground_voxel_count,
        foreground_fraction=connectivity.final_foreground_voxel_count / total,
        initial_connected_component_count=connectivity.initial_component_count,
        connected_component_count=connectivity.final_component_count,
        removed_island_voxel_count=connectivity.removed_island_voxel_count,
        removed_island_fraction=connectivity.removed_island_fraction,
        component_voxel_counts_top20=connectivity.component_voxel_counts_top20,
        mask_volume_um3=connectivity.final_foreground_voxel_count * voxel_volume,
    )
