from __future__ import annotations

import numpy as np
import vtk

from utils.config import VoxelizationConfig
from utils.io import lps_to_ras_affine
from utils.voxel.skeleton import extract_coarse_skeleton
from utils.voxel.voxelize import voxelize_surface


def _sphere(
    radius: float = 6.0, center: tuple[float, float, float] = (0.0, 0.0, 0.0)
) -> vtk.vtkPolyData:
    source = vtk.vtkSphereSource()
    source.SetRadius(radius)
    source.SetCenter(center)
    source.SetThetaResolution(32)
    source.SetPhiResolution(32)
    source.Update()
    return source.GetOutput()


def test_voxelization_and_skeleton_are_nonempty_and_nested() -> None:
    result = voxelize_surface(
        _sphere(), VoxelizationConfig(voxel_size_um=1.0, padding_voxels=2)
    )
    skeleton = extract_coarse_skeleton(result.mask, result.spacing_um)

    assert result.mask.ndim == 3
    assert result.foreground_voxel_count > 0
    assert result.foreground_fraction < 0.5
    assert result.connected_component_count == 1
    assert result.initial_connected_component_count == 1
    assert skeleton.skeleton_voxel_count > 0
    assert skeleton.voxels_outside_mask == 0
    assert np.all(result.mask[skeleton.skeleton])


def test_voxelization_keeps_only_largest_connected_voxel_network() -> None:
    append = vtk.vtkAppendPolyData()
    append.AddInputData(_sphere(6.0, (0.0, 0.0, 0.0)))
    append.AddInputData(_sphere(2.0, (20.0, 0.0, 0.0)))
    append.Update()

    result = voxelize_surface(
        append.GetOutput(), VoxelizationConfig(voxel_size_um=1.0, padding_voxels=2)
    )

    assert result.initial_connected_component_count == 2
    assert result.connected_component_count == 1
    assert result.removed_island_voxel_count > 0
    assert np.any(result.removed_islands_mask)
    assert not np.any(result.mask & result.removed_islands_mask)


def test_lps_to_ras_affine_flips_x_and_y_world_axes() -> None:
    affine = lps_to_ras_affine((-10.0, -20.0, 5.0), (2.0, 3.0, 4.0))
    np.testing.assert_allclose(affine[:3, :3], np.diag([-2.0, -3.0, 4.0]))
    np.testing.assert_allclose(affine[:3, 3], [10.0, 20.0, 5.0])
