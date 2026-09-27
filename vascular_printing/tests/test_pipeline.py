from __future__ import annotations

from pathlib import Path

import vtk

from utils.config import PipelineConfig, VisualizationConfig, VoxelizationConfig
from utils.pipeline import run_pipeline


def _write_sphere(path: Path) -> None:
    source = vtk.vtkSphereSource()
    source.SetRadius(8.0)
    source.SetThetaResolution(28)
    source.SetPhiResolution(28)
    source.Update()
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(path))
    writer.SetFileTypeToBinary()
    writer.SetInputData(source.GetOutput())
    assert writer.Write() == 1


def test_small_end_to_end_pipeline_writes_acceptance_artifacts(tmp_path: Path) -> None:
    input_path = tmp_path / "sphere.stl"
    _write_sphere(input_path)
    result = run_pipeline(
        PipelineConfig(
            input_stl=input_path,
            output_root=tmp_path / "outputs",
            voxel=VoxelizationConfig(voxel_size_um=1.0, padding_voxels=2),
            visualization=VisualizationConfig(enabled=False),
        )
    )

    assert result.status in {"success", "warning"}
    assert result.html_report.is_file()
    assert (result.run_root / "pipeline.log").is_file()
    assert (result.run_root / "meshes" / "cleaned.stl").is_file()
    assert (result.run_root / "volumes" / "voxel_mask.nii.gz").is_file()
    assert (result.run_root / "volumes" / "coarse_skeleton.nii.gz").is_file()
    assert (result.run_root / "reports" / "acceptance_summary.json").is_file()

