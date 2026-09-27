from __future__ import annotations

from pathlib import Path

import numpy as np

from utils.config import HierarchicalGraphConfig
from utils.hierarchical_graph_pipeline import run_hierarchical_graph_pipeline
from utils.io import write_json, write_nifti_mask


def _write_source_run(source_run: Path) -> None:
    volumes = source_run / "volumes"
    reports = source_run / "reports"
    volumes.mkdir(parents=True)
    reports.mkdir(parents=True)
    mask = np.zeros((12, 12, 12), dtype=bool)
    skeleton = np.zeros_like(mask)
    skeleton[2:10, 6, 6] = True
    skeleton[6, 7:10, 6] = True
    mask[1:11, 5:8, 5:8] = True
    mask[5:8, 6:11, 5:8] = True
    origin = (-12.0, -12.0, -12.0)
    spacing = (2.0, 2.0, 2.0)
    write_nifti_mask(mask, volumes / "voxel_mask.nii.gz", origin, spacing, "test mask")
    write_nifti_mask(
        skeleton,
        volumes / "coarse_skeleton.nii.gz",
        origin,
        spacing,
        "test skeleton",
    )
    write_json(
        {
            "array_axis_order": ["x", "y", "z"],
            "source_coordinate_system": "LPS",
            "nifti_world_coordinate_system": "RAS",
            "origin_lps_um": origin,
            "spacing_um": spacing,
            "dimensions_xyz": skeleton.shape,
        },
        volumes / "spatial_metadata.json",
    )
    write_json({"test": True}, source_run / "run_config.json")
    write_json({"skeleton_voxel_count": int(skeleton.sum())}, reports / "skeleton_report.json")


def test_step3_pipeline_writes_graph_reports_and_visualizations(tmp_path: Path) -> None:
    source_run = tmp_path / "outputs" / "sample" / "run_source"
    _write_source_run(source_run)
    result = run_hierarchical_graph_pipeline(
        source_run,
        HierarchicalGraphConfig(
            neighbor_connectivity=6,
            smoothing_enabled=False,
            short_branch_warning_um=0.0,
            high_degree_warning=7,
        ),
    )

    assert result.status == "success"
    assert result.html_report.is_file()
    assert (result.run_root / "hierarchical_graph.log").is_file()
    assert (result.run_root / "graphs" / "hierarchical_vascular_graph.json").is_file()
    assert (result.run_root / "graphs" / "branch_as_node_graph.graphml").is_file()
    assert (result.run_root / "graphs" / "branch_geometry.npz").is_file()
    assert (result.run_root / "tables" / "junction_angles.csv").is_file()
    assert (
        result.run_root / "volumes" / "graph_reconstructed_skeleton.nii.gz"
    ).is_file()
    assert (
        result.run_root / "visualizations" / "representation_fidelity_dashboard.png"
    ).is_file()
    assert (result.run_root / "reports" / "graph_acceptance_summary.json").is_file()
