"""Configuration for the NNE2 Step 1-3 equivalent pipeline."""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class NNE2Config:
    input_dir: Path
    output_root: Path
    stage: str = "all"
    source_run: Path | None = None
    subject_id: str | None = None
    tree_id: int | None = None
    stack_name: str | None = None
    max_stacks: int | None = None
    target_xy_spacing_um: float = 2.0
    gaussian_sigma_um: float = 1.2
    background_sigma_um: float = 18.0
    foreground_quantile: float = 0.92
    min_component_voxels: int = 96
    closing_iterations: int = 1
    graph_connectivity: int = 26
    graph_smoothing_enabled: bool = False
    graph_smoothing_window_points: int = 5
    graph_resample_step_um: float | None = None
    junction_direction_distance_um: float = 10.0
    short_branch_warning_um: float = 6.0
    large_junction_warning_voxels: int = 20
    high_degree_warning: int = 6
    smoothing_deviation_warning_um: float = 2.0
    island_warning_fraction: float = 0.01
    island_fail_fraction: float = 0.10
    max_voxel_count: int = 200_000_000
    registration_patch_um: float = 180.0
    registration_z_search: int = 2
    min_registration_score: float = 0.12
    max_anchor_distance_um: float = 45.0
    io_workers: int = max(1, min(8, (os.cpu_count() or 2)))
    visualizations_enabled: bool = True
    write_nifti: bool = True
    save_graphml: bool = True
    save_vtp: bool = True
    save_npz: bool = True
    use_cache: bool = True

    def validate(self) -> None:
        self.input_dir = Path(self.input_dir).expanduser().resolve()
        self.output_root = Path(self.output_root).expanduser().resolve()
        if self.source_run is not None:
            self.source_run = Path(self.source_run).expanduser().resolve()
        if not self.input_dir.is_dir():
            raise FileNotFoundError(f"NNE2 input directory does not exist: {self.input_dir}")
        if not (self.input_dir / "vdb.mat").is_file():
            raise FileNotFoundError(f"Missing vdb.mat in: {self.input_dir}")
        allowed_stages = {"inventory", "preprocess", "hierarchical-graph", "all"}
        if self.stage not in allowed_stages:
            raise ValueError(f"stage must be one of {sorted(allowed_stages)}")
        if self.stage == "hierarchical-graph":
            if self.source_run is None:
                raise ValueError("source_run is required for stage 'hierarchical-graph'")
            if not self.source_run.is_dir():
                raise FileNotFoundError(f"NNE2 source run does not exist: {self.source_run}")
        if self.tree_id is not None and self.tree_id < 0:
            raise ValueError("tree_id must be non-negative")
        if self.max_stacks is not None and self.max_stacks < 1:
            raise ValueError("max_stacks must be positive")
        if self.target_xy_spacing_um <= 0:
            raise ValueError("target_xy_spacing_um must be positive")
        if self.gaussian_sigma_um < 0 or self.background_sigma_um <= 0:
            raise ValueError("filter scales must be non-negative/positive")
        if not 0.5 < self.foreground_quantile < 1.0:
            raise ValueError("foreground_quantile must be between 0.5 and 1")
        if self.min_component_voxels < 1:
            raise ValueError("min_component_voxels must be positive")
        if self.closing_iterations < 0:
            raise ValueError("closing_iterations must be non-negative")
        if self.graph_connectivity not in {6, 18, 26}:
            raise ValueError("graph_connectivity must be 6, 18, or 26")
        if self.graph_smoothing_window_points < 3 or self.graph_smoothing_window_points % 2 == 0:
            raise ValueError("graph_smoothing_window_points must be an odd integer >= 3")
        if self.graph_resample_step_um is not None and self.graph_resample_step_um <= 0:
            raise ValueError("graph_resample_step_um must be positive when provided")
        if self.junction_direction_distance_um <= 0:
            raise ValueError("junction_direction_distance_um must be positive")
        if self.short_branch_warning_um < 0 or self.smoothing_deviation_warning_um < 0:
            raise ValueError("graph warning distances must be non-negative")
        if self.large_junction_warning_voxels < 1 or self.high_degree_warning < 3:
            raise ValueError("graph warning limits are invalid")
        if not 0 <= self.island_warning_fraction <= self.island_fail_fraction <= 1:
            raise ValueError("island fractions must satisfy 0 <= warning <= fail <= 1")
        if self.max_voxel_count < 1:
            raise ValueError("max_voxel_count must be positive")
        if self.registration_patch_um <= 0 or self.registration_z_search < 0:
            raise ValueError("registration settings are invalid")
        if not -1.0 <= self.min_registration_score <= 1.0:
            raise ValueError("min_registration_score must be in [-1, 1]")
        if self.max_anchor_distance_um <= 0:
            raise ValueError("max_anchor_distance_um must be positive")
        if self.io_workers < 1:
            raise ValueError("io_workers must be positive")

    def report(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["input_dir"] = str(self.input_dir)
        payload["output_root"] = str(self.output_root)
        payload["source_run"] = str(self.source_run) if self.source_run is not None else None
        return payload
