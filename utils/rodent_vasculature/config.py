"""Configuration for the directed rodent vasculature workflow."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class RodentVasculatureConfig:
    input_dir: Path
    output_root: Path
    stage: str = "inventory"
    source_run: Path | None = None
    cohort: str = "raw-total"
    sample_id: str | None = None
    parent_group_id: str | None = None
    split: str | None = None
    max_samples: int | None = None
    spacing_xyz_um: tuple[float, float, float] = (1.0, 1.0, 2.0)
    expected_shape_zyx: tuple[int, int, int] | None = (192, 192, 192)
    smoothing_enabled: bool = False
    smoothing_window_points: int = 5
    resample_step_um: float = 1.0
    strict_nonpositive_radius: bool = False
    analysis_component_id: int | None = None
    visualizations_enabled: bool = True
    max_visualization_samples: int = 3
    max_direction_arrows: int = 600
    figure2a_enabled: bool = False
    figure2a_volume_opacity: float = 0.32
    figure2a_window_size: tuple[int, int] = (1800, 900)
    save_graphml: bool = True
    save_vtp: bool = True
    save_npz: bool = True

    def validate(self) -> None:
        if self.stage not in {"inventory", "preprocess", "hierarchical-graph", "all"}:
            raise ValueError(f"Unsupported stage: {self.stage}")
        if self.cohort not in {"raw-total", "raw-analysis", "train", "all"}:
            raise ValueError(f"Unsupported cohort: {self.cohort}")
        if self.stage == "hierarchical-graph" and self.source_run is None:
            raise ValueError("source_run is required for the hierarchical-graph stage")
        if self.max_samples is not None and self.max_samples <= 0:
            raise ValueError("max_samples must be positive")
        if any(value <= 0 for value in self.spacing_xyz_um):
            raise ValueError("spacing_xyz_um values must be positive")
        if self.expected_shape_zyx is not None and any(
            value <= 0 for value in self.expected_shape_zyx
        ):
            raise ValueError("expected_shape_zyx values must be positive")
        if self.smoothing_window_points < 3 or self.smoothing_window_points % 2 == 0:
            raise ValueError("smoothing_window_points must be an odd integer >= 3")
        if self.resample_step_um <= 0:
            raise ValueError("resample_step_um must be positive")
        if self.analysis_component_id is not None and self.analysis_component_id < 0:
            raise ValueError("analysis_component_id must be non-negative")
        if self.max_visualization_samples < 0:
            raise ValueError("max_visualization_samples cannot be negative")
        if self.max_direction_arrows <= 0:
            raise ValueError("max_direction_arrows must be positive")
        if not 0.0 < self.figure2a_volume_opacity <= 1.0:
            raise ValueError("figure2a_volume_opacity must be in (0, 1]")
        if any(value <= 0 for value in self.figure2a_window_size):
            raise ValueError("figure2a_window_size values must be positive")
        if self.stage != "hierarchical-graph" and not self.input_dir.is_dir():
            raise FileNotFoundError(f"Input dataset directory does not exist: {self.input_dir}")

    def report(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["input_dir"] = str(self.input_dir)
        payload["output_root"] = str(self.output_root)
        payload["source_run"] = str(self.source_run) if self.source_run else None
        payload["flow_direction_rule"] = "SWC parent_id node -> current node"
        payload["flow_direction_is_measured"] = False
        return payload
