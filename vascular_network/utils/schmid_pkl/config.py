"""Configuration for importing the Schmid NW1 pickle dictionaries."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class SchmidPKLConfig:
    input_dir: Path
    output_root: Path
    pressure_tolerance_mmhg: float = 1.0e-9
    zero_flow_tolerance_um3_per_ms: float = 1.0e-9
    endpoint_tolerance_um: float = 1.0e-6
    length_relative_warning: float = 0.05
    smoothing_enabled: bool = True
    smoothing_window_points: int = 5
    resample_step_um: float = 2.0
    keep_largest_component: bool = True
    visualizations_enabled: bool = True
    write_preview_volume: bool = False
    preview_voxel_size_um: float = 4.0
    max_preview_voxel_count: int = 100_000_000
    max_direction_arrows: int = 800

    def validate(self) -> None:
        if not self.input_dir.is_dir():
            raise FileNotFoundError(f"Input directory does not exist: {self.input_dir}")
        for name in ("verticesDict.pkl", "edgesDict.pkl"):
            path = self.input_dir / name
            if not path.is_file():
                raise FileNotFoundError(f"Required Schmid file is missing: {path}")
        if self.pressure_tolerance_mmhg < 0:
            raise ValueError("pressure_tolerance_mmhg must be non-negative")
        if self.zero_flow_tolerance_um3_per_ms < 0:
            raise ValueError("zero_flow_tolerance_um3_per_ms must be non-negative")
        if self.endpoint_tolerance_um < 0:
            raise ValueError("endpoint_tolerance_um must be non-negative")
        if not 0 <= self.length_relative_warning <= 1:
            raise ValueError("length_relative_warning must be between 0 and 1")
        if self.smoothing_window_points < 3 or self.smoothing_window_points % 2 == 0:
            raise ValueError("smoothing_window_points must be an odd integer >= 3")
        if self.resample_step_um <= 0:
            raise ValueError("resample_step_um must be positive")
        if self.preview_voxel_size_um <= 0:
            raise ValueError("preview_voxel_size_um must be positive")
        if self.max_preview_voxel_count <= 0:
            raise ValueError("max_preview_voxel_count must be positive")
        if self.max_direction_arrows <= 0:
            raise ValueError("max_direction_arrows must be positive")

    def report(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["input_dir"] = str(self.input_dir.resolve())
        payload["output_root"] = str(self.output_root.resolve())
        return payload
