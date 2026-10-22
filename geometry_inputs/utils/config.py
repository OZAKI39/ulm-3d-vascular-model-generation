"""Configuration objects for mesh cleanup and voxelization."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class MeshCleanupConfig:
    """Rules used to select the connected surface passed downstream.

    ``main_network_only`` keeps one selected connected component.  The older
    ``conservative`` policy is retained as an explicit comparison mode.
    """

    component_policy: str = "main_network_only"
    main_component_id: int | None = None
    min_component_faces: int = 100
    min_component_area_um2: float = 50.0
    min_component_diagonal_um: float = 10.0
    degenerate_area_epsilon_um2: float = 1.0e-10
    repair_non_manifold: bool = True
    max_repair_face_change_fraction: float = 0.25
    smoothing_iterations: int = 0
    smoothing_pass_band: float = 0.1

    def validate(self) -> None:
        if self.component_policy not in {"main_network_only", "conservative"}:
            raise ValueError(
                "component_policy must be 'main_network_only' or 'conservative'"
            )
        if self.main_component_id is not None and self.main_component_id < 0:
            raise ValueError("main_component_id must be non-negative")
        if self.min_component_faces < 0:
            raise ValueError("min_component_faces must be non-negative")
        if self.min_component_area_um2 < 0:
            raise ValueError("min_component_area_um2 must be non-negative")
        if self.min_component_diagonal_um < 0:
            raise ValueError("min_component_diagonal_um must be non-negative")
        if self.degenerate_area_epsilon_um2 < 0:
            raise ValueError("degenerate_area_epsilon_um2 must be non-negative")
        if not 0 <= self.max_repair_face_change_fraction <= 1:
            raise ValueError("max_repair_face_change_fraction must be in [0, 1]")
        if self.smoothing_iterations < 0:
            raise ValueError("smoothing_iterations must be non-negative")
        if not 0 < self.smoothing_pass_band <= 2:
            raise ValueError("smoothing_pass_band must be in (0, 2]")


@dataclass(slots=True)
class VoxelizationConfig:
    """Settings for the uniform voxel mask and coarse skeleton."""

    voxel_size_um: float = 2.0
    padding_voxels: int = 3
    max_voxel_count: int = 200_000_000
    skeletonize: bool = True
    keep_largest_connected_component: bool = True
    island_warning_fraction: float = 0.01
    island_fail_fraction: float = 0.10

    def validate(self) -> None:
        if self.voxel_size_um <= 0:
            raise ValueError("voxel_size_um must be positive")
        if self.padding_voxels < 1:
            raise ValueError("padding_voxels must be at least 1")
        if self.max_voxel_count <= 0:
            raise ValueError("max_voxel_count must be positive")
        if not 0 <= self.island_warning_fraction <= 1:
            raise ValueError("island_warning_fraction must be in [0, 1]")
        if not 0 <= self.island_fail_fraction <= 1:
            raise ValueError("island_fail_fraction must be in [0, 1]")
        if self.island_warning_fraction > self.island_fail_fraction:
            raise ValueError("island_warning_fraction cannot exceed island_fail_fraction")


@dataclass(slots=True)
class VisualizationConfig:
    enabled: bool = True
    window_width: int = 1400
    window_height: int = 1000
    max_display_faces: int = 350_000
    max_skeleton_points: int = 200_000

    def validate(self) -> None:
        if self.window_width < 400 or self.window_height < 300:
            raise ValueError("visualization window is too small")
        if self.max_display_faces < 1 or self.max_skeleton_points < 1:
            raise ValueError("visualization limits must be positive")


@dataclass(slots=True)
class HierarchicalGraphConfig:
    """Settings for the navigation-quality hierarchical vascular graph."""

    neighbor_connectivity: int = 26
    smoothing_enabled: bool = True
    smoothing_window_points: int = 5
    resample_step_um: float | None = None
    junction_direction_distance_um: float = 10.0
    short_branch_warning_um: float = 6.0
    large_junction_warning_voxels: int = 20
    high_degree_warning: int = 6
    smoothing_deviation_warning_um: float = 2.0
    save_graphml: bool = True
    save_vtp: bool = True
    save_npz: bool = True

    def validate(self) -> None:
        if self.neighbor_connectivity not in {6, 18, 26}:
            raise ValueError("neighbor_connectivity must be 6, 18, or 26")
        if self.smoothing_window_points < 3 or self.smoothing_window_points % 2 == 0:
            raise ValueError("smoothing_window_points must be an odd integer >= 3")
        if self.resample_step_um is not None and self.resample_step_um <= 0:
            raise ValueError("resample_step_um must be positive when provided")
        if self.junction_direction_distance_um <= 0:
            raise ValueError("junction_direction_distance_um must be positive")
        if self.short_branch_warning_um < 0:
            raise ValueError("short_branch_warning_um must be non-negative")
        if self.large_junction_warning_voxels < 1:
            raise ValueError("large_junction_warning_voxels must be positive")
        if self.high_degree_warning < 3:
            raise ValueError("high_degree_warning must be at least 3")
        if self.smoothing_deviation_warning_um < 0:
            raise ValueError("smoothing_deviation_warning_um must be non-negative")


@dataclass(slots=True)
class PipelineConfig:
    input_stl: Path
    output_root: Path
    mesh: MeshCleanupConfig = field(default_factory=MeshCleanupConfig)
    voxel: VoxelizationConfig = field(default_factory=VoxelizationConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    coordinate_system: str = "LPS"
    physical_unit: str = "micrometer"

    def validate(self) -> None:
        self.input_stl = Path(self.input_stl).expanduser().resolve()
        self.output_root = Path(self.output_root).expanduser().resolve()
        if not self.input_stl.is_file():
            raise FileNotFoundError(f"Input STL does not exist: {self.input_stl}")
        if self.input_stl.suffix.lower() != ".stl":
            raise ValueError(f"Input must be an STL file: {self.input_stl}")
        if self.coordinate_system.upper() != "LPS":
            raise ValueError("This pipeline currently expects an LPS STL")
        self.mesh.validate()
        self.voxel.validate()
        self.visualization.validate()

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["input_stl"] = str(self.input_stl)
        payload["output_root"] = str(self.output_root)
        return payload
