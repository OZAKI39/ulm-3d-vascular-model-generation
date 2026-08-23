"""Validated configuration for SWC-ROI to CFD-lumen reconstruction."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(slots=True)
class GeometryConfig:
    source_unit: str = "um"
    export_unit: str = "m"
    centerline_smoothing: bool = False
    smoothing_window_points: int = 3
    resample_radius_fraction: float = 0.35
    min_resample_spacing_um: float = 0.2
    max_resample_spacing_um: float = 1.0
    minimum_edge_length_um: float = 1.0e-8
    tube_sides: int = 32


@dataclass(slots=True)
class LocalJunctionImplicitConfig:
    cells_across_min_diameter: int = 16
    dtype: str = "float32"
    smooth_union: bool = False
    chunk_size: int = 500_000
    k_nearest: int = 16
    max_grid_cells: int = 20_000_000
    direction_clip_epsilon_radius_fraction: float = 0.10


@dataclass(slots=True)
class JunctionRemeshConfig:
    enabled: bool = False


@dataclass(slots=True)
class JunctionConfig:
    enabled: bool = True
    backend: str = "local_implicit"
    radius_scale: float = 1.0
    collar_diameters: float = 2.0
    overlap_diameters: float = 0.5
    bbox_padding_radius: float = 1.5
    separation_tolerance_fraction: float = 0.05
    implicit: LocalJunctionImplicitConfig = field(default_factory=LocalJunctionImplicitConfig)
    remesh: JunctionRemeshConfig = field(default_factory=JunctionRemeshConfig)


@dataclass(slots=True)
class ReconstructionConfig:
    branch_backend: str = "explicit"
    port_backend: str = "explicit"
    junction_backend: str = "local_implicit"


@dataclass(slots=True)
class HybridMergeConfig:
    backend: str = "manifold"
    cleanup_tolerance_um: float = 3.0e-3


@dataclass(slots=True)
class PortTransitionConfig:
    backend: str = "continuous_centerline"
    profile_diameter_offsets: tuple[float, ...] = (-1.0, 0.0, 1.0, 2.0)


@dataclass(slots=True)
class V6Config:
    enabled: bool = True
    transition_backend: str = "continuous_implicit_field"
    legacy_transition_backend: str = "loop_stitch_legacy"
    port_source_fit_points: int = 5
    tangent_blend_diameters: float = 1.5
    tangent_conditioning_threshold_deg: float = 3.0
    radius_blend_diameters: float = 1.5
    radius_slope_conditioning_threshold_um_per_um: float = 0.02
    port_extension_spacing_radius_fraction: float = 0.05
    port_extension_min_spacing_um: float = 0.05
    pure_branch_overlap_diameters: float = 0.5
    convergence_cells_across_min_diameter: tuple[int, ...] = (16, 20, 24)
    interface_plane_tolerance_radius_fraction: float = 0.12
    silhouette_large_corner_deg: float = 25.0
    loop_center_tolerance_radius_fraction: float = 0.08
    loop_radius_tolerance_fraction: float = 0.05
    loop_hausdorff_tolerance_radius_fraction: float = 0.10
    local_patch_cleanup_tolerance_um: float = 5.0e-4


@dataclass(slots=True)
class V7Config:
    enabled: bool = False
    backend: str = "unified_polyball"
    polyball_provider: str = "auto"
    cells_across_min_diameter: int = 16
    convergence_cells_across_min_diameter: tuple[int, ...] = (12, 16, 20, 24)
    extraction_backend: str = "flying_edges"
    compare_marching_cubes: bool = True
    port_tail_diameters: float = 1.0
    k_nearest_segments: int = 32
    field_padding_cells: int = 3
    max_grid_cells: int = 1_000_000_000
    memory_map_threshold_cells: int = 120_000_000
    memory_map_directory: str | None = None
    newton_iterations: int = 3
    projection_tolerance_um: float = 1.0e-5
    remesh_backend: str = "pyacvd"
    remesh_target_edge_spacing_factor: float = 1.5
    remesh_minimum_clusters: int = 1_000
    oriented_grid: bool = True


@dataclass(slots=True)
class V8Config:
    enabled: bool = False
    backend: str = "unified_polyball_smooth_junction"
    k_radius_ratios: tuple[float, ...] = (0.10, 0.20, 0.30)
    blend_length_scale: float = 1.0
    ownership_switch_margin_fraction: float = 0.02
    sawtooth_dihedral_percentile: float = 95.0
    switch_overlap_distance_edge_factors: float = 1.5
    junction_volume_samples_across_diameter: int = 64
    maximum_junction_volume_increase_fraction: float = 0.05
    maximum_radius_p95_error: float = 0.01
    maximum_collar_radius_error_increase: float = 0.01
    maximum_hydraulic_resistance_error: float = 0.05
    minimum_silhouette_roughness_reduction_fraction: float = 0.05


@dataclass(slots=True)
class HybridTransitionConfig:
    backend: str = "loop_stitch"
    fallback_backend: str = "manifold_boolean"
    transition_rings: int = 6
    smoothstep: str = "quintic"
    constrained_smoothing: bool = False
    smoothing_iterations: int = 5
    smoothing_lambda: float = 0.25
    smoothing_mu: float = -0.26
    capture_v4_comparison: bool = True
    minimum_normal_p99_reduction_fraction: float = 0.20
    minimum_port_normal_p99_reduction_fraction: float = 0.20
    minimum_transition_roughness_reduction_fraction: float = 0.10
    maximum_port_area_error_increase: float = 0.002
    maximum_radius_p95_error_increase: float = 0.002
    maximum_collar_error_increase: float = 0.002


@dataclass(slots=True)
class ContextDomainConfig:
    enabled: bool = True
    source_rodent_run: str | None = None
    max_added_global_edges: int = 1_000
    verify_source_geometry: bool = True


@dataclass(slots=True)
class BranchLocalQCConfig:
    enabled: bool = True
    ownership_margin_threshold: float = 0.10
    minimum_ownership_confidence: float = 0.80
    cross_section_cells_across_source_diameter: int = 96
    maximum_cross_section_grid_points: int = 262_144
    severe_area_relative_error: float = 0.25
    minimum_junction_total_area_ratio: float = 0.20
    maximum_junction_normal_jump_p99_deg: float = 120.0


@dataclass(slots=True)
class PortConfig:
    extension_diameters: float = 5.0
    overlap_diameters: float = 0.5
    plane_tolerance_um: float = 0.02
    radial_tolerance_fraction: float = 0.03
    area_relative_tolerance: float = 0.15
    minimum_normal_alignment: float = 0.90


@dataclass(slots=True)
class CollisionConfig:
    enabled: bool = True
    hard_collision_tolerance_um: float = 0.0
    near_contact_tolerance_um: float = 0.2


@dataclass(slots=True)
class BooleanConfig:
    backend: str = "manifold"
    allow_implicit_fallback: bool = True


@dataclass(slots=True)
class ImplicitConfig:
    dtype: str = "float32"
    cells_across_min_diameter: int = 8
    chunk_size: int = 1_000_000
    k_nearest: int = 16
    min_spacing_um: float = 0.15
    max_spacing_um: float = 1.0
    max_grid_cells: int = 80_000_000


@dataclass(slots=True)
class SurfaceQCConfig:
    check_watertight: bool = True
    check_manifold: bool = True
    check_connected: bool = True
    check_port_area: bool = True
    radius_fidelity_samples_per_branch: int = 10
    radius_fidelity_skip_diameters: float = 2.0
    max_radius_p95_error: float | None = None
    require_watertight: bool = True
    require_single_component: bool = True
    require_zero_boundary_edges: bool = True
    require_zero_nonmanifold_edges: bool = True
    require_zero_self_intersections: bool = True
    require_zero_internal_faces: bool = True
    require_zero_internal_caps: bool = True


@dataclass(slots=True)
class ConvergenceConfig:
    enabled: bool = True
    tube_sides: tuple[int, ...] = (16, 24, 32, 48)


@dataclass(slots=True)
class VolumeMeshConfig:
    enabled: bool = False
    characteristic_length_factor: float = 1.0


@dataclass(slots=True)
class OutputConfig:
    visualizations: bool = True
    save_centerlines: bool = True


@dataclass(slots=True)
class CFDLumenConfig:
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    junction: JunctionConfig = field(default_factory=JunctionConfig)
    reconstruction: ReconstructionConfig = field(default_factory=ReconstructionConfig)
    hybrid_merge: HybridMergeConfig = field(default_factory=HybridMergeConfig)
    port_transition: PortTransitionConfig = field(default_factory=PortTransitionConfig)
    v6: V6Config = field(default_factory=V6Config)
    v7: V7Config = field(default_factory=V7Config)
    v8: V8Config = field(default_factory=V8Config)
    hybrid_transition: HybridTransitionConfig = field(default_factory=HybridTransitionConfig)
    context_domain: ContextDomainConfig = field(default_factory=ContextDomainConfig)
    branch_local_qc: BranchLocalQCConfig = field(default_factory=BranchLocalQCConfig)
    ports: PortConfig = field(default_factory=PortConfig)
    collision_qc: CollisionConfig = field(default_factory=CollisionConfig)
    boolean: BooleanConfig = field(default_factory=BooleanConfig)
    implicit_fallback: ImplicitConfig = field(default_factory=ImplicitConfig)
    surface_qc: SurfaceQCConfig = field(default_factory=SurfaceQCConfig)
    convergence: ConvergenceConfig = field(default_factory=ConvergenceConfig)
    volume_mesh: VolumeMeshConfig = field(default_factory=VolumeMeshConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    def validate(self) -> None:
        geometry = self.geometry
        if geometry.source_unit != "um" or geometry.export_unit != "m":
            raise ValueError("The supported reconstruction/export units are exactly um and m")
        if geometry.resample_radius_fraction <= 0:
            raise ValueError("resample_radius_fraction must be positive")
        if not 0 < geometry.min_resample_spacing_um <= geometry.max_resample_spacing_um:
            raise ValueError("resampling spacing bounds are invalid")
        if geometry.minimum_edge_length_um <= 0 or geometry.tube_sides < 8:
            raise ValueError("minimum_edge_length_um must be positive and tube_sides >= 8")
        if self.junction.radius_scale != 1.0:
            raise ValueError("junction.radius_scale must remain 1.0 to preserve the source SWC radius")
        if self.junction.backend not in {"local_implicit", "legacy_sphere"}:
            raise ValueError("junction.backend must be local_implicit or legacy_sphere")
        if min(
            self.junction.collar_diameters,
            self.junction.overlap_diameters,
            self.junction.bbox_padding_radius,
        ) <= 0:
            raise ValueError("junction collar, overlap, and bbox padding must be positive")
        if self.junction.implicit.dtype != "float32":
            raise ValueError("local junction implicit grids are intentionally float32")
        if min(
            self.junction.implicit.cells_across_min_diameter,
            self.junction.implicit.chunk_size,
            self.junction.implicit.k_nearest,
            self.junction.implicit.max_grid_cells,
        ) < 1:
            raise ValueError("local junction implicit parameters must be positive")
        if self.junction.implicit.direction_clip_epsilon_radius_fraction < 0:
            raise ValueError(
                "junction.implicit.direction_clip_epsilon_radius_fraction cannot be negative"
            )
        if (
            self.reconstruction.branch_backend != "explicit"
            or self.reconstruction.port_backend != "explicit"
            or self.reconstruction.junction_backend not in {"local_implicit", "legacy_sphere"}
        ):
            raise ValueError("supported reconstruction is explicit branch/port plus local_implicit or legacy_sphere junction")
        if self.reconstruction.junction_backend != self.junction.backend:
            raise ValueError("reconstruction.junction_backend must match junction.backend")
        if self.hybrid_merge.backend != "manifold":
            raise ValueError("hybrid_merge.backend must be manifold")
        if self.hybrid_merge.cleanup_tolerance_um < 0:
            raise ValueError("hybrid_merge.cleanup_tolerance_um cannot be negative")
        if self.port_transition.backend != "continuous_centerline":
            raise ValueError("port_transition.backend must be continuous_centerline")
        if tuple(self.port_transition.profile_diameter_offsets) != (-1.0, 0.0, 1.0, 2.0):
            raise ValueError(
                "port_transition.profile_diameter_offsets must remain [-1, 0, 1, 2]"
            )
        v6 = self.v6
        if v6.transition_backend != "continuous_implicit_field":
            raise ValueError("v6.transition_backend must be continuous_implicit_field")
        if v6.legacy_transition_backend != "loop_stitch_legacy":
            raise ValueError("v6.legacy_transition_backend must be loop_stitch_legacy")
        if v6.port_source_fit_points < 3:
            raise ValueError("v6.port_source_fit_points must be >= 3")
        if min(
            v6.tangent_blend_diameters,
            v6.tangent_conditioning_threshold_deg,
            v6.radius_blend_diameters,
            v6.radius_slope_conditioning_threshold_um_per_um,
            v6.port_extension_spacing_radius_fraction,
            v6.port_extension_min_spacing_um,
            v6.pure_branch_overlap_diameters,
            v6.interface_plane_tolerance_radius_fraction,
            v6.silhouette_large_corner_deg,
            v6.loop_center_tolerance_radius_fraction,
            v6.loop_radius_tolerance_fraction,
            v6.loop_hausdorff_tolerance_radius_fraction,
            v6.local_patch_cleanup_tolerance_um,
        ) <= 0:
            raise ValueError("v6 physical lengths and tolerances must be positive")
        if tuple(v6.convergence_cells_across_min_diameter) != (16, 20, 24):
            raise ValueError(
                "v6.convergence_cells_across_min_diameter must remain [16, 20, 24]"
            )
        v7 = self.v7
        if v7.backend != "unified_polyball":
            raise ValueError("v7.backend must be unified_polyball")
        if v7.polyball_provider not in {"auto", "vmtk", "ckdtree"}:
            raise ValueError("v7.polyball_provider must be auto, vmtk, or ckdtree")
        if v7.extraction_backend not in {"flying_edges", "marching_cubes"}:
            raise ValueError(
                "v7.extraction_backend must be flying_edges or marching_cubes"
            )
        if v7.remesh_backend not in {"pyacvd", "none"}:
            raise ValueError("v7.remesh_backend must be pyacvd or none")
        if tuple(v7.convergence_cells_across_min_diameter) != (12, 16, 20, 24):
            raise ValueError(
                "v7.convergence_cells_across_min_diameter must remain [12, 16, 20, 24]"
            )
        if min(
            v7.cells_across_min_diameter,
            v7.port_tail_diameters,
            v7.k_nearest_segments,
            v7.field_padding_cells,
            v7.max_grid_cells,
            v7.memory_map_threshold_cells,
            v7.newton_iterations,
            v7.projection_tolerance_um,
            v7.remesh_target_edge_spacing_factor,
            v7.remesh_minimum_clusters,
        ) <= 0:
            raise ValueError("v7 field, projection, tail, and remesh values must be positive")
        v8 = self.v8
        if v8.backend != "unified_polyball_smooth_junction":
            raise ValueError(
                "v8.backend must be unified_polyball_smooth_junction"
            )
        if tuple(v8.k_radius_ratios) != (0.10, 0.20, 0.30):
            raise ValueError("v8.k_radius_ratios must remain [0.10, 0.20, 0.30]")
        if min(
            v8.blend_length_scale,
            v8.ownership_switch_margin_fraction,
            v8.switch_overlap_distance_edge_factors,
            v8.junction_volume_samples_across_diameter,
            v8.maximum_radius_p95_error,
        ) <= 0:
            raise ValueError("v8 field and diagnostic values must be positive")
        if not 0.0 < v8.sawtooth_dihedral_percentile < 100.0:
            raise ValueError("v8.sawtooth_dihedral_percentile must be in (0, 100)")
        for name in (
            "maximum_junction_volume_increase_fraction",
            "maximum_collar_radius_error_increase",
            "maximum_hydraulic_resistance_error",
            "minimum_silhouette_roughness_reduction_fraction",
        ):
            value = float(getattr(v8, name))
            if not 0.0 <= value < 1.0:
                raise ValueError(f"v8.{name} must be in [0, 1)")
        transition = self.hybrid_transition
        if transition.backend not in {"loop_stitch", "manifold_boolean"}:
            raise ValueError(
                "hybrid_transition.backend must be loop_stitch or manifold_boolean"
            )
        if transition.fallback_backend != "manifold_boolean":
            raise ValueError(
                "hybrid_transition.fallback_backend must be manifold_boolean"
            )
        if not transition.capture_v4_comparison:
            raise ValueError(
                "hybrid_transition.capture_v4_comparison must remain true for v5 acceptance"
            )
        if transition.transition_rings < 2:
            raise ValueError("hybrid_transition.transition_rings must be >= 2")
        if transition.smoothstep not in {"cubic", "quintic"}:
            raise ValueError("hybrid_transition.smoothstep must be cubic or quintic")
        if transition.smoothing_iterations < 0:
            raise ValueError("hybrid transition smoothing iterations cannot be negative")
        if transition.smoothing_lambda <= 0 or transition.smoothing_mu >= 0:
            raise ValueError("hybrid transition Taubin lambda must be > 0 and mu < 0")
        for name in (
            "minimum_normal_p99_reduction_fraction",
            "minimum_port_normal_p99_reduction_fraction",
            "minimum_transition_roughness_reduction_fraction",
        ):
            value = float(getattr(transition, name))
            if not 0.0 <= value < 1.0:
                raise ValueError(f"hybrid_transition.{name} must be in [0, 1)")
        for name in (
            "maximum_port_area_error_increase",
            "maximum_radius_p95_error_increase",
            "maximum_collar_error_increase",
        ):
            if float(getattr(transition, name)) < 0.0:
                raise ValueError(f"hybrid_transition.{name} cannot be negative")
        if self.context_domain.max_added_global_edges < 1:
            raise ValueError("context_domain.max_added_global_edges must be positive")
        branch_qc = self.branch_local_qc
        if branch_qc.ownership_margin_threshold < 0:
            raise ValueError("branch_local_qc.ownership_margin_threshold cannot be negative")
        if not 0 < branch_qc.minimum_ownership_confidence <= 1:
            raise ValueError("branch_local_qc.minimum_ownership_confidence must be in (0, 1]")
        if branch_qc.cross_section_cells_across_source_diameter < 16:
            raise ValueError(
                "branch_local_qc.cross_section_cells_across_source_diameter must be >= 16"
            )
        if branch_qc.maximum_cross_section_grid_points < 1_024:
            raise ValueError("branch_local_qc.maximum_cross_section_grid_points must be >= 1024")
        if branch_qc.severe_area_relative_error <= 0:
            raise ValueError("branch_local_qc.severe_area_relative_error must be positive")
        if not 0 < branch_qc.minimum_junction_total_area_ratio < 1:
            raise ValueError("branch_local_qc.minimum_junction_total_area_ratio must be in (0, 1)")
        if not 0 < branch_qc.maximum_junction_normal_jump_p99_deg <= 180:
            raise ValueError(
                "branch_local_qc.maximum_junction_normal_jump_p99_deg must be in (0, 180]"
            )
        if self.ports.extension_diameters <= 0 or self.ports.overlap_diameters < 0:
            raise ValueError("port extension must be positive and overlap non-negative")
        if self.collision_qc.hard_collision_tolerance_um < 0:
            raise ValueError("hard_collision_tolerance_um cannot be negative")
        if self.collision_qc.near_contact_tolerance_um < 0:
            raise ValueError("near_contact_tolerance_um cannot be negative")
        if self.boolean.backend not in {"manifold", "implicit"}:
            raise ValueError("boolean.backend must be manifold or implicit")
        implicit = self.implicit_fallback
        if implicit.dtype != "float32":
            raise ValueError("implicit fallback is intentionally restricted to float32")
        if min(implicit.cells_across_min_diameter, implicit.chunk_size, implicit.k_nearest) < 1:
            raise ValueError("implicit grid parameters must be positive")
        if not 0 < implicit.min_spacing_um <= implicit.max_spacing_um:
            raise ValueError("implicit spacing bounds are invalid")
        if self.surface_qc.radius_fidelity_samples_per_branch < 1:
            raise ValueError("radius_fidelity_samples_per_branch must be positive")
        if any(int(value) < 8 for value in self.convergence.tube_sides):
            raise ValueError("all convergence tube_sides values must be >= 8")

    def report(self) -> dict[str, Any]:
        return asdict(self)


_SECTIONS: dict[str, type[Any]] = {
    "geometry": GeometryConfig,
    "junction": JunctionConfig,
    "reconstruction": ReconstructionConfig,
    "hybrid_merge": HybridMergeConfig,
    "port_transition": PortTransitionConfig,
    "v6": V6Config,
    "v7": V7Config,
    "v8": V8Config,
    "hybrid_transition": HybridTransitionConfig,
    "context_domain": ContextDomainConfig,
    "branch_local_qc": BranchLocalQCConfig,
    "ports": PortConfig,
    "collision_qc": CollisionConfig,
    "boolean": BooleanConfig,
    "implicit_fallback": ImplicitConfig,
    "surface_qc": SurfaceQCConfig,
    "convergence": ConvergenceConfig,
    "volume_mesh": VolumeMeshConfig,
    "output": OutputConfig,
}


def load_cfd_lumen_config(path: Path | None = None) -> CFDLumenConfig:
    """Load a partial YAML configuration over the validated scientific defaults."""

    payload: dict[str, Any] = {}
    if path is not None:
        config_path = Path(path).resolve()
        if not config_path.is_file():
            raise FileNotFoundError(f"CFD lumen configuration not found: {config_path}")
        loaded = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        if loaded is not None and not isinstance(loaded, dict):
            raise ValueError("CFD lumen configuration root must be a mapping")
        payload = loaded or {}
    unknown = sorted(set(payload) - set(_SECTIONS))
    if unknown:
        raise ValueError(f"Unknown CFD lumen config section(s): {', '.join(unknown)}")
    defaults = CFDLumenConfig()
    sections: dict[str, Any] = {}
    for name, section_type in _SECTIONS.items():
        values = asdict(getattr(defaults, name))
        override = payload.get(name, {})
        if not isinstance(override, dict):
            raise ValueError(f"Config section {name!r} must be a mapping")
        values.update(override)
        if name == "convergence" and isinstance(values.get("tube_sides"), list):
            values["tube_sides"] = tuple(int(value) for value in values["tube_sides"])
        if name == "port_transition" and isinstance(
            values.get("profile_diameter_offsets"), list
        ):
            values["profile_diameter_offsets"] = tuple(
                float(value) for value in values["profile_diameter_offsets"]
            )
        if name == "v6" and isinstance(
            values.get("convergence_cells_across_min_diameter"), list
        ):
            values["convergence_cells_across_min_diameter"] = tuple(
                int(value)
                for value in values["convergence_cells_across_min_diameter"]
            )
        if name == "v7" and isinstance(
            values.get("convergence_cells_across_min_diameter"), list
        ):
            values["convergence_cells_across_min_diameter"] = tuple(
                int(value)
                for value in values["convergence_cells_across_min_diameter"]
            )
        if name == "v8" and isinstance(values.get("k_radius_ratios"), list):
            values["k_radius_ratios"] = tuple(
                float(value) for value in values["k_radius_ratios"]
            )
        if name == "junction":
            implicit_values = asdict(LocalJunctionImplicitConfig())
            implicit_values.update(values.pop("implicit", {}) or {})
            remesh_values = asdict(JunctionRemeshConfig())
            remesh_values.update(values.pop("remesh", {}) or {})
            values["implicit"] = LocalJunctionImplicitConfig(**implicit_values)
            values["remesh"] = JunctionRemeshConfig(**remesh_values)
        try:
            sections[name] = section_type(**values)
        except TypeError as exc:
            raise ValueError(f"Invalid key in config section {name!r}: {exc}") from exc
    config = CFDLumenConfig(**sections)
    config.validate()
    return config
