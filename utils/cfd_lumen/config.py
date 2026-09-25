"""Configuration for the single supported Ultraliser reconstruction workflow."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(slots=True)
class UltraliserConfig:
    enabled: bool = True
    ultraliser_root: str = "Ultraliser"
    radius_scale: float = 0.91
    voxels_per_micron: float = 6.0
    packing_algorithm: str = "polylines-with-spheres"
    voxelization_axis: str = "xyz"
    isosurface_technique: str = "dmc"
    solid_voxelization: bool = True
    adaptive_optimization: bool = True
    optimization_iterations: int = 5
    smooth_iterations: int = 5
    laplacian_iterations: int = 10
    threads: int = 8
    export_stl: bool = True


@dataclass(slots=True)
class SurfaceQCConfig:
    radius_fidelity_samples_per_branch: int = 10
    radius_fidelity_skip_diameters: float = 2.0
    max_radius_p95_error: float = 0.05
    require_watertight: bool = True
    require_single_component: bool = True
    require_zero_boundary_edges: bool = True
    require_zero_nonmanifold_edges: bool = True
    require_zero_self_intersections: bool = True
    require_zero_degenerate_triangles: bool = True


@dataclass(slots=True)
class CFDLumenConfig:
    ultraliser: UltraliserConfig = field(default_factory=UltraliserConfig)
    surface_qc: SurfaceQCConfig = field(default_factory=SurfaceQCConfig)

    def validate(self) -> None:
        settings = self.ultraliser
        if not settings.enabled:
            raise ValueError("Ultraliser is the only supported reconstruction backend")
        if settings.radius_scale <= 0.0:
            raise ValueError("ultraliser.radius_scale must be positive")
        if settings.voxels_per_micron <= 0.0:
            raise ValueError("ultraliser.voxels_per_micron must be positive")
        if settings.packing_algorithm != "polylines-with-spheres":
            raise ValueError("only polylines-with-spheres packing is supported")
        if settings.voxelization_axis != "xyz":
            raise ValueError("only xyz voxelization is supported")
        if settings.isosurface_technique != "dmc":
            raise ValueError("only dmc isosurface extraction is supported")
        if not settings.solid_voxelization or not settings.adaptive_optimization:
            raise ValueError("solid voxelization and adaptive optimization must remain enabled")
        if min(
            settings.optimization_iterations,
            settings.smooth_iterations,
            settings.laplacian_iterations,
        ) < 0:
            raise ValueError("Ultraliser iteration counts cannot be negative")
        if settings.threads < 1:
            raise ValueError("ultraliser.threads must be positive")
        if not settings.export_stl:
            raise ValueError("STL export is required by the reconstruction workflow")
        qc = self.surface_qc
        if qc.radius_fidelity_samples_per_branch < 1:
            raise ValueError("radius_fidelity_samples_per_branch must be positive")
        if qc.radius_fidelity_skip_diameters < 0.0:
            raise ValueError("radius_fidelity_skip_diameters cannot be negative")
        if not 0.0 < qc.max_radius_p95_error < 1.0:
            raise ValueError("max_radius_p95_error must be in (0, 1)")

    def report(self) -> dict[str, Any]:
        return asdict(self)


_SECTIONS: dict[str, type[Any]] = {
    "ultraliser": UltraliserConfig,
    "surface_qc": SurfaceQCConfig,
}


def load_cfd_lumen_config(path: Path | None = None) -> CFDLumenConfig:
    """Load a partial YAML configuration over the validated formal defaults."""

    payload: dict[str, Any] = {}
    if path is not None:
        config_path = Path(path).resolve()
        if not config_path.is_file():
            raise FileNotFoundError(f"Ultraliser configuration not found: {config_path}")
        loaded = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        if loaded is not None and not isinstance(loaded, dict):
            raise ValueError("Configuration root must be a mapping")
        payload = loaded or {}
    unknown = sorted(set(payload) - set(_SECTIONS))
    if unknown:
        raise ValueError(f"Unknown configuration section(s): {', '.join(unknown)}")
    defaults = CFDLumenConfig()
    values: dict[str, Any] = {}
    for name, section_type in _SECTIONS.items():
        section_values = asdict(getattr(defaults, name))
        override = payload.get(name, {})
        if not isinstance(override, dict):
            raise ValueError(f"Configuration section {name!r} must be a mapping")
        section_values.update(override)
        try:
            values[name] = section_type(**section_values)
        except TypeError as exc:
            raise ValueError(f"Invalid key in configuration section {name!r}: {exc}") from exc
    config = CFDLumenConfig(**values)
    config.validate()
    return config
