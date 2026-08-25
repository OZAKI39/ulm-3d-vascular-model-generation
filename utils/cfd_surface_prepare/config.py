"""Strict YAML configuration for local CFD surface preparation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml


@dataclass(frozen=True, slots=True)
class PathsConfig:
    cfd_preprocess_run: Path
    output_root: Path


@dataclass(frozen=True, slots=True)
class GeometryConfig:
    create_meter_copy: bool


@dataclass(frozen=True, slots=True)
class LocalCutConfig:
    local_radial_radius_factor: float
    local_axial_back_radius_factor: float
    local_axial_forward_radius_factor: float


@dataclass(frozen=True, slots=True)
class PressureCorrectionConfig:
    allow_negative_gauge_pressure: bool


@dataclass(frozen=True, slots=True)
class SurfaceQCConfig:
    require_single_component: bool
    require_watertight: bool
    require_zero_boundary_edges: bool
    require_zero_nonmanifold_edges: bool
    require_zero_self_intersections: bool
    require_zero_degenerate_triangles: bool
    expected_boundary_count: int
    maximum_cap_planarity_error_um: float
    minimum_normal_dot: float
    maximum_extension_length_error_um: float
    maximum_core_surface_distance_um: float
    maximum_core_surface_p95_distance_um: float
    maximum_equivalent_radius_relative_error: float


@dataclass(frozen=True, slots=True)
class OutputsConfig:
    save_stl: bool
    save_vtp: bool
    save_boundary_stl: bool
    save_figures: bool


@dataclass(frozen=True, slots=True)
class SurfacePrepareConfig:
    source_path: Path
    paths: PathsConfig
    geometry: GeometryConfig
    local_cut: LocalCutConfig
    pressure_correction: PressureCorrectionConfig
    qc: SurfaceQCConfig
    verbose: bool
    outputs: OutputsConfig


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a YAML mapping")
    return dict(value)


def _exact(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    mapping = _mapping(value, label)
    unknown = sorted(set(mapping) - keys)
    missing = sorted(keys - set(mapping))
    if unknown:
        raise ValueError(f"Unknown keys in {label}: {', '.join(unknown)}")
    if missing:
        raise ValueError(f"Missing keys in {label}: {', '.join(missing)}")
    return mapping


def _section(root: Mapping[str, Any], key: str, keys: set[str]) -> dict[str, Any]:
    if key not in root:
        raise ValueError(f"Missing YAML section: {key}")
    return _exact(root[key], keys, key)


def _require(value: Any, expected: Any, label: str) -> None:
    if value != expected:
        raise ValueError(f"{label} must be {expected!r}")


def _boolean(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false")
    return value


def _number(value: Any, label: str, *, positive: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if positive and result <= 0:
        raise ValueError(f"{label} must be positive")
    return result


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    return int(value)


def _path(value: Any, label: str, project_root: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty path string")
    path = Path(value.strip()).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def load_surface_prepare_config(
    path: Path, *, project_root: Path
) -> SurfacePrepareConfig:
    """Load every parameter explicitly and reject unsupported geometry modes."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"CFD surface YAML does not exist: {source}")
    root = _exact(
        yaml.safe_load(source.read_text(encoding="utf-8")),
        {
            "schema_version",
            "paths",
            "geometry",
            "local_cut",
            "extension",
            "pressure_correction",
            "qc",
            "runtime",
            "outputs",
        },
        "YAML root",
    )
    if _integer(root["schema_version"], "schema_version") != 1:
        raise ValueError("schema_version must be 1")
    project_root = Path(project_root).resolve()
    paths = _section(root, "paths", {"cfd_preprocess_run", "output_root"})

    geometry = _section(
        root,
        "geometry",
        {
            "input_unit",
            "preserve_original_surface",
            "create_meter_copy",
            "triangulate",
            "modify_only_local_boundary_neighborhoods",
            "allow_global_reconstruction",
        },
    )
    _require(geometry["input_unit"], "um", "geometry.input_unit")
    for key in (
        "preserve_original_surface",
        "triangulate",
        "modify_only_local_boundary_neighborhoods",
    ):
        _require(_boolean(geometry[key], f"geometry.{key}"), True, f"geometry.{key}")
    _require(
        _boolean(geometry["allow_global_reconstruction"], "geometry.allow_global_reconstruction"),
        False,
        "geometry.allow_global_reconstruction",
    )
    create_meter_copy = _boolean(
        geometry["create_meter_copy"], "geometry.create_meter_copy"
    )

    local = _section(
        root,
        "local_cut",
        {
            "enabled",
            "plane_source",
            "local_radial_radius_factor",
            "local_axial_back_radius_factor",
            "local_axial_forward_radius_factor",
            "require_single_cut_loop",
            "fail_if_other_vessel_intersects_local_surgery_region",
        },
    )
    _require(_boolean(local["enabled"], "local_cut.enabled"), True, "local_cut.enabled")
    _require(
        local["plane_source"],
        "cfd_preprocess_boundary_geometry",
        "local_cut.plane_source",
    )
    _require(
        _boolean(local["require_single_cut_loop"], "local_cut.require_single_cut_loop"),
        True,
        "local_cut.require_single_cut_loop",
    )
    _require(
        _boolean(
            local["fail_if_other_vessel_intersects_local_surgery_region"],
            "local_cut.fail_if_other_vessel_intersects_local_surgery_region",
        ),
        True,
        "local_cut.fail_if_other_vessel_intersects_local_surgery_region",
    )

    extension = _section(
        root,
        "extension",
        {"enabled", "use_existing_extension_plan", "cross_section_mode", "cap_distal_end"},
    )
    for key in ("enabled", "use_existing_extension_plan", "cap_distal_end"):
        _require(_boolean(extension[key], f"extension.{key}"), True, f"extension.{key}")
    _require(
        extension["cross_section_mode"],
        "extrude_actual_cut_loop",
        "extension.cross_section_mode",
    )

    pressure = _section(
        root,
        "pressure_correction",
        {"enabled", "outlet_method", "inlet", "outlet", "allow_negative_gauge_pressure"},
    )
    _require(
        _boolean(pressure["enabled"], "pressure_correction.enabled"),
        True,
        "pressure_correction.enabled",
    )
    _require(
        pressure["outlet_method"],
        "poiseuille_equivalent_area",
        "pressure_correction.outlet_method",
    )
    inlet = _exact(pressure["inlet"], {"preserve_flow_rate"}, "pressure_correction.inlet")
    outlet = _exact(
        pressure["outlet"],
        {"preserve_original_plane_pressure_target"},
        "pressure_correction.outlet",
    )
    _require(
        _boolean(inlet["preserve_flow_rate"], "pressure_correction.inlet.preserve_flow_rate"),
        True,
        "pressure_correction.inlet.preserve_flow_rate",
    )
    _require(
        _boolean(
            outlet["preserve_original_plane_pressure_target"],
            "pressure_correction.outlet.preserve_original_plane_pressure_target",
        ),
        True,
        "pressure_correction.outlet.preserve_original_plane_pressure_target",
    )

    qc = _section(
        root,
        "qc",
        {
            "require_single_component",
            "require_watertight",
            "require_zero_boundary_edges",
            "require_zero_nonmanifold_edges",
            "require_zero_self_intersections",
            "require_zero_degenerate_triangles",
            "expected_boundary_count",
            "maximum_cap_planarity_error_um",
            "minimum_normal_dot",
            "maximum_extension_length_error_um",
            "maximum_core_surface_distance_um",
            "maximum_core_surface_p95_distance_um",
            "maximum_equivalent_radius_relative_error",
        },
    )
    qc_booleans = {
        key: _boolean(qc[key], f"qc.{key}")
        for key in (
            "require_single_component",
            "require_watertight",
            "require_zero_boundary_edges",
            "require_zero_nonmanifold_edges",
            "require_zero_self_intersections",
            "require_zero_degenerate_triangles",
        )
    }
    expected_count = _integer(qc["expected_boundary_count"], "qc.expected_boundary_count")
    if expected_count <= 0:
        raise ValueError("qc.expected_boundary_count must be positive")
    minimum_normal_dot = _number(qc["minimum_normal_dot"], "qc.minimum_normal_dot")
    if minimum_normal_dot > 1.0:
        raise ValueError("qc.minimum_normal_dot must not exceed 1")

    runtime = _section(root, "runtime", {"verbose"})
    outputs = _section(
        root,
        "outputs",
        {"save_stl", "save_vtp", "save_boundary_stl", "save_figures"},
    )
    output_config = OutputsConfig(
        **{key: _boolean(outputs[key], f"outputs.{key}") for key in outputs}
    )
    if not all((output_config.save_stl, output_config.save_vtp, output_config.save_boundary_stl)):
        raise ValueError("STL, VTP, and per-boundary STL outputs are mandatory")

    return SurfacePrepareConfig(
        source_path=source,
        paths=PathsConfig(
            cfd_preprocess_run=_path(
                paths["cfd_preprocess_run"], "paths.cfd_preprocess_run", project_root
            ),
            output_root=_path(paths["output_root"], "paths.output_root", project_root),
        ),
        geometry=GeometryConfig(create_meter_copy=create_meter_copy),
        local_cut=LocalCutConfig(
            local_radial_radius_factor=_number(
                local["local_radial_radius_factor"],
                "local_cut.local_radial_radius_factor",
            ),
            local_axial_back_radius_factor=_number(
                local["local_axial_back_radius_factor"],
                "local_cut.local_axial_back_radius_factor",
            ),
            local_axial_forward_radius_factor=_number(
                local["local_axial_forward_radius_factor"],
                "local_cut.local_axial_forward_radius_factor",
            ),
        ),
        pressure_correction=PressureCorrectionConfig(
            allow_negative_gauge_pressure=_boolean(
                pressure["allow_negative_gauge_pressure"],
                "pressure_correction.allow_negative_gauge_pressure",
            )
        ),
        qc=SurfaceQCConfig(
            **qc_booleans,
            expected_boundary_count=expected_count,
            maximum_cap_planarity_error_um=_number(
                qc["maximum_cap_planarity_error_um"],
                "qc.maximum_cap_planarity_error_um",
            ),
            minimum_normal_dot=minimum_normal_dot,
            maximum_extension_length_error_um=_number(
                qc["maximum_extension_length_error_um"],
                "qc.maximum_extension_length_error_um",
            ),
            maximum_core_surface_distance_um=_number(
                qc["maximum_core_surface_distance_um"],
                "qc.maximum_core_surface_distance_um",
            ),
            maximum_core_surface_p95_distance_um=_number(
                qc["maximum_core_surface_p95_distance_um"],
                "qc.maximum_core_surface_p95_distance_um",
            ),
            maximum_equivalent_radius_relative_error=_number(
                qc["maximum_equivalent_radius_relative_error"],
                "qc.maximum_equivalent_radius_relative_error",
            ),
        ),
        verbose=_boolean(runtime["verbose"], "runtime.verbose"),
        outputs=output_config,
    )
