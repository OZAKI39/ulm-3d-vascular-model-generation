"""Strict configuration for the formal CFD surface-production path."""

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
class BackendConfig:
    method: str


@dataclass(frozen=True, slots=True)
class CoreCollarConfig:
    mode: str
    face_layers: int


@dataclass(frozen=True, slots=True)
class EntityRemeshConfig:
    enabled: bool
    entity_array_name: str
    far_core_entity_id: int
    active_entity_id: int
    exclude_entity_ids: tuple[int, ...]
    core_collar: CoreCollarConfig
    element_size_mode: str
    target_edge_length_um: float
    preserve_boundary_edges: bool


@dataclass(frozen=True, slots=True)
class VmtkConfig:
    environment_python: Path
    runtime_prefix: Path
    official_repository: Path
    interpolation_mode: str
    preserve_cross_section_shape: bool
    extension_mode: str
    sigma: float
    transition_ratio: float
    adaptive_extension_length: bool
    extension_ratio: float
    adaptive_extension_radius: bool
    adaptive_boundary_points: bool
    postprocess_mode: str
    remesh_after_extension: bool
    entity_remesh: EntityRemeshConfig


@dataclass(frozen=True, slots=True)
class GeometryConfig:
    create_meter_copy: bool


@dataclass(frozen=True, slots=True)
class LocalCutConfig:
    local_radial_radius_factor: float
    local_axial_back_radius_factor: float
    local_axial_forward_radius_factor: float


@dataclass(frozen=True, slots=True)
class RefinementConfig:
    local_mesh_sampling_radius_factor: float


@dataclass(frozen=True, slots=True)
class ExtensionMeshConfig:
    refinement: RefinementConfig


@dataclass(frozen=True, slots=True)
class MeshQualityConfig:
    minimum_triangle_angle_deg: float
    maximum_aspect_ratio: float
    maximum_edge_length_to_local_target_ratio: float
    maximum_neighbor_area_ratio: float
    maximum_interface_edge_length_ratio: float
    maximum_bad_triangle_fraction: float


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
class SurfacePrepareConfig:
    source_path: Path
    paths: PathsConfig
    backend: BackendConfig
    vmtk: VmtkConfig
    geometry: GeometryConfig
    local_cut: LocalCutConfig
    extension_mesh: ExtensionMeshConfig
    mesh_quality: MeshQualityConfig
    pressure_correction: PressureCorrectionConfig
    qc: SurfaceQCConfig


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
    """Load the one supported VMTK boundary-normal production configuration."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"CFD surface YAML does not exist: {source}")
    root = _exact(
        yaml.safe_load(source.read_text(encoding="utf-8")),
        {
            "schema_version", "paths", "backend", "vmtk", "geometry",
            "local_cut", "extension_mesh", "mesh_quality",
            "pressure_correction", "qc",
        },
        "YAML root",
    )
    _require(_integer(root["schema_version"], "schema_version"), 1, "schema_version")
    project_root = Path(project_root).resolve()
    paths = _section(root, "paths", {"cfd_preprocess_run", "output_root"})
    backend = _section(root, "backend", {"method"})
    _require(backend["method"], "vmtk_flowextensions", "backend.method")

    vmtk = _section(
        root, "vmtk",
        {
            "environment_python", "runtime_prefix", "official_repository",
            "interpolation_mode", "preserve_cross_section_shape",
            "extension_mode", "sigma", "transition_ratio",
            "adaptive_extension_length", "extension_ratio",
            "adaptive_extension_radius", "adaptive_boundary_points",
            "postprocess_mode", "remesh_after_extension", "entity_remesh",
        },
    )
    fixed_vmtk = {
        "interpolation_mode": "thinplatespline",
        "preserve_cross_section_shape": False,
        "extension_mode": "boundarynormal",
        "sigma": 1.0,
        "transition_ratio": 0.5,
        "adaptive_extension_length": True,
        "extension_ratio": 10.0,
        "adaptive_extension_radius": True,
        "adaptive_boundary_points": True,
        "postprocess_mode": "cross_seam_active_collar_remesh_then_cap",
        "remesh_after_extension": True,
    }
    for key, expected in fixed_vmtk.items():
        value = vmtk[key]
        if isinstance(expected, bool):
            value = _boolean(value, f"vmtk.{key}")
        elif isinstance(expected, float):
            value = _number(value, f"vmtk.{key}")
        _require(value, expected, f"vmtk.{key}")

    entity = _exact(
        vmtk["entity_remesh"],
        {
            "enabled", "entity_array_name", "far_core_entity_id",
            "active_entity_id", "exclude_entity_ids", "core_collar",
            "element_size_mode", "target_edge_length_um",
            "preserve_boundary_edges",
        },
        "vmtk.entity_remesh",
    )
    collar = _exact(
        entity["core_collar"], {"mode", "face_layers"},
        "vmtk.entity_remesh.core_collar",
    )
    fixed_entity = {
        "enabled": True,
        "entity_array_name": "RemeshEntityId",
        "far_core_entity_id": 1,
        "active_entity_id": 2,
        "exclude_entity_ids": [1],
        "element_size_mode": "edgelength",
        "target_edge_length_um": 0.25913916380971913,
        "preserve_boundary_edges": True,
    }
    for key, expected in fixed_entity.items():
        _require(entity[key], expected, f"vmtk.entity_remesh.{key}")
    _require(collar["mode"], "core_face_adjacency_layers", "vmtk.entity_remesh.core_collar.mode")
    _require(collar["face_layers"], 2, "vmtk.entity_remesh.core_collar.face_layers")

    geometry = _section(
        root, "geometry",
        {
            "input_unit", "preserve_original_surface", "create_meter_copy",
            "triangulate", "modify_only_local_boundary_neighborhoods",
            "allow_global_reconstruction",
        },
    )
    _require(geometry["input_unit"], "um", "geometry.input_unit")
    for key in ("preserve_original_surface", "triangulate", "modify_only_local_boundary_neighborhoods"):
        _require(_boolean(geometry[key], f"geometry.{key}"), True, f"geometry.{key}")
    _require(_boolean(geometry["allow_global_reconstruction"], "geometry.allow_global_reconstruction"), False, "geometry.allow_global_reconstruction")

    local = _section(
        root, "local_cut",
        {
            "enabled", "plane_source", "local_radial_radius_factor",
            "local_axial_back_radius_factor", "local_axial_forward_radius_factor",
            "require_single_cut_loop",
            "fail_if_other_vessel_intersects_local_surgery_region",
        },
    )
    _require(_boolean(local["enabled"], "local_cut.enabled"), True, "local_cut.enabled")
    _require(local["plane_source"], "cfd_preprocess_boundary_geometry", "local_cut.plane_source")
    _require(_boolean(local["require_single_cut_loop"], "local_cut.require_single_cut_loop"), True, "local_cut.require_single_cut_loop")
    _require(_boolean(local["fail_if_other_vessel_intersects_local_surgery_region"], "local_cut.fail_if_other_vessel_intersects_local_surgery_region"), True, "local_cut.fail_if_other_vessel_intersects_local_surgery_region")

    extension_mesh = _section(root, "extension_mesh", {"refinement"})
    refinement = _exact(
        extension_mesh["refinement"], {"local_mesh_sampling_radius_factor"},
        "extension_mesh.refinement",
    )
    mesh = _section(
        root, "mesh_quality",
        {
            "enabled", "minimum_triangle_angle_deg", "maximum_aspect_ratio",
            "maximum_edge_length_to_local_target_ratio",
            "maximum_neighbor_area_ratio", "maximum_interface_edge_length_ratio",
            "maximum_bad_triangle_fraction",
        },
    )
    _require(_boolean(mesh["enabled"], "mesh_quality.enabled"), True, "mesh_quality.enabled")
    bad_fraction = _number(mesh["maximum_bad_triangle_fraction"], "mesh_quality.maximum_bad_triangle_fraction", positive=False)
    if not 0 <= bad_fraction <= 1:
        raise ValueError("mesh_quality.maximum_bad_triangle_fraction must be in [0, 1]")

    pressure = _section(
        root, "pressure_correction",
        {"enabled", "outlet_method", "inlet", "outlet", "allow_negative_gauge_pressure"},
    )
    _require(_boolean(pressure["enabled"], "pressure_correction.enabled"), True, "pressure_correction.enabled")
    _require(pressure["outlet_method"], "geometry_cross_sections_20_station_poiseuille", "pressure_correction.outlet_method")
    inlet = _exact(pressure["inlet"], {"preserve_flow_rate"}, "pressure_correction.inlet")
    outlet = _exact(pressure["outlet"], {"preserve_original_plane_pressure_target"}, "pressure_correction.outlet")
    _require(_boolean(inlet["preserve_flow_rate"], "pressure_correction.inlet.preserve_flow_rate"), True, "pressure_correction.inlet.preserve_flow_rate")
    _require(_boolean(outlet["preserve_original_plane_pressure_target"], "pressure_correction.outlet.preserve_original_plane_pressure_target"), True, "pressure_correction.outlet.preserve_original_plane_pressure_target")

    qc = _section(
        root, "qc",
        {
            "require_single_component", "require_watertight",
            "require_zero_boundary_edges", "require_zero_nonmanifold_edges",
            "require_zero_self_intersections", "require_zero_degenerate_triangles",
            "expected_boundary_count", "maximum_cap_planarity_error_um",
            "minimum_normal_dot", "maximum_extension_length_error_um",
            "maximum_core_surface_distance_um",
            "maximum_core_surface_p95_distance_um",
            "maximum_equivalent_radius_relative_error",
        },
    )
    qc_booleans = {
        key: _boolean(qc[key], f"qc.{key}")
        for key in (
            "require_single_component", "require_watertight",
            "require_zero_boundary_edges", "require_zero_nonmanifold_edges",
            "require_zero_self_intersections", "require_zero_degenerate_triangles",
        )
    }
    expected_count = _integer(qc["expected_boundary_count"], "qc.expected_boundary_count")
    minimum_normal_dot = _number(qc["minimum_normal_dot"], "qc.minimum_normal_dot")
    if expected_count <= 0 or minimum_normal_dot > 1:
        raise ValueError("QC boundary count and normal threshold are invalid")

    return SurfacePrepareConfig(
        source_path=source,
        paths=PathsConfig(
            _path(paths["cfd_preprocess_run"], "paths.cfd_preprocess_run", project_root),
            _path(paths["output_root"], "paths.output_root", project_root),
        ),
        backend=BackendConfig("vmtk_flowextensions"),
        vmtk=VmtkConfig(
            environment_python=_path(vmtk["environment_python"], "vmtk.environment_python", project_root),
            runtime_prefix=_path(vmtk["runtime_prefix"], "vmtk.runtime_prefix", project_root),
            official_repository=_path(vmtk["official_repository"], "vmtk.official_repository", project_root),
            interpolation_mode="thinplatespline",
            preserve_cross_section_shape=False,
            extension_mode="boundarynormal",
            sigma=1.0,
            transition_ratio=0.5,
            adaptive_extension_length=True,
            extension_ratio=10.0,
            adaptive_extension_radius=True,
            adaptive_boundary_points=True,
            postprocess_mode="cross_seam_active_collar_remesh_then_cap",
            remesh_after_extension=True,
            entity_remesh=EntityRemeshConfig(
                True, "RemeshEntityId", 1, 2, (1,),
                CoreCollarConfig("core_face_adjacency_layers", 2),
                "edgelength", 0.25913916380971913, True,
            ),
        ),
        geometry=GeometryConfig(_boolean(geometry["create_meter_copy"], "geometry.create_meter_copy")),
        local_cut=LocalCutConfig(
            _number(local["local_radial_radius_factor"], "local_cut.local_radial_radius_factor"),
            _number(local["local_axial_back_radius_factor"], "local_cut.local_axial_back_radius_factor"),
            _number(local["local_axial_forward_radius_factor"], "local_cut.local_axial_forward_radius_factor"),
        ),
        extension_mesh=ExtensionMeshConfig(
            RefinementConfig(_number(refinement["local_mesh_sampling_radius_factor"], "extension_mesh.refinement.local_mesh_sampling_radius_factor"))
        ),
        mesh_quality=MeshQualityConfig(
            _number(mesh["minimum_triangle_angle_deg"], "mesh_quality.minimum_triangle_angle_deg"),
            _number(mesh["maximum_aspect_ratio"], "mesh_quality.maximum_aspect_ratio"),
            _number(mesh["maximum_edge_length_to_local_target_ratio"], "mesh_quality.maximum_edge_length_to_local_target_ratio"),
            _number(mesh["maximum_neighbor_area_ratio"], "mesh_quality.maximum_neighbor_area_ratio"),
            _number(mesh["maximum_interface_edge_length_ratio"], "mesh_quality.maximum_interface_edge_length_ratio"),
            bad_fraction,
        ),
        pressure_correction=PressureCorrectionConfig(_boolean(pressure["allow_negative_gauge_pressure"], "pressure_correction.allow_negative_gauge_pressure")),
        qc=SurfaceQCConfig(
            **qc_booleans,
            expected_boundary_count=expected_count,
            maximum_cap_planarity_error_um=_number(qc["maximum_cap_planarity_error_um"], "qc.maximum_cap_planarity_error_um"),
            minimum_normal_dot=minimum_normal_dot,
            maximum_extension_length_error_um=_number(qc["maximum_extension_length_error_um"], "qc.maximum_extension_length_error_um"),
            maximum_core_surface_distance_um=_number(qc["maximum_core_surface_distance_um"], "qc.maximum_core_surface_distance_um"),
            maximum_core_surface_p95_distance_um=_number(qc["maximum_core_surface_p95_distance_um"], "qc.maximum_core_surface_p95_distance_um"),
            maximum_equivalent_radius_relative_error=_number(qc["maximum_equivalent_radius_relative_error"], "qc.maximum_equivalent_radius_relative_error"),
        ),
    )
