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
class PressureCorrectionConfig:
    allow_negative_gauge_pressure: bool


@dataclass(frozen=True, slots=True)
class RefinementConfig:
    local_mesh_sampling_radius_factor: float
    target_axial_spacing_factor: float
    minimum_ring_count: int
    maximum_ring_count: int


@dataclass(frozen=True, slots=True)
class LoopRegularizationConfig:
    iterations: int
    relaxation: float


@dataclass(frozen=True, slots=True)
class TransitionConfig:
    transition_length_diameters: float
    loop_regularization: LoopRegularizationConfig


@dataclass(frozen=True, slots=True)
class SmoothingConfig:
    iterations: int
    lambda_factor: float
    mu_factor: float
    transition_only: bool


@dataclass(frozen=True, slots=True)
class OptimizationConfig:
    allow_edge_split: bool
    allow_edge_flip: bool


@dataclass(frozen=True, slots=True)
class ExtensionMeshConfig:
    refinement: RefinementConfig
    transition: TransitionConfig
    smoothing: SmoothingConfig
    optimization: OptimizationConfig


@dataclass(frozen=True, slots=True)
class MeshQualityConfig:
    minimum_triangle_angle_deg: float
    maximum_aspect_ratio: float
    maximum_edge_length_to_local_target_ratio: float
    maximum_neighbor_area_ratio: float
    maximum_interface_edge_length_ratio: float
    maximum_bad_triangle_fraction: float


@dataclass(frozen=True, slots=True)
class ManualReviewConfig:
    previous_surface_run: Path
    generate_before_after_figures: bool


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
    backend: BackendConfig
    vmtk: VmtkConfig
    geometry: GeometryConfig
    local_cut: LocalCutConfig
    extension_mesh: ExtensionMeshConfig
    mesh_quality: MeshQualityConfig
    manual_review: ManualReviewConfig
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
            "backend",
            "vmtk",
            "geometry",
            "local_cut",
            "extension",
            "extension_mesh",
            "mesh_quality",
            "manual_review",
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
    backend = _section(root, "backend", {"method"})
    _require(backend["method"], "vmtk_flowextensions", "backend.method")
    vmtk = _section(
        root,
        "vmtk",
        {
            "environment_python",
            "runtime_prefix",
            "official_repository",
            "interpolation_mode",
            "preserve_cross_section_shape",
            "extension_mode",
            "sigma",
            "transition_ratio",
            "adaptive_extension_length",
            "extension_ratio",
            "adaptive_extension_radius",
            "adaptive_boundary_points",
            "postprocess_mode",
            "remesh_after_extension",
            "entity_remesh",
        },
    )
    _require(
        vmtk["interpolation_mode"],
        "thinplatespline",
        "vmtk.interpolation_mode",
    )
    _require(
        _boolean(
            vmtk["preserve_cross_section_shape"],
            "vmtk.preserve_cross_section_shape",
        ),
        False,
        "vmtk.preserve_cross_section_shape",
    )
    extension_mode = vmtk["extension_mode"]
    if extension_mode not in {"centerlinedirection", "boundarynormal"}:
        raise ValueError(
            "vmtk.extension_mode must be 'centerlinedirection' or 'boundarynormal'"
        )
    _require(_number(vmtk["sigma"], "vmtk.sigma"), 1.0, "vmtk.sigma")
    _require(
        _number(vmtk["transition_ratio"], "vmtk.transition_ratio"),
        0.5,
        "vmtk.transition_ratio",
    )
    for key in (
        "adaptive_extension_length",
        "adaptive_extension_radius",
        "adaptive_boundary_points",
    ):
        _require(_boolean(vmtk[key], f"vmtk.{key}"), True, f"vmtk.{key}")
    _require(
        vmtk["postprocess_mode"],
        "cross_seam_active_collar_remesh_then_cap",
        "vmtk.postprocess_mode",
    )
    _require(
        _boolean(vmtk["remesh_after_extension"], "vmtk.remesh_after_extension"),
        True,
        "vmtk.remesh_after_extension",
    )
    _require(
        _number(vmtk["extension_ratio"], "vmtk.extension_ratio"),
        10.0,
        "vmtk.extension_ratio",
    )
    entity_remesh = _exact(
        vmtk["entity_remesh"],
        {
            "enabled",
            "entity_array_name",
            "far_core_entity_id",
            "active_entity_id",
            "exclude_entity_ids",
            "core_collar",
            "element_size_mode",
            "target_edge_length_um",
            "preserve_boundary_edges",
        },
        "vmtk.entity_remesh",
    )
    _require(
        _boolean(entity_remesh["enabled"], "vmtk.entity_remesh.enabled"),
        True,
        "vmtk.entity_remesh.enabled",
    )
    _require(
        entity_remesh["entity_array_name"],
        "RemeshEntityId",
        "vmtk.entity_remesh.entity_array_name",
    )
    far_core_entity_id = _integer(
        entity_remesh["far_core_entity_id"],
        "vmtk.entity_remesh.far_core_entity_id",
    )
    active_entity_id = _integer(
        entity_remesh["active_entity_id"],
        "vmtk.entity_remesh.active_entity_id",
    )
    _require(far_core_entity_id, 1, "vmtk.entity_remesh.far_core_entity_id")
    _require(active_entity_id, 2, "vmtk.entity_remesh.active_entity_id")
    excluded = entity_remesh["exclude_entity_ids"]
    if not isinstance(excluded, list) or any(
        isinstance(value, bool) or not isinstance(value, int) for value in excluded
    ):
        raise ValueError("vmtk.entity_remesh.exclude_entity_ids must be an integer list")
    exclude_entity_ids = tuple(int(value) for value in excluded)
    _require(exclude_entity_ids, (1,), "vmtk.entity_remesh.exclude_entity_ids")
    core_collar = _exact(
        entity_remesh["core_collar"],
        {"mode", "face_layers"},
        "vmtk.entity_remesh.core_collar",
    )
    _require(
        core_collar["mode"],
        "core_face_adjacency_layers",
        "vmtk.entity_remesh.core_collar.mode",
    )
    core_collar_face_layers = _integer(
        core_collar["face_layers"],
        "vmtk.entity_remesh.core_collar.face_layers",
    )
    _require(
        core_collar_face_layers,
        2,
        "vmtk.entity_remesh.core_collar.face_layers",
    )
    _require(
        entity_remesh["element_size_mode"],
        "edgelength",
        "vmtk.entity_remesh.element_size_mode",
    )
    target_edge_length_um = _number(
        entity_remesh["target_edge_length_um"],
        "vmtk.entity_remesh.target_edge_length_um",
    )
    _require(
        target_edge_length_um,
        0.25913916380971913,
        "vmtk.entity_remesh.target_edge_length_um",
    )
    _require(
        _boolean(
            entity_remesh["preserve_boundary_edges"],
            "vmtk.entity_remesh.preserve_boundary_edges",
        ),
        True,
        "vmtk.entity_remesh.preserve_boundary_edges",
    )

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

    extension_mesh = _section(
        root,
        "extension_mesh",
        {"refinement", "transition", "smoothing", "optimization"},
    )
    refinement = _exact(
        extension_mesh["refinement"],
        {
            "enabled",
            "target_edge_length_source",
            "local_mesh_sampling_radius_factor",
            "target_axial_spacing_factor",
            "minimum_ring_count",
            "maximum_ring_count",
            "require_intermediate_rings",
        },
        "extension_mesh.refinement",
    )
    _require(_boolean(refinement["enabled"], "extension_mesh.refinement.enabled"), True, "extension_mesh.refinement.enabled")
    _require(refinement["target_edge_length_source"], "local_original_mesh", "extension_mesh.refinement.target_edge_length_source")
    _require(_boolean(refinement["require_intermediate_rings"], "extension_mesh.refinement.require_intermediate_rings"), True, "extension_mesh.refinement.require_intermediate_rings")
    minimum_rings = _integer(refinement["minimum_ring_count"], "extension_mesh.refinement.minimum_ring_count")
    maximum_rings = _integer(refinement["maximum_ring_count"], "extension_mesh.refinement.maximum_ring_count")
    if minimum_rings < 3 or maximum_rings < minimum_rings:
        raise ValueError("extension ring limits must satisfy 3 <= minimum <= maximum")

    transition = _exact(
        extension_mesh["transition"],
        {"enabled", "transition_length_diameters", "loop_regularization"},
        "extension_mesh.transition",
    )
    _require(_boolean(transition["enabled"], "extension_mesh.transition.enabled"), True, "extension_mesh.transition.enabled")
    regularization = _exact(
        transition["loop_regularization"],
        {"enabled", "iterations", "relaxation", "preserve_area", "preserve_centroid", "preserve_plane"},
        "extension_mesh.transition.loop_regularization",
    )
    for key in ("enabled", "preserve_area", "preserve_centroid", "preserve_plane"):
        _require(_boolean(regularization[key], f"extension_mesh.transition.loop_regularization.{key}"), True, f"extension_mesh.transition.loop_regularization.{key}")
    regularization_iterations = _integer(regularization["iterations"], "extension_mesh.transition.loop_regularization.iterations")
    relaxation = _number(regularization["relaxation"], "extension_mesh.transition.loop_regularization.relaxation")
    if regularization_iterations <= 0 or relaxation >= 1.0:
        raise ValueError("loop regularization requires iterations > 0 and 0 < relaxation < 1")

    smoothing = _exact(
        extension_mesh["smoothing"],
        {"enabled", "method", "iterations", "lambda", "mu", "lock_proximal_ring", "lock_distal_ring", "lock_distal_cap", "lock_original_core", "transition_only"},
        "extension_mesh.smoothing",
    )
    _require(_boolean(smoothing["enabled"], "extension_mesh.smoothing.enabled"), True, "extension_mesh.smoothing.enabled")
    _require(smoothing["method"], "taubin_constrained", "extension_mesh.smoothing.method")
    for key in ("lock_proximal_ring", "lock_distal_ring", "lock_distal_cap", "lock_original_core"):
        _require(_boolean(smoothing[key], f"extension_mesh.smoothing.{key}"), True, f"extension_mesh.smoothing.{key}")
    smoothing_iterations = _integer(smoothing["iterations"], "extension_mesh.smoothing.iterations")
    lambda_factor = _number(smoothing["lambda"], "extension_mesh.smoothing.lambda")
    mu_factor = _number(smoothing["mu"], "extension_mesh.smoothing.mu", positive=False)
    if smoothing_iterations <= 0 or mu_factor >= 0:
        raise ValueError("constrained Taubin smoothing requires iterations > 0 and mu < 0")

    optimization = _exact(
        extension_mesh["optimization"],
        {"enabled", "improve_triangle_diagonals", "allow_edge_collapse", "allow_edge_split", "allow_edge_flip"},
        "extension_mesh.optimization",
    )
    for key in ("enabled", "improve_triangle_diagonals"):
        _require(_boolean(optimization[key], f"extension_mesh.optimization.{key}"), True, f"extension_mesh.optimization.{key}")
    _require(_boolean(optimization["allow_edge_collapse"], "extension_mesh.optimization.allow_edge_collapse"), False, "extension_mesh.optimization.allow_edge_collapse")

    mesh_quality = _section(
        root,
        "mesh_quality",
        {"enabled", "minimum_triangle_angle_deg", "maximum_aspect_ratio", "maximum_edge_length_to_local_target_ratio", "maximum_neighbor_area_ratio", "maximum_interface_edge_length_ratio", "maximum_bad_triangle_fraction"},
    )
    _require(_boolean(mesh_quality["enabled"], "mesh_quality.enabled"), True, "mesh_quality.enabled")
    bad_fraction = _number(mesh_quality["maximum_bad_triangle_fraction"], "mesh_quality.maximum_bad_triangle_fraction", positive=False)
    if not 0 <= bad_fraction <= 1:
        raise ValueError("mesh_quality.maximum_bad_triangle_fraction must be in [0, 1]")

    manual_review = _section(root, "manual_review", {"previous_surface_run", "generate_before_after_figures"})
    _require(_boolean(manual_review["generate_before_after_figures"], "manual_review.generate_before_after_figures"), True, "manual_review.generate_before_after_figures")

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
        "geometry_cross_sections_20_station_poiseuille",
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
        backend=BackendConfig(method="vmtk_flowextensions"),
        vmtk=VmtkConfig(
            environment_python=_path(
                vmtk["environment_python"],
                "vmtk.environment_python",
                project_root,
            ),
            runtime_prefix=_path(
                vmtk["runtime_prefix"],
                "vmtk.runtime_prefix",
                project_root,
            ),
            official_repository=_path(
                vmtk["official_repository"],
                "vmtk.official_repository",
                project_root,
            ),
            interpolation_mode="thinplatespline",
            preserve_cross_section_shape=False,
            extension_mode=extension_mode,
            sigma=1.0,
            transition_ratio=0.5,
            adaptive_extension_length=True,
            extension_ratio=10.0,
            adaptive_extension_radius=True,
            adaptive_boundary_points=True,
            postprocess_mode="cross_seam_active_collar_remesh_then_cap",
            remesh_after_extension=True,
            entity_remesh=EntityRemeshConfig(
                enabled=True,
                entity_array_name="RemeshEntityId",
                far_core_entity_id=far_core_entity_id,
                active_entity_id=active_entity_id,
                exclude_entity_ids=exclude_entity_ids,
                core_collar=CoreCollarConfig(
                    mode="core_face_adjacency_layers",
                    face_layers=core_collar_face_layers,
                ),
                element_size_mode="edgelength",
                target_edge_length_um=target_edge_length_um,
                preserve_boundary_edges=True,
            ),
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
        extension_mesh=ExtensionMeshConfig(
            refinement=RefinementConfig(
                local_mesh_sampling_radius_factor=_number(refinement["local_mesh_sampling_radius_factor"], "extension_mesh.refinement.local_mesh_sampling_radius_factor"),
                target_axial_spacing_factor=_number(refinement["target_axial_spacing_factor"], "extension_mesh.refinement.target_axial_spacing_factor"),
                minimum_ring_count=minimum_rings,
                maximum_ring_count=maximum_rings,
            ),
            transition=TransitionConfig(
                transition_length_diameters=_number(transition["transition_length_diameters"], "extension_mesh.transition.transition_length_diameters"),
                loop_regularization=LoopRegularizationConfig(iterations=regularization_iterations, relaxation=relaxation),
            ),
            smoothing=SmoothingConfig(
                iterations=smoothing_iterations,
                lambda_factor=lambda_factor,
                mu_factor=mu_factor,
                transition_only=_boolean(smoothing["transition_only"], "extension_mesh.smoothing.transition_only"),
            ),
            optimization=OptimizationConfig(
                allow_edge_split=_boolean(optimization["allow_edge_split"], "extension_mesh.optimization.allow_edge_split"),
                allow_edge_flip=_boolean(optimization["allow_edge_flip"], "extension_mesh.optimization.allow_edge_flip"),
            ),
        ),
        mesh_quality=MeshQualityConfig(
            minimum_triangle_angle_deg=_number(mesh_quality["minimum_triangle_angle_deg"], "mesh_quality.minimum_triangle_angle_deg"),
            maximum_aspect_ratio=_number(mesh_quality["maximum_aspect_ratio"], "mesh_quality.maximum_aspect_ratio"),
            maximum_edge_length_to_local_target_ratio=_number(mesh_quality["maximum_edge_length_to_local_target_ratio"], "mesh_quality.maximum_edge_length_to_local_target_ratio"),
            maximum_neighbor_area_ratio=_number(mesh_quality["maximum_neighbor_area_ratio"], "mesh_quality.maximum_neighbor_area_ratio"),
            maximum_interface_edge_length_ratio=_number(mesh_quality["maximum_interface_edge_length_ratio"], "mesh_quality.maximum_interface_edge_length_ratio"),
            maximum_bad_triangle_fraction=bad_fraction,
        ),
        manual_review=ManualReviewConfig(
            previous_surface_run=_path(manual_review["previous_surface_run"], "manual_review.previous_surface_run", project_root),
            generate_before_after_figures=True,
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
