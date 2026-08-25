"""Strict YAML configuration for CFD boundary preprocessing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml


@dataclass(frozen=True, slots=True)
class PathsConfig:
    sampling_run: Path
    rodent_run: Path | None
    model_run: Path | None
    model_output_root: Path
    output_root: Path


@dataclass(frozen=True, slots=True)
class SelectionConfig:
    roi_anchor: int | None
    roi_id: str | None


@dataclass(frozen=True, slots=True)
class FluidConfig:
    density_kg_m3: float
    kinematic_viscosity_m2_s: float

    @property
    def dynamic_viscosity_pa_s(self) -> float:
        return self.density_kg_m3 * self.kinematic_viscosity_m2_s


@dataclass(frozen=True, slots=True)
class InletConfig:
    boundary_type: str
    mean_velocity_mm_s: float | None
    flow_rate_m3_s: float | None


@dataclass(frozen=True, slots=True)
class SolverConfig:
    relative_mass_tolerance: float
    relative_node_residual_tolerance: float
    reverse_flow_tolerance_m3_s: float


@dataclass(frozen=True, slots=True)
class Global1DConfig:
    inlet: InletConfig
    leaf_pressure_pa: float
    solver: SolverConfig


@dataclass(frozen=True, slots=True)
class ReadinessConfig:
    minimum_cut_ports: int
    require_zero_true_terminals: bool
    required_assumed_inlet_count: int
    minimum_assumed_outlet_count: int
    require_connected_roi: bool
    require_cycle_rank_zero: bool


@dataclass(frozen=True, slots=True)
class TransferConfig:
    maximum_cut_position_error_um: float
    maximum_cut_radius_relative_error: float
    relative_port_mass_tolerance: float


@dataclass(frozen=True, slots=True)
class ExtensionConfig:
    enabled: bool
    inlet_length_diameters: float
    outlet_length_diameters: float


@dataclass(frozen=True, slots=True)
class OutputConfig:
    save_global_vtp: bool
    save_port_vtp: bool
    save_figures: bool


@dataclass(frozen=True, slots=True)
class CFDPreprocessConfig:
    source_path: Path
    paths: PathsConfig
    selection: SelectionConfig
    fluid: FluidConfig
    global_1d: Global1DConfig
    readiness: ReadinessConfig
    transfer: TransferConfig
    extension: ExtensionConfig
    verbose: bool
    outputs: OutputConfig


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a YAML mapping")
    return dict(value)


def _exact(mapping: Mapping[str, Any], keys: set[str], label: str) -> dict[str, Any]:
    values = dict(mapping)
    unknown = sorted(set(values) - keys)
    missing = sorted(keys - set(values))
    if unknown:
        raise ValueError(f"Unknown keys in {label}: {', '.join(unknown)}")
    if missing:
        raise ValueError(f"Missing keys in {label}: {', '.join(missing)}")
    return values


def _section(root: Mapping[str, Any], key: str, keys: set[str]) -> dict[str, Any]:
    if key not in root:
        raise ValueError(f"Missing YAML section: {key}")
    return _exact(_mapping(root[key], key), keys, key)


def _bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false")
    return value


def _number(
    value: Any,
    label: str,
    *,
    positive: bool = False,
    nullable: bool = False,
) -> float | None:
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        suffix = " or null" if nullable else ""
        raise ValueError(f"{label} must be numeric{suffix}")
    result = float(value)
    if positive and result <= 0:
        raise ValueError(f"{label} must be positive")
    return result


def _integer(value: Any, label: str, *, nullable: bool = False) -> int | None:
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        suffix = " or null" if nullable else ""
        raise ValueError(f"{label} must be an integer{suffix}")
    return int(value)


def _string(value: Any, label: str, *, nullable: bool = False) -> str | None:
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not value.strip():
        suffix = " or null" if nullable else ""
        raise ValueError(f"{label} must be a non-empty string{suffix}")
    return value.strip()


def _path(value: Any, label: str, root: Path, *, nullable: bool = False) -> Path | None:
    text = _string(value, label, nullable=nullable)
    if text is None:
        return None
    candidate = Path(text).expanduser()
    return (
        candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()
    )


def _require(value: Any, expected: Any, label: str) -> None:
    if value != expected:
        raise ValueError(f"{label} must be {expected!r}")


def load_cfd_preprocess_config(
    path: Path, *, project_root: Path
) -> CFDPreprocessConfig:
    """Load a complete configuration and reject all silent scientific fallbacks."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"CFD preprocessing YAML does not exist: {source}")
    root = _mapping(yaml.safe_load(source.read_text(encoding="utf-8")), "YAML root")
    top_keys = {
        "schema_version",
        "paths",
        "selection",
        "units",
        "simulation_direction",
        "fluid",
        "global_1d",
        "roi_readiness",
        "boundary_transfer",
        "port_geometry",
        "physics_baseline",
        "runtime",
        "outputs",
    }
    root = _exact(root, top_keys, "YAML root")
    if _integer(root["schema_version"], "schema_version") != 1:
        raise ValueError("schema_version must be 1")
    project_root = Path(project_root).resolve()

    paths = _section(
        root,
        "paths",
        {
            "sampling_run",
            "rodent_run",
            "model_run",
            "model_output_root",
            "output_root",
        },
    )
    sampling_run = _path(paths["sampling_run"], "paths.sampling_run", project_root)
    rodent_run = _path(
        paths["rodent_run"], "paths.rodent_run", project_root, nullable=True
    )
    model_run = _path(
        paths["model_run"], "paths.model_run", project_root, nullable=True
    )
    model_output_root = _path(
        paths["model_output_root"], "paths.model_output_root", project_root
    )
    output_root = _path(paths["output_root"], "paths.output_root", project_root)
    assert sampling_run and model_output_root and output_root

    selection = _section(root, "selection", {"roi_anchor", "roi_id"})
    roi_anchor = _integer(
        selection["roi_anchor"], "selection.roi_anchor", nullable=True
    )
    roi_id = _string(selection["roi_id"], "selection.roi_id", nullable=True)
    if (roi_anchor is None) == (roi_id is None):
        raise ValueError("Exactly one ROI selector must be non-null")
    if roi_anchor is not None and roi_anchor < 0:
        raise ValueError("selection.roi_anchor must be non-negative")

    units = _section(
        root,
        "units",
        {"coordinate_unit", "swc_radius_unit", "pressure_unit", "flow_rate_unit"},
    )
    for key, expected in {
        "coordinate_unit": "um",
        "swc_radius_unit": "um",
        "pressure_unit": "Pa",
        "flow_rate_unit": "m3/s",
    }.items():
        _require(units[key], expected, f"units.{key}")

    direction = _section(
        root, "simulation_direction", {"mode", "measured_flow_direction"}
    )
    _require(direction["mode"], "swc_parent_to_current", "simulation_direction.mode")
    _require(
        _bool(
            direction["measured_flow_direction"],
            "simulation_direction.measured_flow_direction",
        ),
        False,
        "simulation_direction.measured_flow_direction",
    )

    fluid = _section(
        root,
        "fluid",
        {"model", "density_kg_m3", "kinematic_viscosity_m2_s"},
    )
    _require(fluid["model"], "newtonian", "fluid.model")
    density = _number(fluid["density_kg_m3"], "fluid.density_kg_m3", positive=True)
    nu = _number(
        fluid["kinematic_viscosity_m2_s"],
        "fluid.kinematic_viscosity_m2_s",
        positive=True,
    )
    assert density is not None and nu is not None

    global_1d = _section(
        root,
        "global_1d",
        {"enabled", "network", "inlet", "outlets", "resistance", "solver"},
    )
    _require(
        _bool(global_1d["enabled"], "global_1d.enabled"), True, "global_1d.enabled"
    )
    network = _exact(
        _mapping(global_1d["network"], "global_1d.network"),
        {"use_source_swc_edges"},
        "global_1d.network",
    )
    _require(
        _bool(
            network["use_source_swc_edges"], "global_1d.network.use_source_swc_edges"
        ),
        True,
        "global_1d.network.use_source_swc_edges",
    )
    inlet = _exact(
        _mapping(global_1d["inlet"], "global_1d.inlet"),
        {"root_mode", "boundary_type", "mean_velocity_mm_s", "flow_rate_m3_s"},
        "global_1d.inlet",
    )
    _require(inlet["root_mode"], "single_structural_root", "global_1d.inlet.root_mode")
    boundary_type = _string(inlet["boundary_type"], "global_1d.inlet.boundary_type")
    mean_velocity = _number(
        inlet["mean_velocity_mm_s"],
        "global_1d.inlet.mean_velocity_mm_s",
        positive=True,
        nullable=True,
    )
    flow_rate = _number(
        inlet["flow_rate_m3_s"],
        "global_1d.inlet.flow_rate_m3_s",
        positive=True,
        nullable=True,
    )
    if boundary_type == "mean_velocity":
        if mean_velocity is None or flow_rate is not None:
            raise ValueError("mean_velocity mode requires only mean_velocity_mm_s")
    elif boundary_type == "flow_rate":
        if flow_rate is None or mean_velocity is not None:
            raise ValueError("flow_rate mode requires only flow_rate_m3_s")
    else:
        raise ValueError(
            "global_1d.inlet.boundary_type must be mean_velocity or flow_rate"
        )
    outlets = _exact(
        _mapping(global_1d["outlets"], "global_1d.outlets"),
        {"leaf_pressure_pa"},
        "global_1d.outlets",
    )
    leaf_pressure = _number(
        outlets["leaf_pressure_pa"], "global_1d.outlets.leaf_pressure_pa"
    )
    resistance = _exact(
        _mapping(global_1d["resistance"], "global_1d.resistance"),
        {"radius_interpolation", "integration"},
        "global_1d.resistance",
    )
    _require(
        resistance["radius_interpolation"], "linear", "resistance.radius_interpolation"
    )
    _require(resistance["integration"], "exact_linear_radius", "resistance.integration")
    solver = _exact(
        _mapping(global_1d["solver"], "global_1d.solver"),
        {
            "backend",
            "relative_mass_tolerance",
            "relative_node_residual_tolerance",
            "reverse_flow_tolerance_m3_s",
        },
        "global_1d.solver",
    )
    _require(solver["backend"], "scipy_sparse", "global_1d.solver.backend")
    mass_tol = _number(
        solver["relative_mass_tolerance"],
        "global_1d.solver.relative_mass_tolerance",
        positive=True,
    )
    residual_tol = _number(
        solver["relative_node_residual_tolerance"],
        "global_1d.solver.relative_node_residual_tolerance",
        positive=True,
    )
    reverse_tol = _number(
        solver["reverse_flow_tolerance_m3_s"],
        "global_1d.solver.reverse_flow_tolerance_m3_s",
        positive=True,
    )
    assert leaf_pressure is not None and mass_tol and residual_tol and reverse_tol

    ready = _section(
        root,
        "roi_readiness",
        {
            "minimum_cut_ports",
            "require_zero_true_terminals",
            "required_assumed_inlet_count",
            "minimum_assumed_outlet_count",
            "require_connected_roi",
            "require_cycle_rank_zero",
        },
    )
    min_ports = _integer(ready["minimum_cut_ports"], "roi_readiness.minimum_cut_ports")
    inlet_count = _integer(
        ready["required_assumed_inlet_count"],
        "roi_readiness.required_assumed_inlet_count",
    )
    outlet_count = _integer(
        ready["minimum_assumed_outlet_count"],
        "roi_readiness.minimum_assumed_outlet_count",
    )
    if (
        min_ports is None
        or min_ports < 1
        or inlet_count != 1
        or outlet_count is None
        or outlet_count < 1
    ):
        raise ValueError(
            "ROI readiness port counts must define one inlet and at least one outlet"
        )

    transfer_root = _section(
        root,
        "boundary_transfer",
        {
            "inlet",
            "outlet",
            "pressure_reference",
            "maximum_cut_position_error_um",
            "maximum_cut_radius_relative_error",
            "relative_port_mass_tolerance",
        },
    )
    transfer_inlet = _exact(
        _mapping(transfer_root["inlet"], "boundary_transfer.inlet"),
        {"type", "velocity_profile"},
        "boundary_transfer.inlet",
    )
    _require(transfer_inlet["type"], "flow_rate_1d", "boundary_transfer.inlet.type")
    _require(
        transfer_inlet["velocity_profile"],
        "parabolic",
        "boundary_transfer.inlet.velocity_profile",
    )
    transfer_outlet = _exact(
        _mapping(transfer_root["outlet"], "boundary_transfer.outlet"),
        {"type"},
        "boundary_transfer.outlet",
    )
    _require(
        transfer_outlet["type"], "direct_1d_pressure", "boundary_transfer.outlet.type"
    )
    pressure_ref = _exact(
        _mapping(
            transfer_root["pressure_reference"], "boundary_transfer.pressure_reference"
        ),
        {"type"},
        "boundary_transfer.pressure_reference",
    )
    _require(
        pressure_ref["type"],
        "global_structural_leaves_zero_gauge",
        "boundary_transfer.pressure_reference.type",
    )
    position_tol = _number(
        transfer_root["maximum_cut_position_error_um"],
        "boundary_transfer.maximum_cut_position_error_um",
        positive=True,
    )
    radius_tol = _number(
        transfer_root["maximum_cut_radius_relative_error"],
        "boundary_transfer.maximum_cut_radius_relative_error",
        positive=True,
    )
    port_mass_tol = _number(
        transfer_root["relative_port_mass_tolerance"],
        "boundary_transfer.relative_port_mass_tolerance",
        positive=True,
    )
    assert position_tol and radius_tol and port_mass_tol

    geometry = _section(
        root, "port_geometry", {"modify_surface", "plane", "extension_plan"}
    )
    _require(
        _bool(geometry["modify_surface"], "port_geometry.modify_surface"),
        False,
        "port_geometry.modify_surface",
    )
    plane = _exact(
        _mapping(geometry["plane"], "port_geometry.plane"),
        {"normal_source"},
        "port_geometry.plane",
    )
    _require(
        plane["normal_source"],
        "local_centerline_tangent",
        "port_geometry.plane.normal_source",
    )
    extension = _exact(
        _mapping(geometry["extension_plan"], "port_geometry.extension_plan"),
        {"enabled", "inlet_length_diameters", "outlet_length_diameters"},
        "port_geometry.extension_plan",
    )
    extension_enabled = _bool(
        extension["enabled"], "port_geometry.extension_plan.enabled"
    )
    inlet_d = _number(
        extension["inlet_length_diameters"],
        "port_geometry.extension_plan.inlet_length_diameters",
        positive=True,
    )
    outlet_d = _number(
        extension["outlet_length_diameters"],
        "port_geometry.extension_plan.outlet_length_diameters",
        positive=True,
    )
    assert inlet_d and outlet_d

    physics = _section(
        root,
        "physics_baseline",
        {
            "steady",
            "incompressible",
            "rigid_wall",
            "no_slip_wall",
            "turbulence_model",
            "pulsatility",
        },
    )
    for key in ("steady", "incompressible", "rigid_wall", "no_slip_wall"):
        _require(
            _bool(physics[key], f"physics_baseline.{key}"),
            True,
            f"physics_baseline.{key}",
        )
    _require(physics["turbulence_model"], "none", "physics_baseline.turbulence_model")
    _require(
        _bool(physics["pulsatility"], "physics_baseline.pulsatility"),
        False,
        "physics_baseline.pulsatility",
    )
    runtime = _section(root, "runtime", {"verbose"})
    output_values = _section(
        root,
        "outputs",
        {"save_global_vtp", "save_port_vtp", "save_figures"},
    )

    return CFDPreprocessConfig(
        source_path=source,
        paths=PathsConfig(
            sampling_run=sampling_run,
            rodent_run=rodent_run,
            model_run=model_run,
            model_output_root=model_output_root,
            output_root=output_root,
        ),
        selection=SelectionConfig(roi_anchor, roi_id),
        fluid=FluidConfig(density, nu),
        global_1d=Global1DConfig(
            inlet=InletConfig(str(boundary_type), mean_velocity, flow_rate),
            leaf_pressure_pa=leaf_pressure,
            solver=SolverConfig(mass_tol, residual_tol, reverse_tol),
        ),
        readiness=ReadinessConfig(
            minimum_cut_ports=min_ports,
            require_zero_true_terminals=_bool(
                ready["require_zero_true_terminals"],
                "roi_readiness.require_zero_true_terminals",
            ),
            required_assumed_inlet_count=inlet_count,
            minimum_assumed_outlet_count=outlet_count,
            require_connected_roi=_bool(
                ready["require_connected_roi"], "roi_readiness.require_connected_roi"
            ),
            require_cycle_rank_zero=_bool(
                ready["require_cycle_rank_zero"],
                "roi_readiness.require_cycle_rank_zero",
            ),
        ),
        transfer=TransferConfig(position_tol, radius_tol, port_mass_tol),
        extension=ExtensionConfig(extension_enabled, inlet_d, outlet_d),
        verbose=_bool(runtime["verbose"], "runtime.verbose"),
        outputs=OutputConfig(
            _bool(output_values["save_global_vtp"], "outputs.save_global_vtp"),
            _bool(output_values["save_port_vtp"], "outputs.save_port_vtp"),
            _bool(output_values["save_figures"], "outputs.save_figures"),
        ),
    )
