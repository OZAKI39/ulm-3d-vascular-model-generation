"""Schema-v2 configuration for the validated production CFD route."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml

from .validated_contract import (
    BASE_DX_M,
    BULK_VISCOSITY_M2_S,
    FULL_TIMESTEP_IDENTITY_GATE,
    KINEMATIC_VISCOSITY_M2_S,
    LATTICE_CS2,
    MAXIMUM_LATTICE_SPEED,
    OUTLET_GAUGE_PRESSURES_PA,
    RHO0_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_VOLUME_FLOW_M3_S,
)


SCHEMA_VERSION = 2
METHOD = "VALIDATED_TAU1_ADAPTIVE_FLUX_CONTINUOUS_Q_STEADY_LBM"
EXECUTION_MODES = {"FRESH_STEADY", "VALIDATED_BASE_PROMOTION_REPLAY"}
SEEDER_COMMIT = "667109df6fafdcb39f4409e3f5d90f04d75cd33c"
MUSUBI_UPSTREAM_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
MUSUBI_SCHEME_COMMIT = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
TREELM_COMMIT = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
PROTEUS_COMMIT = "fc24fcbcf9623006e8f79989c0126c45a76bc23b"
PLANE_CONTRACT_SHA256 = "ffaa49bdb6e43fb7208ff29df07a90d4e92ef9bfa4b96ca4f997d4f453a7f005"


@dataclass(frozen=True, slots=True)
class PathsConfig:
    source_surface_run: Path
    output_root: Path
    frozen_base_mesh: Path
    accepted_base_restart_header: Path
    accepted_base_restart_binary: Path
    accepted_base_qc: Path
    accepted_base_checkpoint_history: Path
    accepted_base_full_v2: Path
    physical_plane_contract: Path
    coarse_base_grid_evidence: Path


@dataclass(frozen=True, slots=True)
class ExecutionConfig:
    mode: str
    solver_smoke_iterations: int
    fresh_maximum_iterations: int
    run_solver_smoke: bool


@dataclass(frozen=True, slots=True)
class ApesConfig:
    execution_environment: str
    wsl_distribution: str
    seeder_repository: str
    seeder_commit: str
    seeder_executable: str
    seeder_expected_sha256: str
    seeder_patch_sha256: str
    musubi_repository: str
    musubi_upstream_commit: str
    musubi_scheme_commit: str
    treelm_commit: str
    musubi_executable: str
    musubi_expected_sha256: str
    musubi_patch_sha256: str
    musubi_patched_source_sha256: str
    harvesting_executable: str
    harvesting_expected_sha256: str
    mpi_launcher: str
    mpi_ranks: int
    proteus_repository: str
    proteus_reference_commit: str


@dataclass(frozen=True, slots=True)
class MeshConfig:
    dx_m: float
    bounding_margin_cells: int
    uniform_cartesian: bool
    expected_cells: int
    elemlist_sha256: str
    bnd_sha256: str
    qval_sha256: str

    @property
    def dx_target_m(self) -> float:
        return self.dx_m

    @property
    def dx_target_um(self) -> float:
        return self.dx_m * 1.0e6


@dataclass(frozen=True, slots=True)
class PhysicsConfig:
    density_kg_m3: float
    kinematic_viscosity_m2_s: float
    bulk_viscosity_m2_s: float
    lattice_cs_squared: float
    lattice_kinematic_viscosity: float


@dataclass(frozen=True, slots=True)
class BoundaryConfig:
    physical_inlet_requirement: str
    inlet_boundary: str
    target_mass_flow_kg_s: float
    target_volume_flow_m3_s: float
    wall_boundary: str
    outlet_boundary: str
    outlet_gauge_pressures_pa: tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class SolverConfig:
    kind: str
    layout: str
    relaxation: str
    maximum_lattice_speed: float
    controller_target_error: float
    controller_controlled_flux_error: float
    full_timestep_identity_gate: float


@dataclass(frozen=True, slots=True)
class SteadyAcceptanceConfig:
    short_window_s: float
    long_window_s: float
    mass_residual: float
    physical_volume_closure: float
    velocity_residual: float
    pressure_residual: float
    inlet_residual: float
    flow_fraction_drift: float
    rho_min: float
    rho_max: float


@dataclass(frozen=True, slots=True)
class ResourcesConfig:
    maximum_available_ram_fraction: float
    estimated_bytes_per_fluid_cell: int
    wallclock_limit_s: int


@dataclass(frozen=True, slots=True)
class VisualizationConfig:
    width_px: int
    height_px: int
    dpi: int
    pressure_field: str
    velocity_units: str


@dataclass(frozen=True, slots=True)
class FlowConfig:
    source_path: Path
    schema_version: int
    method: str
    paths: PathsConfig
    execution: ExecutionConfig
    apes: ApesConfig
    mesh: MeshConfig
    physics: PhysicsConfig
    boundary: BoundaryConfig
    solver: SolverConfig
    steady_acceptance: SteadyAcceptanceConfig
    resources: ResourcesConfig
    visualization: VisualizationConfig


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a YAML mapping")
    return dict(value)


def _exact(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    result = _mapping(value, label)
    missing = sorted(keys - set(result))
    unknown = sorted(set(result) - keys)
    if missing or unknown:
        raise ValueError(f"{label}: missing={missing}, unknown={unknown}")
    return result


def _path(value: Any, label: str, project_root: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty path")
    path = Path(value.strip()).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def _string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")
    return value.strip()


def _bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false")
    return value


def _int(value: Any, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _number(value: Any, label: str, *, positive: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if positive and result <= 0.0:
        raise ValueError(f"{label} must be positive")
    return result


def _fixed(actual: Any, expected: Any, label: str) -> None:
    if isinstance(expected, float):
        value = _number(actual, label, positive=expected > 0.0)
        if value != expected:
            raise ValueError(f"{label} must be {expected!r}")
    elif actual != expected:
        raise ValueError(f"{label} must be {expected!r}")


def _sha(value: Any, label: str) -> str:
    result = _string(value, label).lower()
    if len(result) != 64 or any(char not in "0123456789abcdef" for char in result):
        raise ValueError(f"{label} must be a 64-character SHA-256")
    return result


def load_cfd_flow_config(path: Path, *, project_root: Path) -> FlowConfig:
    """Load the strict production schema; schema-v1 cannot silently migrate."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"CFD flow YAML does not exist: {source}")
    root = _exact(
        yaml.safe_load(source.read_text(encoding="utf-8")),
        {
            "schema_version", "method", "paths", "execution", "apes", "mesh",
            "physics", "boundary_conditions", "solver", "steady_acceptance",
            "resources", "visualization",
        },
        "YAML root",
    )
    _fixed(root["schema_version"], SCHEMA_VERSION, "schema_version")
    _fixed(root["method"], METHOD, "method")
    project = Path(project_root).resolve()

    path_keys = {
        "source_surface_run", "output_root", "frozen_base_mesh",
        "accepted_base_restart_header", "accepted_base_restart_binary",
        "accepted_base_qc", "accepted_base_checkpoint_history",
        "accepted_base_full_v2", "physical_plane_contract",
        "coarse_base_grid_evidence",
    }
    paths = _exact(root["paths"], path_keys, "paths")
    parsed_paths = {key: _path(value, f"paths.{key}", project) for key, value in paths.items()}

    execution = _exact(
        root["execution"],
        {"mode", "solver_smoke_iterations", "fresh_maximum_iterations", "run_solver_smoke"},
        "execution",
    )
    mode = _string(execution["mode"], "execution.mode")
    if mode not in EXECUTION_MODES:
        raise ValueError(f"execution.mode must be one of {sorted(EXECUTION_MODES)}")
    smoke_iterations = _int(execution["solver_smoke_iterations"], "execution.solver_smoke_iterations")
    if smoke_iterations > 5000:
        raise ValueError("execution.solver_smoke_iterations cannot exceed 5000")

    apes_keys = {
        "execution_environment", "wsl_distribution", "seeder_repository", "seeder_commit",
        "seeder_executable", "seeder_expected_sha256", "seeder_patch_sha256",
        "musubi_repository", "musubi_upstream_commit", "musubi_scheme_commit", "treelm_commit",
        "musubi_executable", "musubi_expected_sha256", "musubi_patch_sha256",
        "musubi_patched_source_sha256", "harvesting_executable",
        "harvesting_expected_sha256", "mpi_launcher", "mpi_ranks",
        "proteus_repository", "proteus_reference_commit",
    }
    apes = _exact(root["apes"], apes_keys, "apes")
    for key, expected in {
        "execution_environment": "WSL2", "seeder_commit": SEEDER_COMMIT,
        "musubi_upstream_commit": MUSUBI_UPSTREAM_COMMIT,
        "musubi_scheme_commit": MUSUBI_SCHEME_COMMIT, "treelm_commit": TREELM_COMMIT,
        "proteus_reference_commit": PROTEUS_COMMIT,
    }.items():
        _fixed(apes[key], expected, f"apes.{key}")

    mesh = _exact(
        root["mesh"],
        {"dx_m", "bounding_margin_cells", "uniform_cartesian", "expected_cells", "elemlist_sha256", "bnd_sha256", "qval_sha256"},
        "mesh",
    )
    _fixed(mesh["dx_m"], BASE_DX_M, "mesh.dx_m")
    _fixed(mesh["uniform_cartesian"], True, "mesh.uniform_cartesian")
    _fixed(mesh["expected_cells"], 182_320, "mesh.expected_cells")

    physics = _exact(
        root["physics"],
        {"density_kg_m3", "kinematic_viscosity_m2_s", "bulk_viscosity_m2_s", "lattice_cs_squared", "lattice_kinematic_viscosity"},
        "physics",
    )
    for key, expected in {
        "density_kg_m3": RHO0_KG_M3,
        "kinematic_viscosity_m2_s": KINEMATIC_VISCOSITY_M2_S,
        "bulk_viscosity_m2_s": BULK_VISCOSITY_M2_S,
        "lattice_cs_squared": LATTICE_CS2,
        "lattice_kinematic_viscosity": 1.0 / 6.0,
    }.items():
        _fixed(physics[key], expected, f"physics.{key}")

    boundary = _exact(
        root["boundary_conditions"],
        {"physical_inlet_requirement", "inlet_boundary", "target_mass_flow_kg_s", "target_volume_flow_m3_s", "wall_boundary", "outlet_boundary", "outlet_gauge_pressures_pa"},
        "boundary_conditions",
    )
    for key, expected in {
        "physical_inlet_requirement": "GLOBAL_TARGET_VOLUMETRIC_FLOW",
        "inlet_boundary": "adaptive_flux_pressure", "wall_boundary": "wall_libb",
        "outlet_boundary": "pressure_eq", "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
        "target_volume_flow_m3_s": TARGET_VOLUME_FLOW_M3_S,
    }.items():
        _fixed(boundary[key], expected, f"boundary_conditions.{key}")
    gauges = _mapping(boundary["outlet_gauge_pressures_pa"], "boundary_conditions.outlet_gauge_pressures_pa")
    if set(gauges) != set(OUTLET_GAUGE_PRESSURES_PA):
        raise ValueError("outlet gauge labels must be outlet_01/outlet_02/outlet_03")
    for key, expected in OUTLET_GAUGE_PRESSURES_PA.items():
        _fixed(gauges[key], expected, f"boundary_conditions.outlet_gauge_pressures_pa.{key}")

    solver = _exact(
        root["solver"],
        {"kind", "layout", "relaxation", "maximum_lattice_speed", "controller_target_error", "controller_controlled_flux_error", "full_timestep_identity_gate"},
        "solver",
    )
    for key, expected in {
        "kind": "fluid", "layout": "d3q19", "relaxation": "bgk",
        "maximum_lattice_speed": MAXIMUM_LATTICE_SPEED,
        "controller_target_error": 1.0e-8, "controller_controlled_flux_error": 1.0e-8,
        "full_timestep_identity_gate": FULL_TIMESTEP_IDENTITY_GATE,
    }.items():
        _fixed(solver[key], expected, f"solver.{key}")

    steady = _exact(
        root["steady_acceptance"],
        {"short_window_s", "long_window_s", "mass_residual", "physical_volume_closure", "velocity_residual", "pressure_residual", "inlet_residual", "flow_fraction_drift", "rho_min", "rho_max"},
        "steady_acceptance",
    )
    expected_steady = {
        "short_window_s": 0.0002441406727828746, "long_window_s": 0.0004882813455657492,
        "mass_residual": 0.01, "physical_volume_closure": 0.01,
        "velocity_residual": 0.01, "pressure_residual": 0.005,
        "inlet_residual": 0.01, "flow_fraction_drift": 0.01,
        "rho_min": 0.9, "rho_max": 1.1,
    }
    for key, expected in expected_steady.items():
        _fixed(steady[key], expected, f"steady_acceptance.{key}")

    resources = _exact(root["resources"], {"maximum_available_ram_fraction", "estimated_bytes_per_fluid_cell", "wallclock_limit_s"}, "resources")
    visualization = _exact(root["visualization"], {"width_px", "height_px", "dpi", "pressure_field", "velocity_units"}, "visualization")
    _fixed(visualization["pressure_field"], "pressure_gauge_pa", "visualization.pressure_field")
    _fixed(visualization["velocity_units"], "mm/s", "visualization.velocity_units")

    return FlowConfig(
        source_path=source,
        schema_version=SCHEMA_VERSION,
        method=METHOD,
        paths=PathsConfig(**parsed_paths),
        execution=ExecutionConfig(
            mode=mode,
            solver_smoke_iterations=smoke_iterations,
            fresh_maximum_iterations=_int(execution["fresh_maximum_iterations"], "execution.fresh_maximum_iterations", minimum=1),
            run_solver_smoke=_bool(execution["run_solver_smoke"], "execution.run_solver_smoke"),
        ),
        apes=ApesConfig(
            execution_environment="WSL2",
            wsl_distribution=_string(apes["wsl_distribution"], "apes.wsl_distribution"),
            seeder_repository=_string(apes["seeder_repository"], "apes.seeder_repository"),
            seeder_commit=SEEDER_COMMIT,
            seeder_executable=_string(apes["seeder_executable"], "apes.seeder_executable"),
            seeder_expected_sha256=_sha(apes["seeder_expected_sha256"], "apes.seeder_expected_sha256"),
            seeder_patch_sha256=_sha(apes["seeder_patch_sha256"], "apes.seeder_patch_sha256"),
            musubi_repository=_string(apes["musubi_repository"], "apes.musubi_repository"),
            musubi_upstream_commit=MUSUBI_UPSTREAM_COMMIT,
            musubi_scheme_commit=MUSUBI_SCHEME_COMMIT,
            treelm_commit=TREELM_COMMIT,
            musubi_executable=_string(apes["musubi_executable"], "apes.musubi_executable"),
            musubi_expected_sha256=_sha(apes["musubi_expected_sha256"], "apes.musubi_expected_sha256"),
            musubi_patch_sha256=_sha(apes["musubi_patch_sha256"], "apes.musubi_patch_sha256"),
            musubi_patched_source_sha256=_sha(apes["musubi_patched_source_sha256"], "apes.musubi_patched_source_sha256"),
            harvesting_executable=_string(apes["harvesting_executable"], "apes.harvesting_executable"),
            harvesting_expected_sha256=_sha(apes["harvesting_expected_sha256"], "apes.harvesting_expected_sha256"),
            mpi_launcher=_string(apes["mpi_launcher"], "apes.mpi_launcher"),
            mpi_ranks=_int(apes["mpi_ranks"], "apes.mpi_ranks", minimum=1),
            proteus_repository=_string(apes["proteus_repository"], "apes.proteus_repository"),
            proteus_reference_commit=PROTEUS_COMMIT,
        ),
        mesh=MeshConfig(
            dx_m=BASE_DX_M,
            bounding_margin_cells=_int(mesh["bounding_margin_cells"], "mesh.bounding_margin_cells", minimum=1),
            uniform_cartesian=True,
            expected_cells=182_320,
            elemlist_sha256=_sha(mesh["elemlist_sha256"], "mesh.elemlist_sha256"),
            bnd_sha256=_sha(mesh["bnd_sha256"], "mesh.bnd_sha256"),
            qval_sha256=_sha(mesh["qval_sha256"], "mesh.qval_sha256"),
        ),
        physics=PhysicsConfig(RHO0_KG_M3, KINEMATIC_VISCOSITY_M2_S, BULK_VISCOSITY_M2_S, LATTICE_CS2, 1.0 / 6.0),
        boundary=BoundaryConfig(
            "GLOBAL_TARGET_VOLUMETRIC_FLOW", "adaptive_flux_pressure",
            TARGET_MASS_FLOW_KG_S, TARGET_VOLUME_FLOW_M3_S, "wall_libb", "pressure_eq",
            tuple(OUTLET_GAUGE_PRESSURES_PA[key] for key in sorted(OUTLET_GAUGE_PRESSURES_PA)),
        ),
        solver=SolverConfig("fluid", "d3q19", "bgk", MAXIMUM_LATTICE_SPEED, 1.0e-8, 1.0e-8, FULL_TIMESTEP_IDENTITY_GATE),
        steady_acceptance=SteadyAcceptanceConfig(**{key: float(value) for key, value in steady.items()}),
        resources=ResourcesConfig(
            _number(resources["maximum_available_ram_fraction"], "resources.maximum_available_ram_fraction"),
            _int(resources["estimated_bytes_per_fluid_cell"], "resources.estimated_bytes_per_fluid_cell", minimum=1),
            _int(resources["wallclock_limit_s"], "resources.wallclock_limit_s", minimum=1),
        ),
        visualization=VisualizationConfig(
            _int(visualization["width_px"], "visualization.width_px", minimum=1200),
            _int(visualization["height_px"], "visualization.height_px", minimum=800),
            _int(visualization["dpi"], "visualization.dpi", minimum=72),
            "pressure_gauge_pa", "mm/s",
        ),
    )
