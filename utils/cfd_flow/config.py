"""Strict YAML configuration for the single production CFD-flow route."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml


METHOD = "PROTEUS_COMPATIBLE_SEEDER_MUSUBI_STEADY_LBM_BASELINE"
SEEDER_COMMIT = "667109df6fafdcb39f4409e3f5d90f04d75cd33c"
MUSUBI_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
PROTEUS_COMMIT = "fc24fcbcf9623006e8f79989c0126c45a76bc23b"
BULK_VISCOSITY_SOURCE = "MUSUBI_D3Q19_REQUIRED_EXPLICIT_PARAMETER"
BULK_VISCOSITY_POLICY = "OFFICIAL_MUSUBI_BASELINE_TWO_THIRDS_KINEMATIC_VISCOSITY"


@dataclass(frozen=True, slots=True)
class PathsConfig:
    source_surface_run: Path
    output_root: Path


@dataclass(frozen=True, slots=True)
class ExecutionConfig:
    mode: str
    parent_failed_run: Path
    frozen_seeder_run: Path
    prior_cumulative_seeder_launches: int
    prior_cumulative_musubi_launches: int
    prior_musubi_runs_reaching_iteration_1: int


@dataclass(frozen=True, slots=True)
class ApesConfig:
    execution_environment: str
    wsl_distribution: str
    seeder_repository: str
    seeder_commit: str
    musubi_repository: str
    musubi_commit: str
    proteus_repository: str
    proteus_reference_commit: str
    seeder_executable: str
    musubi_executable: str
    harvesting_executable: str
    mpi_launcher: str
    mpi_ranks: int | None


@dataclass(frozen=True, slots=True)
class MeshConfig:
    dx_target_um: float
    bounding_margin_cells: int
    uniform_cartesian: bool
    automatic_resolution_sweep: bool

    @property
    def dx_target_m(self) -> float:
        return self.dx_target_um * 1.0e-6


@dataclass(frozen=True, slots=True)
class PhysicsConfig:
    reference_dx_m: float
    reference_dt_s: float
    lattice_cs_squared: float
    bulk_viscosity_source: str
    bulk_viscosity_policy: str


@dataclass(frozen=True, slots=True)
class SolverConfig:
    kind: str
    layout: str
    relaxation: str
    wall_boundary: str
    inlet_boundary: str
    outlet_boundary: str
    maximum_iterations: int
    wallclock_limit_s: int
    convergence_interval_iterations: int
    convergence_nvals: int
    velocity_absolute_threshold_m_s: float
    pressure_absolute_threshold_pa: float
    maximum_lattice_mach: float


@dataclass(frozen=True, slots=True)
class ResourcesConfig:
    maximum_available_ram_fraction: float
    estimated_bytes_per_fluid_cell: int


@dataclass(frozen=True, slots=True)
class QCConfig:
    maximum_mass_conservation_error: float
    maximum_inlet_flow_relative_error: float
    fluid_inside_tolerance_cells: float


@dataclass(frozen=True, slots=True)
class FlowConfig:
    source_path: Path
    method: str
    paths: PathsConfig
    execution: ExecutionConfig
    apes: ApesConfig
    mesh: MeshConfig
    physics: PhysicsConfig
    solver: SolverConfig
    resources: ResourcesConfig
    qc: QCConfig


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a YAML mapping")
    return dict(value)


def _exact(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    mapping = _mapping(value, label)
    missing = sorted(keys - set(mapping))
    unknown = sorted(set(mapping) - keys)
    if missing:
        raise ValueError(f"Missing keys in {label}: {', '.join(missing)}")
    if unknown:
        raise ValueError(f"Unknown keys in {label}: {', '.join(unknown)}")
    return mapping


def _section(root: Mapping[str, Any], key: str, keys: set[str]) -> dict[str, Any]:
    return _exact(root.get(key), keys, key)


def _require(value: Any, expected: Any, label: str) -> None:
    if value != expected:
        raise ValueError(f"{label} must be {expected!r}")


def _number(value: Any, label: str, *, positive: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if positive and result <= 0.0:
        raise ValueError(f"{label} must be positive")
    return result


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    return int(value)


def _boolean(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false")
    return value


def _path(value: Any, label: str, project_root: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty path")
    path = Path(value.strip()).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def _auto_or_path(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be 'auto' or a non-empty executable path")
    return value.strip()


def load_cfd_flow_config(path: Path, *, project_root: Path) -> FlowConfig:
    """Load and enforce the immutable first-baseline configuration."""

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"CFD flow YAML does not exist: {source}")
    root = _exact(
        yaml.safe_load(source.read_text(encoding="utf-8")),
        {"schema_version", "method", "paths", "execution", "apes", "mesh", "physics", "solver", "resources", "qc"},
        "YAML root",
    )
    _require(_integer(root["schema_version"], "schema_version"), 1, "schema_version")
    _require(root["method"], METHOD, "method")
    project = Path(project_root).resolve()

    paths = _section(root, "paths", {"source_surface_run", "output_root"})
    execution = _section(
        root,
        "execution",
        {
            "mode", "parent_failed_run", "frozen_seeder_run",
            "prior_cumulative_seeder_launches", "prior_cumulative_musubi_launches",
            "prior_musubi_runs_reaching_iteration_1",
        },
    )
    _require(execution["mode"], "MUSUBI_ONLY_RECOVERY", "execution.mode")
    _require(_integer(execution["prior_cumulative_seeder_launches"], "execution.prior_cumulative_seeder_launches"), 2, "execution.prior_cumulative_seeder_launches")
    _require(_integer(execution["prior_cumulative_musubi_launches"], "execution.prior_cumulative_musubi_launches"), 1, "execution.prior_cumulative_musubi_launches")
    _require(_integer(execution["prior_musubi_runs_reaching_iteration_1"], "execution.prior_musubi_runs_reaching_iteration_1"), 0, "execution.prior_musubi_runs_reaching_iteration_1")
    apes = _section(
        root,
        "apes",
        {
            "execution_environment", "wsl_distribution", "seeder_repository",
            "seeder_commit", "musubi_repository", "musubi_commit",
            "proteus_repository", "proteus_reference_commit", "seeder_executable",
            "musubi_executable", "harvesting_executable", "mpi_launcher", "mpi_ranks",
        },
    )
    fixed_apes = {
        "execution_environment": "WSL2",
        "seeder_repository": "https://github.com/apes-suite/seeder",
        "seeder_commit": SEEDER_COMMIT,
        "musubi_repository": "https://github.com/apes-suite/musubi",
        "musubi_commit": MUSUBI_COMMIT,
        "proteus_repository": "https://github.com/PROTEUS-SIM/PROTEUS",
        "proteus_reference_commit": PROTEUS_COMMIT,
    }
    for key, expected in fixed_apes.items():
        _require(apes[key], expected, f"apes.{key}")
    if not isinstance(apes["wsl_distribution"], str) or not apes["wsl_distribution"].strip():
        raise ValueError("apes.wsl_distribution must be non-empty")
    ranks_raw = apes["mpi_ranks"]
    if ranks_raw == "auto":
        ranks = None
    else:
        ranks = _integer(ranks_raw, "apes.mpi_ranks")
        if not 1 <= ranks <= 8:
            raise ValueError("apes.mpi_ranks must be in [1, 8] or 'auto'")

    mesh = _section(
        root,
        "mesh",
        {"dx_target_um", "bounding_margin_cells", "uniform_cartesian", "automatic_resolution_sweep"},
    )
    _require(_number(mesh["dx_target_um"], "mesh.dx_target_um"), 0.20, "mesh.dx_target_um")
    _require(_integer(mesh["bounding_margin_cells"], "mesh.bounding_margin_cells"), 4, "mesh.bounding_margin_cells")
    _require(_boolean(mesh["uniform_cartesian"], "mesh.uniform_cartesian"), True, "mesh.uniform_cartesian")
    _require(_boolean(mesh["automatic_resolution_sweep"], "mesh.automatic_resolution_sweep"), False, "mesh.automatic_resolution_sweep")

    physics = _section(
        root,
        "physics",
        {
            "reference_dx_m", "reference_dt_s", "lattice_cs_squared",
            "bulk_viscosity_source", "bulk_viscosity_policy",
        },
    )
    fixed_physics = {"reference_dx_m": 1.28e-7, "reference_dt_s": 1.0e-8, "lattice_cs_squared": 1.0 / 3.0}
    for key, expected in fixed_physics.items():
        actual = _number(physics[key], f"physics.{key}")
        if abs(actual - expected) > abs(expected) * 1.0e-12:
            raise ValueError(f"physics.{key} must be {expected!r}")
    _require(physics["bulk_viscosity_source"], BULK_VISCOSITY_SOURCE, "physics.bulk_viscosity_source")
    _require(physics["bulk_viscosity_policy"], BULK_VISCOSITY_POLICY, "physics.bulk_viscosity_policy")

    solver = _section(
        root,
        "solver",
        {
            "kind", "layout", "relaxation", "wall_boundary", "inlet_boundary",
            "outlet_boundary", "maximum_iterations", "wallclock_limit_s",
            "convergence_interval_iterations", "convergence_nvals",
            "velocity_absolute_threshold_m_s", "pressure_absolute_threshold_pa",
            "maximum_lattice_mach",
        },
    )
    fixed_solver = {
        "kind": "fluid", "layout": "d3q19", "relaxation": "bgk",
        "wall_boundary": "wall_libb", "inlet_boundary": "mfr_eq",
        "outlet_boundary": "pressure_eq", "maximum_iterations": 1_000_000,
        "wallclock_limit_s": 3600, "convergence_interval_iterations": 100,
        "convergence_nvals": 100, "velocity_absolute_threshold_m_s": 1.0e-9,
        "pressure_absolute_threshold_pa": 1.0e-3, "maximum_lattice_mach": 0.05,
    }
    for key, expected in fixed_solver.items():
        value = solver[key]
        if isinstance(expected, int):
            value = _integer(value, f"solver.{key}")
        elif isinstance(expected, float):
            value = _number(value, f"solver.{key}")
        _require(value, expected, f"solver.{key}")

    resources = _section(root, "resources", {"maximum_available_ram_fraction", "estimated_bytes_per_fluid_cell"})
    ram_fraction = _number(resources["maximum_available_ram_fraction"], "resources.maximum_available_ram_fraction")
    _require(ram_fraction, 0.60, "resources.maximum_available_ram_fraction")
    bytes_per_cell = _integer(resources["estimated_bytes_per_fluid_cell"], "resources.estimated_bytes_per_fluid_cell")
    if bytes_per_cell < 1:
        raise ValueError("resources.estimated_bytes_per_fluid_cell must be positive")

    qc = _section(root, "qc", {"maximum_mass_conservation_error", "maximum_inlet_flow_relative_error", "fluid_inside_tolerance_cells"})
    mass_error = _number(qc["maximum_mass_conservation_error"], "qc.maximum_mass_conservation_error")
    inlet_error = _number(qc["maximum_inlet_flow_relative_error"], "qc.maximum_inlet_flow_relative_error")
    inside_tolerance = _number(qc["fluid_inside_tolerance_cells"], "qc.fluid_inside_tolerance_cells")
    _require(mass_error, 0.01, "qc.maximum_mass_conservation_error")
    _require(inlet_error, 0.01, "qc.maximum_inlet_flow_relative_error")
    _require(inside_tolerance, 1.0, "qc.fluid_inside_tolerance_cells")

    return FlowConfig(
        source_path=source,
        method=METHOD,
        paths=PathsConfig(
            source_surface_run=_path(paths["source_surface_run"], "paths.source_surface_run", project),
            output_root=_path(paths["output_root"], "paths.output_root", project),
        ),
        execution=ExecutionConfig(
            mode="MUSUBI_ONLY_RECOVERY",
            parent_failed_run=_path(execution["parent_failed_run"], "execution.parent_failed_run", project),
            frozen_seeder_run=_path(execution["frozen_seeder_run"], "execution.frozen_seeder_run", project),
            prior_cumulative_seeder_launches=2,
            prior_cumulative_musubi_launches=1,
            prior_musubi_runs_reaching_iteration_1=0,
        ),
        apes=ApesConfig(
            execution_environment=str(apes["execution_environment"]),
            wsl_distribution=str(apes["wsl_distribution"]),
            seeder_repository=str(apes["seeder_repository"]),
            seeder_commit=str(apes["seeder_commit"]),
            musubi_repository=str(apes["musubi_repository"]),
            musubi_commit=str(apes["musubi_commit"]),
            proteus_repository=str(apes["proteus_repository"]),
            proteus_reference_commit=str(apes["proteus_reference_commit"]),
            seeder_executable=_auto_or_path(apes["seeder_executable"], "apes.seeder_executable"),
            musubi_executable=_auto_or_path(apes["musubi_executable"], "apes.musubi_executable"),
            harvesting_executable=_auto_or_path(apes["harvesting_executable"], "apes.harvesting_executable"),
            mpi_launcher=_auto_or_path(apes["mpi_launcher"], "apes.mpi_launcher"),
            mpi_ranks=ranks,
        ),
        mesh=MeshConfig(0.20, 4, True, False),
        physics=PhysicsConfig(
            1.28e-7,
            1.0e-8,
            1.0 / 3.0,
            BULK_VISCOSITY_SOURCE,
            BULK_VISCOSITY_POLICY,
        ),
        solver=SolverConfig(
            kind="fluid", layout="d3q19", relaxation="bgk",
            wall_boundary="wall_libb", inlet_boundary="mfr_eq",
            outlet_boundary="pressure_eq", maximum_iterations=1_000_000,
            wallclock_limit_s=3600, convergence_interval_iterations=100,
            convergence_nvals=100, velocity_absolute_threshold_m_s=1.0e-9,
            pressure_absolute_threshold_pa=1.0e-3, maximum_lattice_mach=0.05,
        ),
        resources=ResourcesConfig(ram_fraction, bytes_per_cell),
        qc=QCConfig(mass_error, inlet_error, inside_tolerance),
    )
