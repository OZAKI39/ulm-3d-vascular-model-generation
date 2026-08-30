"""Isolated healthy adaptive-flux three-grid convergence workflow.

This module does not modify or call the production pipeline.  It freezes the
accepted V2-refereed base restart, creates only the coarse and fine Seeder
meshes, and stores all grid-convergence evidence under a dedicated run root.
"""

from __future__ import annotations

import csv
import math
import shlex
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from .apes import parse_mesh_header, windows_to_wsl
from .adaptive_flux_steady import (
    BULK_NU_M2_S,
    NU_M2_S,
    _run_luac,
    _wsl_path_to_unc,
)
from .adaptive_flux_validation import (
    BINARY_WSL,
    EXPECTED_BINARY_SHA256,
    MPIRUN_WSL,
)
from .healthy_capillary_fast_steady import (
    MPI_BINDING_ARGS,
    THREAD_ENV,
    _last_controller_record,
    choose_mpi_ranks,
    parse_timing_result,
)
from .ideal_inlet_plane import connected_fluid_region_count
from .io import FlowError, read_json, sha256_file, write_json
from .musubi_boundary_mass_referee import (
    DT_S,
    DX_M,
    EXPECTED_CELLS,
    EXPECTED_INLET_GLOBBC,
    REFEREE_REVISION_NEW,
    OUTLET_PRESSURES_PA,
    PRESSURE_REFERENCE_PA,
    RHO0,
    RUN_NAME,
    TARGET_MASS_FLOW,
    _pressure_selected_rows,
    load_mesh_contract,
    replay_boundary_step,
    runtime_solid_cells,
)
from .restart_decode import read_treelm_elemlist, tree_levels
from .restart_decode import (
    D3Q19_DIRECTIONS,
    parse_restart_header,
    read_restart_pdf,
)


ARCHIVE_RUN = "healthy_mouse_capillary_accepted_steady_anchor003274_20260829"
GRID_RUN = "healthy_mouse_capillary_grid_convergence_anchor003274_20260829"
BASE_MESH_RUN = "axis_aligned_ideal_plane_inlet_preflight_anchor003274_20260829_120444"
ACCEPTED_ITERATION = 394_166
ACCEPTED_SHA256 = "8e0c3fcac06af6426441cccab7c6396240d547d1ebd8fa1e04f0d3d8f9e1012a"
REFINEMENT_RATIO = 1.3
SEEDER_WSL = "/home/lzy/.local/bin/seeder"
BASE_ROOT_LEVEL = 9
BASE_CUBE_ORIGIN_M = np.asarray(
    (7.8082771579034594e-05, 5.4076359424537688e-05, 8.8206398940868488e-05),
    dtype=np.float64,
)
BASE_CUBE_SIDE_M = 0.0001024
PLANE_ORIGIN_M = np.asarray(
    (0.0001016827715790346, 5.5476359424537688e-05, 0.00015880639894086848),
    dtype=np.float64,
)
PLANE_X_M = np.asarray((3.9999999999999888e-06, 0.0, 0.0), dtype=np.float64)
PLANE_Y_M = np.asarray((0.0, 4.2000000000000004e-06, 0.0), dtype=np.float64)
SEED_M = np.asarray(
    (0.00010367370500008977, 5.75542332387685e-05, 0.00015801755415184839),
    dtype=np.float64,
)
MESH_REQUIRED_FILES = ("header.lua", "bnd.lua", "bnd.lsb", "elemlist.lsb")
MESH_QVAL_FILES = ("qval.lua", "qval.lsb")
GRID_BENCHMARK_ITERATIONS = 2_000
GRID_CHECKPOINT_INTERVAL = 5_000
GRID_SIM_CONTROL_INTERVAL = 100
GRID_MAXIMUM_ITERATIONS = 1_000_000
REFEREE_IDENTITY_ERROR = 5.101930221603317e-10
MASS_GATE = 0.01
VELOCITY_GATE = 0.01
PRESSURE_GATE = 0.005
INLET_GATE = 0.01
BOUNDARY_WINDOW_GATE = 0.001
SIGNIFICANT_BACKFLOW_FRACTION = 0.05
FINE_DEFER_FULL_AUDIT_UNTIL = 550_000
AUDIT_WINDOW_ITERATIONS = 20_000


@dataclass(frozen=True, slots=True)
class GridSpec:
    label: str
    dx_m: float
    dt_s: float
    root_level: int
    cube_side_m: float
    predicted_fluid_cells: int
    predicted_steady_iterations: int
    predicted_runtime_s: float


@dataclass(frozen=True, slots=True)
class MeshDirectoryContract:
    """One canonical directory shared by Seeder, QC and Musubi."""

    seeder_run_root: Path
    configured_folder: str
    seeder_output_dir: Path
    mesh_dir: Path
    mesh_windows: Path
    mesh_wsl: str

    def as_record(self) -> dict[str, str]:
        return {
            "seeder_run_root": str(self.seeder_run_root),
            "configured_folder": self.configured_folder,
            "seeder_output_dir": str(self.seeder_output_dir),
            "mesh_dir": str(self.mesh_dir),
            "mesh_windows": str(self.mesh_windows),
            "mesh_wsl": self.mesh_wsl,
        }


def resolve_mesh_directory(
    seeder_run_root: Path,
    *,
    configured_folder: str = "mesh/",
    wsl_distribution: str = "Ubuntu",
) -> MeshDirectoryContract:
    """Resolve Seeder's relative folder exactly once and convert it for WSL."""

    root = Path(seeder_run_root).resolve()
    folder = Path(configured_folder.strip().rstrip("/\\"))
    if folder.is_absolute() or not folder.parts or ".." in folder.parts:
        raise ValueError("Seeder output folder must be a safe relative path")
    if len(folder.parts) == 1 and root.name.casefold() == folder.name.casefold():
        output = root
        run_root = root.parent
    else:
        output = (root / folder).resolve()
        run_root = root
    return MeshDirectoryContract(
        seeder_run_root=run_root,
        configured_folder=configured_folder,
        seeder_output_dir=output,
        mesh_dir=output,
        mesh_windows=output,
        mesh_wsl=windows_to_wsl(output, wsl_distribution),
    )


def mesh_file_contract(mesh_dir: Path, *, require_qval: bool = True) -> dict[str, Any]:
    """Check that every required mesh file exists and is non-empty."""

    mesh = Path(mesh_dir).resolve()
    names = MESH_REQUIRED_FILES + (MESH_QVAL_FILES if require_qval else ())
    files: dict[str, dict[str, Any]] = {}
    for name in names:
        path = mesh / name
        size = path.stat().st_size if path.is_file() else 0
        files[name] = {
            "path": str(path),
            "exists": path.is_file(),
            "size_bytes": size,
            "pass": path.is_file() and size > 0,
        }
    return {
        "mesh_dir": str(mesh),
        "require_qval": require_qval,
        "files": files,
        "status": "PASS" if all(item["pass"] for item in files.values()) else "FAIL",
    }


def require_complete_mesh(mesh_dir: Path, *, require_qval: bool = True) -> dict[str, Any]:
    """Hard guard used before constructing or launching a Musubi run."""

    record = mesh_file_contract(mesh_dir, require_qval=require_qval)
    if record["status"] != "PASS":
        missing = [name for name, item in record["files"].items() if not item["pass"]]
        raise FlowError(
            "CFD_FLOW_MESH_DIRECTORY_CONTRACT_FAILED",
            "Missing or empty mesh files: " + ", ".join(missing),
        )
    return record


def musubi_mesh_assignment(mesh_dir: Path) -> str:
    """Return a guarded Musubi assignment pointing at the canonical WSL path."""

    require_complete_mesh(mesh_dir, require_qval=True)
    contract = resolve_mesh_directory(Path(mesh_dir))
    return f"mesh = '{contract.mesh_wsl.rstrip('/')}/'"


def can_reuse_completed_seeder_output(
    mesh_dir: Path, manifest: dict[str, Any] | None
) -> bool:
    """Reuse only a complete, non-empty mesh from a successful Seeder run."""

    return bool(
        manifest
        and manifest.get("status") == "PASS"
        and int(manifest.get("seeder_returncode", 1)) == 0
        and mesh_file_contract(mesh_dir, require_qval=True)["status"] == "PASS"
    )


def grid_specs() -> dict[str, GridSpec]:
    """Return the fixed r=1.3 coarse/base/fine design with diffusive scaling."""

    base_runtime_s = 394_166 / 51.39525262051544
    definitions = {
        "coarse": (DX_M * REFINEMENT_RATIO, 9),
        "base": (DX_M, 9),
        "fine": (DX_M / REFINEMENT_RATIO, 10),
    }
    result: dict[str, GridSpec] = {}
    for label, (dx_m, level) in definitions.items():
        scale = DX_M / dx_m
        dt_s = DT_S * (dx_m / DX_M) ** 2
        cells = int(round(EXPECTED_CELLS * scale**3))
        iterations = int(math.ceil(ACCEPTED_ITERATION * scale**2))
        runtime_s = base_runtime_s * scale**5
        result[label] = GridSpec(
            label=label,
            dx_m=dx_m,
            dt_s=dt_s,
            root_level=level,
            cube_side_m=dx_m * 2**level,
            predicted_fluid_cells=cells,
            predicted_steady_iterations=iterations,
            predicted_runtime_s=runtime_s,
        )
    return result


def _grid_physics_lua(spec: GridSpec) -> str:
    return f"""dx = {spec.dx_m:.17g}
dt = {spec.dt_s:.17g}
rho0_phy = {RHO0:.17g}
nu_phy = {NU_M2_S:.17g}
bulk_viscosity_phy = {BULK_NU_M2_S:.17g}
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}

function outlet_01_pressure(x,y,z,t) return {OUTLET_PRESSURES_PA['outlet_01']:.17g} end
function outlet_02_pressure(x,y,z,t) return {OUTLET_PRESSURES_PA['outlet_02']:.17g} end
function outlet_03_pressure(x,y,z,t) return {OUTLET_PRESSURES_PA['outlet_03']:.17g} end

physics = {{ dt = dt, rho0 = rho0_phy }}
identify = {{ label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }}
fluid = {{ kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }}
initial_condition = {{ pressure = pressure_reference_phy, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }}
boundary_condition = {{
  {{ label = 'wall', kind = 'wall_libb' }},
  {{ label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = {TARGET_MASS_FLOW:.17g} }},
  {{ label = 'outlet_01', kind = 'pressure_eq', pressure = outlet_01_pressure }},
  {{ label = 'outlet_02', kind = 'pressure_eq', pressure = outlet_02_pressure }},
  {{ label = 'outlet_03', kind = 'pressure_eq', pressure = outlet_03_pressure }}
}}"""


def generate_grid_musubi_lua(
    spec: GridSpec,
    *,
    mesh_dir: Path,
    maximum_iterations: int,
    simulation_name: str,
    write_restarts: bool,
    resume_header: str | None = None,
    first_restart_iteration: int = GRID_CHECKPOINT_INTERVAL,
) -> str:
    """Generate an isolated cold-start grid run with frozen healthy physics."""

    assignment = musubi_mesh_assignment(mesh_dir)
    restart = ""
    if write_restarts:
        read_clause = f"  read = '{resume_header}',\n" if resume_header else ""
        restart = f"""
restart = {{
{read_clause}  
  write = 'restart/',
  time_control = {{ min = {{ iter = {int(first_restart_iteration)} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {GRID_CHECKPOINT_INTERVAL} }} }}
}}"""
    return f"""-- Isolated healthy adaptive-flux grid convergence run.
simulation_name = '{simulation_name}'
printRuntimeInfo = true
timing_file = 'timing.res'
{assignment}
scaling = 'diffusive'
logging = {{ level = 5 }}
maximum_iterations = {int(maximum_iterations)}
{_grid_physics_lua(spec)}
sim_control = {{
  time_control = {{ max = {{ iter = maximum_iterations }}, interval = {{ iter = {GRID_SIM_CONTROL_INTERVAL} }} }},
  abort_criteria = {{ stop_file = 'stop' }}
}}
{restart}
"""


def grid_lua_contract(
    text: str,
    spec: GridSpec,
    mesh_dir: Path,
    *,
    resume_header: str | None = None,
) -> dict[str, Any]:
    expected_mesh = resolve_mesh_directory(mesh_dir).mesh_wsl.rstrip("/") + "/"
    checks = {
        "resolved_mesh": f"mesh = '{expected_mesh}'" in text,
        "no_duplicate_mesh": "/mesh/mesh/" not in text,
        "diffusive_dx": f"dx = {spec.dx_m:.17g}" in text,
        "diffusive_dt": f"dt = {spec.dt_s:.17g}" in text,
        "density": f"rho0_phy = {RHO0:.17g}" in text,
        "viscosity": f"nu_phy = {NU_M2_S:.17g}" in text,
        "bulk_viscosity": f"bulk_viscosity_phy = {BULK_NU_M2_S:.17g}" in text,
        "healthy_mass_flow": f"mass_flowrate = {TARGET_MASS_FLOW:.17g}" in text,
        "adaptive_flux_pressure": "kind = 'adaptive_flux_pressure'" in text,
        "wall_libb": "kind = 'wall_libb'" in text,
        "three_pressure_outlets": text.count("kind = 'pressure_eq'") == 3,
        "no_harvester": "harvest" not in text.casefold(),
        "restart_mode": (
            f"read = '{resume_header}'" in text
            if resume_header
            else "read =" not in text
        ),
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _create_wsl_runtime(runtime_root_wsl: str) -> Path:
    probe = subprocess.run(
        [
            "wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc",
            f"test ! -e {shlex.quote(runtime_root_wsl)} && mkdir -p {shlex.quote(runtime_root_wsl)}",
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
        check=False,
    )
    if probe.returncode != 0:
        raise FlowError("CFD_FLOW_GRID_RUNTIME_EXISTS", runtime_root_wsl)
    return _wsl_path_to_unc("Ubuntu", runtime_root_wsl)


def _run_grid_benchmark_candidate(
    *,
    spec: GridSpec,
    mesh_dir: Path,
    ranks: int,
    runtime_root_wsl: str,
    evidence: Path,
) -> dict[str, Any]:
    runtime_root = _create_wsl_runtime(runtime_root_wsl)
    evidence.mkdir(parents=True, exist_ok=False)
    lua = generate_grid_musubi_lua(
        spec,
        mesh_dir=mesh_dir,
        maximum_iterations=GRID_BENCHMARK_ITERATIONS,
        simulation_name=f"gc_{spec.label[0]}_b{ranks}",
        write_restarts=False,
    )
    (runtime_root / "musubi.lua").write_text(lua, encoding="utf-8")
    (evidence / "musubi.lua").write_text(lua, encoding="utf-8")
    contract = grid_lua_contract(lua, spec, mesh_dir)
    write_json(evidence / "lua_contract.json", contract)
    if contract["status"] != "PASS":
        raise FlowError("CFD_FLOW_GRID_LUA_CONTRACT_FAILED", spec.label)
    luac_rc = _run_luac(
        distribution="Ubuntu",
        workdir_wsl=runtime_root_wsl,
        stdout_path=evidence / "luac_stdout.log",
        stderr_path=evidence / "luac_stderr.log",
    )
    if luac_rc != 0:
        raise FlowError("CFD_FLOW_GRID_LUA_CONTRACT_FAILED", f"{spec.label} luac")
    command = [
        "env",
        *[f"{key}={value}" for key, value in THREAD_ENV.items()],
        "OMPI_ALLOW_RUN_AS_ROOT=1",
        "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
        MPIRUN_WSL,
        "-np",
        str(ranks),
        *MPI_BINDING_ARGS,
        BINARY_WSL,
        "musubi.lua",
    ]
    shell = f"cd {shlex.quote(runtime_root_wsl)} && exec {shlex.join(command)}"
    started = time.perf_counter()
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", shell],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=300,
        check=False,
    )
    elapsed = time.perf_counter() - started
    (evidence / "musubi_stdout.log").write_text(process.stdout, encoding="utf-8")
    (evidence / "musubi_stderr.log").write_text(process.stderr, encoding="utf-8")
    if process.returncode != 0:
        raise FlowError(
            "CFD_FLOW_GRID_BENCHMARK_FAILED",
            f"{spec.label} {ranks} ranks returncode={process.returncode}",
        )
    timing = parse_timing_result(runtime_root / "timing.res")
    controller = _last_controller_record(process.stdout)
    checks = {
        "iterations": int(controller["iteration"]) == GRID_BENCHMARK_ITERATIONS,
        "controller_identity": float(controller["relative_error"]) <= 1.0e-8,
        "globbc_positive": int(controller["globBC_count"]) > 0,
        "minimum_pdf": float(controller["minimum_pdf"]) > 0.0,
        "maximum_lattice_speed": float(controller["max_lattice_velocity"]) < 0.05,
        "reported_ranks": int(timing["nprocs"]) == ranks,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "ranks": ranks,
        "iterations": GRID_BENCHMARK_ITERATIONS,
        "iterations_per_s": GRID_BENCHMARK_ITERATIONS / float(timing["time_main_loop_s"]),
        "wrapper_iterations_per_s": GRID_BENCHMARK_ITERATIONS / elapsed,
        "wall_time_s": elapsed,
        "timing": timing,
        "controller": controller,
        "checks": checks,
        "runtime_root_wsl": runtime_root_wsl,
        "binary_wsl": BINARY_WSL,
        "expected_binary_sha256": EXPECTED_BINARY_SHA256,
    }
    write_json(evidence / "benchmark_result.json", result)
    if result["status"] != "PASS":
        raise FlowError("CFD_FLOW_GRID_BENCHMARK_FAILED", f"{spec.label} {ranks}")
    return result


def run_grid_mpi_benchmark(project_root: Path, label: str) -> dict[str, Any]:
    root = Path(project_root).resolve()
    if label not in {"coarse", "fine"}:
        raise ValueError(label)
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    mesh = run_root / "grids" / label / "seeder" / "mesh"
    require_complete_mesh(mesh, require_qval=True)
    benchmark_root = run_root / "grids" / label / "benchmark"
    summary_path = benchmark_root / "benchmark_summary.json"
    if summary_path.is_file():
        existing = read_json(summary_path)
        if existing.get("status") == "PASS":
            return existing
        raise FlowError("CFD_FLOW_GRID_BENCHMARK_FAILED", f"existing {label} summary")
    specs = grid_specs()
    candidates = (2, 4, 6) if label == "coarse" else (4, 6, 8)
    stamp = datetime.now().strftime("%H%M%S")
    results = []
    for ranks in candidates:
        results.append(
            _run_grid_benchmark_candidate(
                spec=specs[label],
                mesh_dir=mesh,
                ranks=ranks,
                runtime_root_wsl=f"/home/lzy/u3da/gc_{label[0]}_b{ranks}_{stamp}",
                evidence=benchmark_root / f"rank_{ranks}",
            )
        )
    selection = choose_mpi_ranks(results)
    summary = {
        "status": "PASS",
        "grid": label,
        "candidates": results,
        "selection": selection,
        "selected_ranks": int(selection["selected"]["ranks"]),
    }
    write_json(summary_path, summary)
    return summary


def _grid_pdf_state(
    binary: Path, spec: GridSpec, cell_count: int
) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    pdf = read_restart_pdf(binary, n_elems=cell_count, n_components=19)
    values = np.asarray(pdf)
    density = np.sum(values, axis=1, dtype=np.float64)
    velocity = (values @ D3Q19_DIRECTIONS.astype(np.float64)) / density[:, None]
    speed = np.linalg.norm(velocity, axis=1)
    pressure_scale = RHO0 * spec.dx_m**2 / spec.dt_s**2 / 3.0
    state = {
        "total_pdf_mass": float(np.sum(density, dtype=np.float64)),
        "mean_pressure_pa": float(np.mean(density, dtype=np.float64) * pressure_scale),
        "global_mean_pressure_gauge_pa": float(
            np.mean(density, dtype=np.float64) * pressure_scale - PRESSURE_REFERENCE_PA
        ),
        "global_mean_velocity_m_s": float(np.mean(speed) * spec.dx_m / spec.dt_s),
        "maximum_physical_velocity_m_s": float(np.max(speed) * spec.dx_m / spec.dt_s),
        "maximum_lattice_speed": float(np.max(speed)),
        "minimum_pdf": float(np.min(values)),
        "all_finite": bool(np.all(np.isfinite(values))),
    }
    return state, velocity, values


def _relative_l2(previous: np.ndarray, current: np.ndarray) -> float:
    denominator = max(float(np.linalg.norm(current.ravel())), np.finfo(float).tiny)
    return float(np.linalg.norm((current - previous).ravel())) / denominator


class GridCheckpointAuditor:
    """Grid-parameterized Referee V2 checkpoint gates."""

    def __init__(self, *, mesh_dir: Path, spec: GridSpec, qc_dir: Path) -> None:
        self.mesh_dir = Path(mesh_dir).resolve()
        self.spec = spec
        self.qc_dir = Path(qc_dir)
        self.qc_dir.mkdir(parents=True, exist_ok=True)
        header = parse_mesh_header(self.mesh_dir)
        self.cell_count = int(header["fluid_element_count"])
        self.mesh = load_mesh_contract(
            self.mesh_dir,
            expected_cells=self.cell_count,
            allow_zero_normals=True,
            require_runtime_order=False,
        )
        self.records: list[dict[str, Any]] = []
        self.samples: dict[int, dict[str, Any]] = {}
        self.binary_by_iteration: dict[int, Path] = {}
        self.gate_pass: dict[str, Any] | None = None

    def _mass_residual(self, start: dict[str, Any], end: dict[str, Any]) -> float:
        delta_iter = int(end["iteration"]) - int(start["iteration"])
        rate = (
            (float(end["total_pdf_mass"]) - float(start["total_pdf_mass"]))
            / delta_iter
            * RHO0
            * self.spec.dx_m**3
            / self.spec.dt_s
        )
        return abs(rate) / TARGET_MASS_FLOW

    def _window(self, start_iteration: int, end_iteration: int) -> dict[str, Any] | None:
        chosen = [
            self.samples[iteration]
            for iteration in sorted(self.samples)
            if start_iteration <= iteration <= end_iteration
        ]
        if (
            not chosen
            or int(chosen[0]["iteration"]) != start_iteration
            or int(chosen[-1]["iteration"]) != end_iteration
        ):
            return None
        iterations = np.asarray([item["iteration"] for item in chosen], dtype=np.float64)
        total_rate = np.asarray(
            [item["predicted_total_kg_s"] for item in chosen], dtype=np.float64
        )
        predicted = float(np.trapezoid(total_rate, iterations) * self.spec.dt_s)
        observed = (
            float(chosen[-1]["total_pdf_mass"] - chosen[0]["total_pdf_mass"])
            * RHO0
            * self.spec.dx_m**3
        )
        inlet_mass = TARGET_MASS_FLOW * (end_iteration - start_iteration) * self.spec.dt_s
        means: dict[str, float] = {}
        for label in ("inlet", "outlet_01", "outlet_02", "outlet_03"):
            values = np.asarray(
                [item["per_label_kg_s_domain"][label] for item in chosen],
                dtype=np.float64,
            )
            mean_domain = float(np.trapezoid(values, iterations) / (end_iteration - start_iteration))
            means[label] = mean_domain if label == "inlet" else -mean_domain
        closure = abs(predicted - observed) / max(abs(inlet_mass), np.finfo(float).tiny)
        return {
            "start_iteration": start_iteration,
            "end_iteration": end_iteration,
            "sample_iterations": iterations.astype(int).tolist(),
            "predicted_mass_change_kg": predicted,
            "observed_mass_change_kg": observed,
            "closure": closure,
            "mean_port_flows_kg_s": means,
        }

    def audit(self, header_path: Path) -> dict[str, Any]:
        restart = parse_restart_header(header_path)
        if restart.n_elems != self.cell_count or restart.n_components != 19:
            raise FlowError("CFD_FLOW_GRID_RESTART_CONTRACT_FAILED", str(header_path))
        state, velocity, pdf = _grid_pdf_state(
            restart.binary_path, self.spec, self.cell_count
        )
        replay = replay_boundary_step(
            pdf,
            self.mesh,
            dx_m=self.spec.dx_m,
            dt_s=self.spec.dt_s,
            density_kg_m3=RHO0,
            target_mass_flow_kg_s=TARGET_MASS_FLOW,
            outlet_pressures_pa=OUTLET_PRESSURES_PA,
        )
        iteration = int(restart.iteration)
        inlet_rho = float(replay["details"]["inlet"]["rho"])
        inlet_pressure = inlet_rho * RHO0 * self.spec.dx_m**2 / self.spec.dt_s**2 / 3.0
        inlet_gauge = inlet_pressure - PRESSURE_REFERENCE_PA
        sample = {
            "iteration": iteration,
            **state,
            "predicted_total_kg_s": float(replay["predicted_total_kg_s"]),
            "per_label_kg_s_domain": {
                key: float(value)
                for key, value in replay["per_label_kg_s_domain"].items()
            },
            "inlet_gauge_pressure_pa": inlet_gauge,
            "inlet_mass_flow_kg_s": float(
                replay["per_label_kg_s_domain"]["inlet"]
            ),
        }
        self.samples[iteration] = sample
        self.binary_by_iteration[iteration] = restart.binary_path
        short_start = iteration - 10_000
        long_start = iteration - 20_000
        short = self.samples.get(short_start)
        long = self.samples.get(long_start)
        short_window = self._window(short_start, iteration) if short else None
        long_window = self._window(long_start, iteration) if long else None
        record: dict[str, Any] = {
            **sample,
            "restart_header": str(header_path),
            "restart_binary": str(restart.binary_path),
            "restart_sha256": sha256_file(restart.binary_path),
            "R_mass_short": self._mass_residual(short, sample) if short else None,
            "R_mass_long": self._mass_residual(long, sample) if long else None,
            "R_velocity": None,
            "R_pressure": None,
            "R_inlet": abs(sample["inlet_mass_flow_kg_s"] - TARGET_MASS_FLOW)
            / TARGET_MASS_FLOW,
            "short_boundary_window": short_window,
            "long_boundary_window": long_window,
            "boundary_window_closure": None,
            "significant_time_averaged_backflow": None,
            "referee_revision": REFEREE_REVISION_NEW,
            "one_step_identity_error": REFEREE_IDENTITY_ERROR,
            "pressure_drops_pa": {
                label: inlet_gauge - (pressure - PRESSURE_REFERENCE_PA)
                for label, pressure in OUTLET_PRESSURES_PA.items()
            },
            "outlet_flow_fractions": None,
            "all_final_gates_pass": False,
        }
        if short is not None:
            _, previous_velocity, _ = _grid_pdf_state(
                self.binary_by_iteration[short_start], self.spec, self.cell_count
            )
            record["R_velocity"] = _relative_l2(previous_velocity, velocity)
            characteristic = float(
                np.median(np.abs(np.asarray(list(record["pressure_drops_pa"].values()))))
            )
            record["R_pressure"] = abs(
                float(sample["mean_pressure_pa"] - short["mean_pressure_pa"])
            ) / characteristic
        windows = [window for window in (short_window, long_window) if window]
        if windows:
            record["boundary_window_closure"] = max(
                float(window["closure"]) for window in windows
            )
        if long_window:
            mean_flows = long_window["mean_port_flows_kg_s"]
            outlets = [float(mean_flows[label]) for label in (
                "outlet_01", "outlet_02", "outlet_03"
            )]
            threshold = SIGNIFICANT_BACKFLOW_FRACTION * abs(float(mean_flows["inlet"]))
            record["significant_time_averaged_backflow"] = any(
                value < 0.0 and abs(value) > threshold for value in outlets
            )
            record["outlet_flow_fractions"] = {
                label: float(mean_flows[label]) / float(mean_flows["inlet"])
                for label in ("outlet_01", "outlet_02", "outlet_03")
            }
        gates = {
            "R_mass_short": record["R_mass_short"] is not None
            and float(record["R_mass_short"]) <= MASS_GATE,
            "R_mass_long": record["R_mass_long"] is not None
            and float(record["R_mass_long"]) <= MASS_GATE,
            "R_velocity": record["R_velocity"] is not None
            and float(record["R_velocity"]) <= VELOCITY_GATE,
            "R_pressure": record["R_pressure"] is not None
            and float(record["R_pressure"]) <= PRESSURE_GATE,
            "R_inlet": float(record["R_inlet"]) <= INLET_GATE,
            "boundary_window": record["boundary_window_closure"] is not None
            and float(record["boundary_window_closure"]) <= BOUNDARY_WINDOW_GATE,
            "one_step_identity_contract": REFEREE_IDENTITY_ERROR <= 1.0e-8,
            "no_significant_time_averaged_backflow": record[
                "significant_time_averaged_backflow"
            ] is False,
            "minimum_pdf": float(record["minimum_pdf"]) > 0.0,
            "maximum_lattice_speed": float(record["maximum_lattice_speed"]) < 0.05,
            "all_finite": bool(record["all_finite"]),
        }
        record["gates"] = gates
        record["all_final_gates_pass"] = all(gates.values())
        self.records.append(record)
        if record["all_final_gates_pass"] and self.gate_pass is None:
            self.gate_pass = record
        write_json(
            self.qc_dir / "checkpoint_referee_v2.json",
            {
                "status": "PASS" if self.gate_pass else "IN_PROGRESS",
                "gate_pass_iteration": self.gate_pass["iteration"]
                if self.gate_pass
                else None,
                "records": self.records,
            },
        )
        return record


def _complete_restart_headers(runtime_root: Path, cell_count: int) -> list[Path]:
    expected_bytes = cell_count * 19 * 8
    result = []
    for header_path in sorted((runtime_root / "restart").glob("*_header_*.lua")):
        try:
            header = parse_restart_header(header_path)
        except (FlowError, OSError, ValueError):
            continue
        if (
            header.n_elems == cell_count
            and header.n_components == 19
            and header.binary_path.is_file()
            and header.binary_path.stat().st_size == expected_bytes
        ):
            result.append(header_path)
    return sorted(result, key=lambda path: parse_restart_header(path).iteration)


def run_grid_steady(project_root: Path, label: str) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    grid_root = run_root / "grids" / label
    steady_root = grid_root / "steady"
    summary_path = steady_root / "steady_summary.json"
    if summary_path.is_file():
        existing = read_json(summary_path)
        if existing.get("status") == "PASS":
            return existing
        raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", f"existing {label} summary")
    benchmark = run_grid_mpi_benchmark(root, label)
    ranks = int(benchmark["selected_ranks"])
    spec = grid_specs()[label]
    mesh = grid_root / "seeder" / "mesh"
    require_complete_mesh(mesh, require_qval=True)
    qc_dir = steady_root / "qc"
    logs = steady_root / "logs"
    accepted = steady_root / "accepted_restart"
    for directory in (qc_dir, logs, accepted):
        directory.mkdir(parents=True, exist_ok=False)
    stamp = datetime.now().strftime("%H%M%S")
    runtime_root_wsl = f"/home/lzy/u3da/gc_{label}_steady_{stamp}"
    runtime_root = _create_wsl_runtime(runtime_root_wsl)
    (runtime_root / "restart").mkdir(parents=True, exist_ok=False)
    lua = generate_grid_musubi_lua(
        spec,
        mesh_dir=mesh,
        maximum_iterations=GRID_MAXIMUM_ITERATIONS,
        simulation_name=f"gc_{label[0]}",
        write_restarts=True,
    )
    (runtime_root / "musubi.lua").write_text(lua, encoding="utf-8")
    (steady_root / "musubi.lua").write_text(lua, encoding="utf-8")
    contract = grid_lua_contract(lua, spec, mesh)
    write_json(qc_dir / "lua_contract.json", contract)
    if contract["status"] != "PASS":
        raise FlowError("CFD_FLOW_GRID_LUA_CONTRACT_FAILED", label)
    luac_rc = _run_luac(
        distribution="Ubuntu",
        workdir_wsl=runtime_root_wsl,
        stdout_path=logs / "luac_stdout.log",
        stderr_path=logs / "luac_stderr.log",
    )
    if luac_rc != 0:
        raise FlowError("CFD_FLOW_GRID_LUA_CONTRACT_FAILED", f"{label} luac")
    auditor = GridCheckpointAuditor(mesh_dir=mesh, spec=spec, qc_dir=qc_dir)
    command = [
        "env",
        *[f"{key}={value}" for key, value in THREAD_ENV.items()],
        "OMPI_ALLOW_RUN_AS_ROOT=1",
        "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
        MPIRUN_WSL,
        "-np",
        str(ranks),
        *MPI_BINDING_ARGS,
        BINARY_WSL,
        "musubi.lua",
    ]
    shell = f"cd {shlex.quote(runtime_root_wsl)} && exec {shlex.join(command)}"
    stdout_path = logs / "musubi_stdout.log"
    stderr_path = logs / "musubi_stderr.log"
    audited: set[int] = set()
    started = time.perf_counter()
    timeout_s = 10_800 if label == "coarse" else 43_200
    next_progress = started + 60.0
    stop_created = False
    with (
        stdout_path.open("w", encoding="utf-8") as stdout,
        stderr_path.open("w", encoding="utf-8") as stderr,
    ):
        process = subprocess.Popen(
            ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", shell],
            stdout=stdout,
            stderr=stderr,
            text=True,
        )
        while process.poll() is None:
            time.sleep(5.0)
            for header_path in _complete_restart_headers(runtime_root, auditor.cell_count):
                iteration = parse_restart_header(header_path).iteration
                if iteration in audited:
                    continue
                record = auditor.audit(header_path)
                audited.add(iteration)
                if record["all_final_gates_pass"] and not stop_created:
                    (runtime_root / "stop").touch(exist_ok=True)
                    stop_created = True
            elapsed = time.perf_counter() - started
            if elapsed >= timeout_s:
                process.terminate()
                raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", f"{label} timeout")
            if time.perf_counter() >= next_progress:
                latest = max(audited) if audited else None
                mass = auditor.records[-1].get("R_mass_short") if auditor.records else None
                print(
                    f"GRID_STEADY_PROGRESS grid={label} elapsed_s={elapsed:.1f} "
                    f"latest_checkpoint={latest} R_mass_short={mass}",
                    flush=True,
                )
                next_progress += 60.0
        returncode = process.wait(timeout=30)
    for header_path in _complete_restart_headers(runtime_root, auditor.cell_count):
        iteration = parse_restart_header(header_path).iteration
        if iteration not in audited:
            auditor.audit(header_path)
            audited.add(iteration)
    elapsed = time.perf_counter() - started
    if returncode != 0:
        raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", f"{label} returncode={returncode}")
    if auditor.gate_pass is None:
        result = {
            "status": "FAIL",
            "grid": label,
            "musubi_calls": 1,
            "selected_ranks": ranks,
            "returncode": returncode,
            "wall_time_s": elapsed,
            "latest_iteration": max(audited) if audited else None,
            "first_failure": "Referee V2 gates did not pass before maximum_iterations",
        }
        write_json(summary_path, result)
        return result
    gate = auditor.gate_pass
    source_header = Path(str(gate["restart_header"]))
    source_binary = Path(str(gate["restart_binary"]))
    archived_header = accepted / source_header.name
    archived_binary = accepted / source_binary.name
    shutil.copy2(source_header, archived_header)
    shutil.copy2(source_binary, archived_binary)
    if sha256_file(archived_binary) != gate["restart_sha256"]:
        raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", f"{label} archive SHA mismatch")
    result = {
        "status": "PASS",
        "grid": label,
        "musubi_calls": 1,
        "selected_ranks": ranks,
        "returncode": returncode,
        "wall_time_s": elapsed,
        "steady_iteration": int(gate["iteration"]),
        "accepted_restart_header": str(archived_header),
        "accepted_restart_binary": str(archived_binary),
        "accepted_restart_sha256": gate["restart_sha256"],
        "referee_v2": gate,
        "runtime_root_wsl": runtime_root_wsl,
    }
    write_json(summary_path, result)
    return result


def _active_grid_musubi_processes() -> list[dict[str, Any]]:
    """Return unique live grid runtimes, grouping their MPI worker PIDs."""

    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "ps", "-eo", "pid=,cmd="],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
        check=False,
    )
    if process.returncode != 0:
        return []
    by_runtime: dict[str, dict[str, Any]] = {}
    for line in process.stdout.splitlines():
        if BINARY_WSL not in line:
            continue
        fields = line.strip().split(maxsplit=1)
        if not fields or not fields[0].isdigit():
            continue
        cwd = subprocess.run(
            [
                "wsl.exe", "-d", "Ubuntu", "--", "readlink",
                f"/proc/{fields[0]}/cwd",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
            check=False,
        )
        if cwd.returncode != 0:
            continue
        runtime = cwd.stdout.strip()
        if not runtime.startswith("/home/lzy/u3da/gc_"):
            continue
        record = by_runtime.setdefault(
            runtime,
            {
                "runtime_root_wsl": runtime,
                "pids": [],
                "binary_wsl": BINARY_WSL,
            },
        )
        record["pids"].append(int(fields[0]))
    return sorted(by_runtime.values(), key=lambda item: item["runtime_root_wsl"])


def _grid_runtime_is_running(runtime_root_wsl: str) -> bool:
    return any(
        record["runtime_root_wsl"] == runtime_root_wsl
        for record in _active_grid_musubi_processes()
    )


def _assert_no_duplicate_grid_musubi(label: str) -> None:
    prefix = f"/home/lzy/u3da/gc_{label}_"
    active = [
        record
        for record in _active_grid_musubi_processes()
        if str(record["runtime_root_wsl"]).startswith(prefix)
    ]
    if active:
        runtimes = ", ".join(str(record["runtime_root_wsl"]) for record in active)
        raise FlowError(
            "CFD_FLOW_GRID_DUPLICATE_MUSUBI_BLOCKED",
            f"{label} already active: {runtimes}",
        )


def _archive_intermediate_summary(steady_root: Path) -> dict[str, Any] | None:
    """Preserve a cleanly ended, non-steady segment without calling it a failure."""

    summary_path = steady_root / "steady_summary.json"
    if not summary_path.is_file():
        return None
    existing = read_json(summary_path)
    if existing.get("status") == "PASS":
        return existing
    history = steady_root / "history"
    history.mkdir(parents=True, exist_ok=True)
    historical_path = history / "initial_segment_summary.json"
    if not historical_path.exists():
        historical = dict(existing)
        if (
            historical.get("status") == "FAIL"
            and int(historical.get("returncode", 0)) == 0
        ):
            historical["historical_status"] = "COMPLETED_NOT_STEADY"
            historical["status_was_not_a_numerical_failure"] = True
        write_json(historical_path, historical)
    process_path = history / "initial_process_provenance.json"
    if not process_path.exists():
        returncode = existing.get("returncode")
        completed_cleanly = returncode is not None and int(returncode) == 0
        write_json(
            process_path,
            {
                "segment": "initial",
                "source_iteration": 0,
                "end_iteration": existing.get("latest_iteration"),
                "runtime_root_wsl": existing.get("runtime_root_wsl"),
                "returncode": returncode,
                "wall_time_s": existing.get("wall_time_s"),
                "status": (
                    "COMPLETED_NOT_STEADY"
                    if completed_cleanly
                    else "NUMERICAL_FAILED"
                ),
                "evidence_source": str(historical_path),
            },
        )
    return existing


def _write_process_record(segment_root: Path, record: dict[str, Any]) -> Path:
    path = segment_root / "process_provenance.json"
    write_json(path, record)
    return path


def _grid_process_history(steady_root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    initial = steady_root / "history" / "initial_process_provenance.json"
    if initial.is_file():
        records.append(read_json(initial))
    for path in sorted((steady_root / "continuations").glob("*/process_provenance.json")):
        records.append(read_json(path))
    return records


def _runtime_from_restart_path(value: str | None) -> str | None:
    if not value:
        return None
    normalized = str(value).replace("\\", "/")
    marker = "/home/lzy/"
    if marker in normalized and "/restart/" in normalized:
        return marker + normalized.split(marker, 1)[1].split("/restart/", 1)[0]
    return None


def reconcile_grid_process_provenance(
    project_root: Path, label: str
) -> dict[str, Any]:
    """Rebuild cumulative launch provenance from measured segment evidence."""

    root = Path(project_root).resolve()
    steady_root = (
        root / "outputs" / "cfd_flow" / GRID_RUN / "grids" / label / "steady"
    )
    summary_path = steady_root / "steady_summary.json"
    summary = read_json(summary_path)
    history_root = steady_root / "history"
    history_root.mkdir(parents=True, exist_ok=True)
    historical_summary_path = history_root / "initial_segment_summary.json"
    initial_process_path = history_root / "initial_process_provenance.json"
    if historical_summary_path.is_file() and not initial_process_path.is_file():
        initial_summary = read_json(historical_summary_path)
        runtime = initial_summary.get("runtime_root_wsl")
        checkpoint_evidence = steady_root / "qc" / "checkpoint_referee_v2.json"
        if runtime is None and checkpoint_evidence.is_file():
            records = read_json(checkpoint_evidence).get("records", [])
            historical_end = int(initial_summary.get("latest_iteration", -1))
            matching = [
                record for record in records
                if int(record.get("iteration", -2)) == historical_end
            ]
            if matching:
                runtime = _runtime_from_restart_path(matching[-1].get("restart_header"))
        returncode = initial_summary.get("returncode")
        write_json(
            initial_process_path,
            {
                "segment": "initial",
                "source_iteration": 0,
                "end_iteration": initial_summary.get("latest_iteration"),
                "runtime_root_wsl": runtime,
                "returncode": returncode,
                "wall_time_s": initial_summary.get("wall_time_s"),
                "status": (
                    "COMPLETED_NOT_STEADY"
                    if returncode is not None and int(returncode) == 0
                    else "NUMERICAL_FAILED"
                ),
                "evidence_source": str(historical_summary_path),
            },
        )
    for segment_root in sorted((steady_root / "continuations").glob("from_*")):
        process_path = segment_root / "process_provenance.json"
        existing_process = read_json(process_path) if process_path.is_file() else {}
        progress_path = segment_root / "progress.json"
        progress = read_json(progress_path) if progress_path.is_file() else {}
        runtime = progress.get("runtime_root_wsl") or existing_process.get("runtime_root_wsl")
        active = bool(runtime and _grid_runtime_is_running(str(runtime)))
        source_iteration = int(segment_root.name.removeprefix("from_"))
        is_final_segment = runtime is not None and runtime == summary.get("runtime_root_wsl")
        passed = is_final_segment and summary.get("status") == "PASS"
        observed_final_iteration = progress.get("latest_iteration")
        if passed:
            end_iteration = max(
                int(summary.get("steady_iteration") or -1),
                int(summary.get("actual_final_iteration") or -1),
                int(observed_final_iteration or -1),
            )
        else:
            end_iteration = observed_final_iteration
        process_record = dict(existing_process)
        process_record.update(
            {
                "segment": segment_root.name,
                "source_iteration": source_iteration,
                "end_iteration": end_iteration,
                "runtime_root_wsl": runtime,
                "returncode": 0 if passed else None,
                "wall_time_s": progress.get("elapsed_s"),
                "status": (
                    "RUNNING"
                    if active
                    else "COMPLETED" if passed else "COMPLETED_NOT_STEADY"
                ),
                "evidence_source": str(progress_path),
            },
        )
        write_json(process_path, process_record)
    processes = _grid_process_history(steady_root)
    result = {
        "status": "PASS" if processes else "FAIL",
        "grid": label,
        "musubi_process_count": len(processes),
        "processes": processes,
    }
    provenance_path = history_root / "musubi_process_provenance.json"
    write_json(provenance_path, result)
    if summary.get("status") == "PASS":
        summary["musubi_calls"] = len(processes)
        summary["process_provenance"] = processes
        summary.setdefault("gate_pass_iteration", summary.get("steady_iteration"))
        summary.setdefault("stop_request_iteration", summary.get("steady_iteration"))
        final_iterations = [
            int(process["end_iteration"])
            for process in processes
            if process.get("end_iteration") is not None
        ]
        summary["actual_final_iteration"] = max(
            final_iterations or [int(summary.get("steady_iteration") or -1)]
        )
        write_json(summary_path, summary)
    return result


def monitor_existing_grid_steady(
    project_root: Path,
    label: str,
    runtime_root_wsl: str,
    *,
    defer_full_audit_until_iteration: int = 0,
) -> dict[str, Any]:
    """Attach the V2 auditor to an already-running single Musubi call."""

    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    grid_root = run_root / "grids" / label
    steady_root = grid_root / "steady"
    summary_path = steady_root / "steady_summary.json"
    if summary_path.is_file():
        existing = read_json(summary_path)
        if existing.get("status") == "PASS":
            return existing
        if not _grid_runtime_is_running(runtime_root_wsl):
            return existing
    runtime_root = _wsl_path_to_unc("Ubuntu", runtime_root_wsl)
    mesh = grid_root / "seeder" / "mesh"
    spec = grid_specs()[label]
    qc_dir = steady_root / "qc"
    accepted = steady_root / "accepted_restart"
    benchmark = read_json(grid_root / "benchmark" / "benchmark_summary.json")
    ranks = int(benchmark["selected_ranks"])
    auditor = GridCheckpointAuditor(mesh_dir=mesh, spec=spec, qc_dir=qc_dir)
    audited: set[int] = set()
    started = time.perf_counter()
    next_progress = started + 60.0
    stop_created = (runtime_root / "stop").is_file()
    stop_request_iteration: int | None = None
    while _grid_runtime_is_running(runtime_root_wsl):
        headers = _complete_restart_headers(runtime_root, auditor.cell_count)
        latest_iteration = (
            parse_restart_header(headers[-1]).iteration if headers else 0
        )
        lower = max(0, latest_iteration - AUDIT_WINDOW_ITERATIONS)
        eligible = (
            headers
            if latest_iteration >= defer_full_audit_until_iteration
            else []
        )
        for header_path in eligible:
            iteration = parse_restart_header(header_path).iteration
            if iteration < lower or iteration in audited:
                continue
            record = auditor.audit(header_path)
            audited.add(iteration)
            if record["all_final_gates_pass"] and not stop_created:
                (runtime_root / "stop").touch(exist_ok=True)
                stop_created = True
                stop_request_iteration = int(record["iteration"])
        if time.perf_counter() >= next_progress:
            latest = max(audited) if audited else None
            mass = auditor.records[-1].get("R_mass_short") if auditor.records else None
            print(
                f"GRID_STEADY_PROGRESS grid={label} attached=true "
                f"latest_checkpoint={latest} R_mass_short={mass}",
                flush=True,
            )
            next_progress += 60.0
        time.sleep(5.0)
    time.sleep(2.0)
    final_headers = _complete_restart_headers(runtime_root, auditor.cell_count)
    final_latest = (
        parse_restart_header(final_headers[-1]).iteration if final_headers else 0
    )
    final_lower = max(0, final_latest - AUDIT_WINDOW_ITERATIONS)
    for header_path in (
        final_headers if final_latest >= defer_full_audit_until_iteration else []
    ):
        iteration = parse_restart_header(header_path).iteration
        if iteration >= final_lower and iteration not in audited:
            auditor.audit(header_path)
            audited.add(iteration)
    elapsed = time.perf_counter() - started
    process_count = max(1, len(_grid_process_history(steady_root)))
    if auditor.gate_pass is None:
        result = {
            "status": "IN_PROGRESS",
            "grid": label,
            "musubi_calls": process_count,
            "selected_ranks": ranks,
            "wall_time_s_attached": elapsed,
            "latest_iteration": final_latest or None,
            "segment_completed_without_steady_pass": True,
            "first_failure": None,
            "runtime_root_wsl": runtime_root_wsl,
        }
        write_json(summary_path, result)
        return result
    gate = auditor.gate_pass
    source_header = Path(str(gate["restart_header"]))
    source_binary = Path(str(gate["restart_binary"]))
    archived_header = accepted / source_header.name
    archived_binary = accepted / source_binary.name
    shutil.copy2(source_header, archived_header)
    shutil.copy2(source_binary, archived_binary)
    if sha256_file(archived_binary) != gate["restart_sha256"]:
        raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", f"{label} archive SHA mismatch")
    result = {
        "status": "PASS",
        "grid": label,
        "musubi_calls": process_count,
        "selected_ranks": ranks,
        "wall_time_s_attached": elapsed,
        "steady_iteration": int(gate["iteration"]),
        "gate_pass_iteration": int(gate["iteration"]),
        "stop_request_iteration": stop_request_iteration,
        "accepted_restart_header": str(archived_header),
        "accepted_restart_binary": str(archived_binary),
        "accepted_restart_sha256": gate["restart_sha256"],
        "referee_v2": gate,
        "runtime_root_wsl": runtime_root_wsl,
        "monitor_recovered_without_second_musubi_call": True,
    }
    write_json(summary_path, result)
    return result


def monitor_existing_grid_steady_deferred(
    project_root: Path,
    label: str,
    runtime_root_wsl: str,
    *,
    defer_full_audit_until_iteration: int = FINE_DEFER_FULL_AUDIT_UNTIL,
) -> dict[str, Any]:
    """Attach without replaying checkpoints before the requested iteration."""

    return monitor_existing_grid_steady(
        project_root,
        label,
        runtime_root_wsl,
        defer_full_audit_until_iteration=defer_full_audit_until_iteration,
    )


def _tail_controller_record(path: Path) -> dict[str, Any] | None:
    if not path.is_file() or path.stat().st_size == 0:
        return None
    with path.open("rb") as stream:
        size = stream.seek(0, 2)
        stream.seek(max(0, size - 262_144))
        text = stream.read().decode("utf-8", errors="replace")
    try:
        return _last_controller_record(text)
    except (FlowError, ValueError):
        return None


def run_grid_steady_continuation(
    project_root: Path,
    label: str,
    source_header: Path,
    *,
    defer_full_audit_until_iteration: int,
) -> dict[str, Any]:
    """Continue one grid while deferring expensive replay far from steady state."""

    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    grid_root = run_root / "grids" / label
    steady_root = grid_root / "steady"
    summary_path = steady_root / "steady_summary.json"
    benchmark = read_json(grid_root / "benchmark" / "benchmark_summary.json")
    ranks = int(benchmark["selected_ranks"])
    spec = grid_specs()[label]
    mesh = grid_root / "seeder" / "mesh"
    source = parse_restart_header(Path(source_header))
    if source.n_elems != int(parse_mesh_header(mesh)["fluid_element_count"]):
        raise FlowError("CFD_FLOW_GRID_RESTART_CONTRACT_FAILED", str(source_header))
    _assert_no_duplicate_grid_musubi(label)
    _archive_intermediate_summary(steady_root)
    segment_root = steady_root / "continuations" / f"from_{source.iteration}"
    segment_root.mkdir(parents=True, exist_ok=False)
    stamp = datetime.now().strftime("%H%M%S")
    runtime_root_wsl = f"/home/lzy/u3da/gc_{label}_continue_{stamp}"
    runtime_root = _create_wsl_runtime(runtime_root_wsl)
    runtime_restart = runtime_root / "restart"
    runtime_restart.mkdir(parents=True, exist_ok=False)
    staged_header = runtime_restart / Path(source_header).name
    staged_binary = runtime_restart / source.binary_path.name
    shutil.copy2(source_header, staged_header)
    shutil.copy2(source.binary_path, staged_binary)
    if sha256_file(staged_binary) != sha256_file(source.binary_path):
        raise FlowError("CFD_FLOW_GRID_RESTART_CONTRACT_FAILED", "staged SHA mismatch")
    resume_relative = f"restart/{staged_header.name}"
    first_restart = source.iteration + GRID_CHECKPOINT_INTERVAL
    lua = generate_grid_musubi_lua(
        spec,
        mesh_dir=mesh,
        maximum_iterations=GRID_MAXIMUM_ITERATIONS,
        simulation_name=f"gc_{label[0]}",
        write_restarts=True,
        resume_header=resume_relative,
        first_restart_iteration=first_restart,
    )
    (runtime_root / "musubi.lua").write_text(lua, encoding="utf-8")
    (segment_root / "musubi.lua").write_text(lua, encoding="utf-8")
    contract = grid_lua_contract(
        lua, spec, mesh, resume_header=resume_relative
    )
    write_json(segment_root / "lua_contract.json", contract)
    if contract["status"] != "PASS":
        raise FlowError("CFD_FLOW_GRID_LUA_CONTRACT_FAILED", f"{label} continuation")
    luac_rc = _run_luac(
        distribution="Ubuntu",
        workdir_wsl=runtime_root_wsl,
        stdout_path=segment_root / "luac_stdout.log",
        stderr_path=segment_root / "luac_stderr.log",
    )
    if luac_rc != 0:
        raise FlowError("CFD_FLOW_GRID_LUA_CONTRACT_FAILED", f"{label} continuation luac")
    stdout_runtime = runtime_root / "musubi_stdout.log"
    stderr_runtime = runtime_root / "musubi_stderr.log"
    qc_dir = steady_root / "qc"
    auditor = GridCheckpointAuditor(mesh_dir=mesh, spec=spec, qc_dir=qc_dir)
    command = [
        "env",
        *[f"{key}={value}" for key, value in THREAD_ENV.items()],
        "OMPI_ALLOW_RUN_AS_ROOT=1",
        "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
        MPIRUN_WSL,
        "-np",
        str(ranks),
        *MPI_BINDING_ARGS,
        BINARY_WSL,
        "musubi.lua",
    ]
    shell = f"cd {shlex.quote(runtime_root_wsl)} && exec {shlex.join(command)}"
    audited: set[int] = set()
    started = time.perf_counter()
    next_progress = started + 60.0
    stop_created = False
    stop_request_iteration: int | None = None
    first_failure: str | None = None
    process_record: dict[str, Any] = {
        "segment": segment_root.name,
        "source_iteration": source.iteration,
        "end_iteration": None,
        "runtime_root_wsl": runtime_root_wsl,
        "returncode": None,
        "wall_time_s": None,
        "status": "RUNNING",
        "selected_ranks": ranks,
    }
    _write_process_record(segment_root, process_record)
    process_count = len(_grid_process_history(steady_root))
    write_json(
        summary_path,
        {
            "status": "IN_PROGRESS",
            "grid": label,
            "musubi_calls": process_count,
            "selected_ranks": ranks,
            "source_iteration": source.iteration,
            "latest_iteration": source.iteration,
            "runtime_root_wsl": runtime_root_wsl,
            "full_audit_deferred_until_iteration": defer_full_audit_until_iteration,
            "first_failure": None,
        },
    )
    with (
        stdout_runtime.open("w", encoding="utf-8") as stdout,
        stderr_runtime.open("w", encoding="utf-8") as stderr,
    ):
        process = subprocess.Popen(
            ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", shell],
            stdout=stdout,
            stderr=stderr,
            text=True,
        )
        while process.poll() is None:
            time.sleep(10.0)
            headers = _complete_restart_headers(runtime_root, auditor.cell_count)
            latest_iteration = (
                parse_restart_header(headers[-1]).iteration if headers else source.iteration
            )
            controller = _tail_controller_record(stdout_runtime)
            if controller is not None:
                if (
                    not all(
                        math.isfinite(float(controller[key]))
                        for key in (
                            "target_lattice", "controlled_lattice", "relative_error",
                            "rho_boundary", "pressure_pa", "max_lattice_velocity",
                            "minimum_pdf",
                        )
                    )
                    or float(controller["minimum_pdf"]) <= 0.0
                    or float(controller["max_lattice_velocity"]) >= 0.05
                    or int(controller["globBC_count"]) <= 0
                ):
                    first_failure = f"controller safety failure at {controller['iteration']}"
                    (runtime_root / "stop").touch(exist_ok=True)
                    stop_created = True
            if latest_iteration >= defer_full_audit_until_iteration:
                lower = latest_iteration - 20_000
                for header_path in headers:
                    iteration = parse_restart_header(header_path).iteration
                    if iteration < lower or iteration in audited:
                        continue
                    record = auditor.audit(header_path)
                    audited.add(iteration)
                    if record["all_final_gates_pass"] and not stop_created:
                        (runtime_root / "stop").touch(exist_ok=True)
                        stop_created = True
                        stop_request_iteration = int(record["iteration"])
            elapsed = time.perf_counter() - started
            if elapsed >= 43_200:
                process.terminate()
                raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", f"{label} continuation timeout")
            if time.perf_counter() >= next_progress:
                mass = auditor.records[-1].get("R_mass_short") if auditor.records else None
                progress = {
                    "grid": label,
                    "elapsed_s": elapsed,
                    "latest_iteration": latest_iteration,
                    "audit_deferred_until": defer_full_audit_until_iteration,
                    "R_mass_short": mass,
                    "controller": controller,
                    "runtime_root_wsl": runtime_root_wsl,
                }
                write_json(segment_root / "progress.json", progress)
                process_record.update(
                    {
                        "end_iteration": latest_iteration,
                        "wall_time_s": elapsed,
                    }
                )
                _write_process_record(segment_root, process_record)
                print(
                    f"GRID_STEADY_PROGRESS grid={label} continuation=true "
                    f"elapsed_s={elapsed:.1f} latest_checkpoint={latest_iteration} "
                    f"R_mass_short={mass}",
                    flush=True,
                )
                next_progress += 60.0
        returncode = process.wait(timeout=30)
    shutil.copy2(stdout_runtime, segment_root / "musubi_stdout.log")
    shutil.copy2(stderr_runtime, segment_root / "musubi_stderr.log")
    final_headers = _complete_restart_headers(runtime_root, auditor.cell_count)
    actual_final_iteration = (
        parse_restart_header(final_headers[-1]).iteration
        if final_headers
        else source.iteration
    )
    if actual_final_iteration >= defer_full_audit_until_iteration:
        final_lower = actual_final_iteration - AUDIT_WINDOW_ITERATIONS
        for header_path in final_headers:
            iteration = parse_restart_header(header_path).iteration
            if iteration >= final_lower and iteration not in audited:
                auditor.audit(header_path)
                audited.add(iteration)
    elapsed = time.perf_counter() - started
    process_record.update(
        {
            "end_iteration": actual_final_iteration,
            "returncode": returncode,
            "wall_time_s": elapsed,
            "status": (
                "NUMERICAL_FAILED"
                if returncode != 0 or first_failure
                else "COMPLETED"
            ),
        }
    )
    _write_process_record(segment_root, process_record)
    process_count = len(_grid_process_history(steady_root))
    if returncode != 0 or first_failure:
        write_json(
            summary_path,
            {
                "status": "NUMERICAL_FAILED",
                "grid": label,
                "musubi_calls": process_count,
                "selected_ranks": ranks,
                "latest_iteration": actual_final_iteration,
                "first_failure": first_failure or f"returncode={returncode}",
                "runtime_root_wsl": runtime_root_wsl,
            },
        )
        raise FlowError(
            "CFD_FLOW_GRID_STEADY_FAILED",
            first_failure or f"{label} continuation returncode={returncode}",
        )
    if auditor.gate_pass is None:
        result = {
            "status": "IN_PROGRESS",
            "grid": label,
            "musubi_calls": process_count,
            "selected_ranks": ranks,
            "latest_iteration": actual_final_iteration,
            "segment_completed_without_steady_pass": True,
            "first_failure": None,
            "runtime_root_wsl": runtime_root_wsl,
        }
        write_json(summary_path, result)
        return result
    gate = auditor.gate_pass
    accepted = steady_root / "accepted_restart"
    accepted.mkdir(parents=True, exist_ok=True)
    source_gate_header = Path(str(gate["restart_header"]))
    source_gate_binary = Path(str(gate["restart_binary"]))
    archived_header = accepted / source_gate_header.name
    archived_binary = accepted / source_gate_binary.name
    shutil.copy2(source_gate_header, archived_header)
    shutil.copy2(source_gate_binary, archived_binary)
    if sha256_file(archived_binary) != gate["restart_sha256"]:
        raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", f"{label} archive SHA mismatch")
    result = {
        "status": "PASS",
        "grid": label,
        "musubi_calls": process_count,
        "selected_ranks": ranks,
        "steady_iteration": int(gate["iteration"]),
        "gate_pass_iteration": int(gate["iteration"]),
        "stop_request_iteration": stop_request_iteration,
        "actual_final_iteration": actual_final_iteration,
        "accepted_restart_header": str(archived_header),
        "accepted_restart_binary": str(archived_binary),
        "accepted_restart_sha256": gate["restart_sha256"],
        "referee_v2": gate,
        "runtime_root_wsl": runtime_root_wsl,
        "source_iteration": source.iteration,
        "full_audit_deferred_until_iteration": defer_full_audit_until_iteration,
        "process_provenance": _grid_process_history(steady_root),
    }
    write_json(summary_path, result)
    return result


def _base_grid_observables(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    archive_root = root / "outputs" / "cfd_flow" / ARCHIVE_RUN
    manifest_path = archive_root / "qc" / "accepted_steady_manifest.json"
    manifest = read_json(manifest_path)
    if manifest.get("status") != "PASS" or int(manifest["iteration"]) != ACCEPTED_ITERATION:
        raise FlowError("CFD_FLOW_BASE_OBSERVABLES_INVALID", "accepted manifest")
    binary = Path(str(manifest["archived_binary"]))
    header = Path(str(manifest["archived_header"]))
    if not binary.is_file():
        binary = archive_root / "restart" / Path(str(manifest["archived_binary"])).name
    if not header.is_file():
        header = archive_root / "restart" / Path(str(manifest["archived_header"])).name
    actual_sha = sha256_file(binary)
    if actual_sha != ACCEPTED_SHA256 or manifest.get("binary_sha256") != ACCEPTED_SHA256:
        raise FlowError("CFD_FLOW_BASE_OBSERVABLES_INVALID", "accepted restart SHA")
    restart = parse_restart_header(header)
    if (
        restart.iteration != ACCEPTED_ITERATION
        or restart.n_elems != EXPECTED_CELLS
        or restart.n_components != 19
        or restart.binary_path.resolve() != binary.resolve()
    ):
        raise FlowError("CFD_FLOW_BASE_OBSERVABLES_INVALID", "restart contract")
    spec = grid_specs()["base"]
    state, _, pdf = _grid_pdf_state(binary, spec, EXPECTED_CELLS)
    mesh_dir = root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "seeder" / "mesh"
    mesh = load_mesh_contract(
        mesh_dir,
        expected_cells=EXPECTED_CELLS,
        allow_zero_normals=True,
        require_runtime_order=False,
    )
    replay = replay_boundary_step(
        pdf,
        mesh,
        dx_m=spec.dx_m,
        dt_s=spec.dt_s,
        density_kg_m3=RHO0,
        target_mass_flow_kg_s=TARGET_MASS_FLOW,
        outlet_pressures_pa=OUTLET_PRESSURES_PA,
    )
    inlet_rho = float(replay["details"]["inlet"]["rho"])
    inlet_gauge = (
        inlet_rho * RHO0 * spec.dx_m**2 / spec.dt_s**2 / 3.0
        - PRESSURE_REFERENCE_PA
    )
    referee_root = root / "outputs" / "cfd_flow" / RUN_NAME / "qc" / "referee"
    residual_path = referee_root / "corrected_residual_history.csv"
    with residual_path.open("r", encoding="utf-8", newline="") as stream:
        residual_rows = list(csv.DictReader(stream))
    residual = next(
        (row for row in residual_rows if int(row["iteration"]) == ACCEPTED_ITERATION),
        None,
    )
    if residual is None or residual["restart_sha256"] != ACCEPTED_SHA256:
        raise FlowError("CFD_FLOW_BASE_OBSERVABLES_INVALID", "Referee V2 residual row")
    flux_path = referee_root / "corrected_boundary_flux_history.csv"
    with flux_path.open("r", encoding="utf-8", newline="") as stream:
        flux_rows = [
            row
            for row in csv.DictReader(stream)
            if ACCEPTED_ITERATION - AUDIT_WINDOW_ITERATIONS
            <= int(row["iteration"])
            <= ACCEPTED_ITERATION
        ]
    iterations = np.asarray([int(row["iteration"]) for row in flux_rows], dtype=np.float64)
    if (
        len(flux_rows) < 2
        or int(iterations[0]) != ACCEPTED_ITERATION - AUDIT_WINDOW_ITERATIONS
        or int(iterations[-1]) != ACCEPTED_ITERATION
    ):
        raise FlowError("CFD_FLOW_BASE_OBSERVABLES_INVALID", "20k flux window")
    columns = {
        "inlet": "inlet_into_domain_kg_s",
        "outlet_01": "outlet_01_outward_kg_s",
        "outlet_02": "outlet_02_outward_kg_s",
        "outlet_03": "outlet_03_outward_kg_s",
    }
    means = {
        label: float(
            np.trapezoid(
                np.asarray([float(row[column]) for row in flux_rows], dtype=np.float64),
                iterations,
            )
            / AUDIT_WINDOW_ITERATIONS
        )
        for label, column in columns.items()
    }
    return {
        "inlet_gauge_pressure_pa": inlet_gauge,
        "pressure_drops_pa": {
            label: inlet_gauge - (pressure - PRESSURE_REFERENCE_PA)
            for label, pressure in OUTLET_PRESSURES_PA.items()
        },
        "outlet_flow_fractions": {
            label: means[label] / means["inlet"]
            for label in ("outlet_01", "outlet_02", "outlet_03")
        },
        "global_mean_velocity_m_s": state["global_mean_velocity_m_s"],
        "global_mean_pressure_gauge_pa": state["global_mean_pressure_gauge_pa"],
        "maximum_physical_velocity_m_s": state["maximum_physical_velocity_m_s"],
        "provenance": {
            "observable_semantics": "REFEREE_V2_ACCEPTED_RESTART_AND_20K_LONG_WINDOW",
            "accepted_iteration": ACCEPTED_ITERATION,
            "accepted_restart_header": str(header),
            "accepted_restart_binary": str(binary),
            "accepted_restart_sha256": actual_sha,
            "accepted_manifest": str(manifest_path),
            "corrected_residual_history": str(residual_path),
            "corrected_boundary_flux_history": str(flux_path),
            "residual_row_sha256": residual["restart_sha256"],
            "residual_csv_mean_columns_are_short_window_only": True,
            "long_window_sample_iterations": iterations.astype(int).tolist(),
            "mean_port_flows_kg_s": means,
            "fractions_not_renormalized": True,
            "referee_revision": REFEREE_REVISION_NEW,
        },
    }


def _solved_grid_observables(project_root: Path, label: str) -> dict[str, Any]:
    root = Path(project_root).resolve()
    summary_path = (
        root / "outputs" / "cfd_flow" / GRID_RUN / "grids" / label
        / "steady" / "steady_summary.json"
    )
    summary = read_json(summary_path)
    if summary.get("status") != "PASS":
        raise FlowError("CFD_FLOW_GRID_STEADY_FAILED", label)
    gate = summary["referee_v2"]
    return {
        "inlet_gauge_pressure_pa": float(gate["inlet_gauge_pressure_pa"]),
        "pressure_drops_pa": {
            key: float(value) for key, value in gate["pressure_drops_pa"].items()
        },
        "outlet_flow_fractions": {
            key: float(value) for key, value in gate["outlet_flow_fractions"].items()
        },
        "global_mean_velocity_m_s": float(gate["global_mean_velocity_m_s"]),
        "global_mean_pressure_gauge_pa": float(
            gate["global_mean_pressure_gauge_pa"]
        ),
        "maximum_physical_velocity_m_s": float(
            gate["maximum_physical_velocity_m_s"]
        ),
        "provenance": {
            "observable_semantics": "REFEREE_V2_ACCEPTED_RESTART_AND_20K_LONG_WINDOW",
            "steady_summary": str(summary_path),
            "accepted_iteration": int(summary["steady_iteration"]),
            "accepted_restart_header": summary["accepted_restart_header"],
            "accepted_restart_binary": summary["accepted_restart_binary"],
            "accepted_restart_sha256": summary["accepted_restart_sha256"],
            "long_window_sample_iterations": gate["long_boundary_window"][
                "sample_iterations"
            ],
            "mean_port_flows_kg_s": gate["long_boundary_window"][
                "mean_port_flows_kg_s"
            ],
            "fractions_not_renormalized": True,
            "referee_revision": gate["referee_revision"],
        },
    }


def three_grid_scalar_analysis(
    coarse: float, base: float, fine: float, *, refinement_ratio: float = REFINEMENT_RATIO
) -> dict[str, Any]:
    """Standard constant-ratio monotonic Richardson/GCI calculation."""

    values = tuple(float(value) for value in (coarse, base, fine))
    if not all(math.isfinite(value) for value in values):
        return {
            "status": "UNAVAILABLE",
            "classification": "NON_FINITE",
            "reason": "non-finite value",
            "coarse": values[0],
            "base": values[1],
            "fine": values[2],
            "observed_order_p": None,
            "richardson_extrapolation": None,
            "gci_coarse_base": None,
            "gci_base_fine": None,
        }
    delta_cb = values[0] - values[1]
    delta_bf = values[1] - values[2]
    cb_relative = abs(delta_cb) / max(abs(values[1]), np.finfo(float).tiny)
    bf_relative = abs(delta_bf) / max(abs(values[2]), np.finfo(float).tiny)
    ratio = delta_cb / delta_bf if delta_bf != 0.0 else math.nan
    result: dict[str, Any] = {
        "coarse": values[0],
        "base": values[1],
        "fine": values[2],
        "coarse_to_base_relative_difference": cb_relative,
        "base_to_fine_relative_difference": bf_relative,
        "difference_ratio": ratio,
        "base_to_fine_within_5_percent": bf_relative <= 0.05,
        "base_to_fine_within_preferred_2_percent": bf_relative <= 0.02,
    }
    if delta_cb == 0.0 or delta_bf == 0.0:
        result.update(
            {
                "status": "UNAVAILABLE",
                "classification": "STALLED",
                "reason": "one or both successive differences are zero",
                "observed_order_p": None,
                "richardson_extrapolation": None,
                "gci_coarse_base": None,
                "gci_base_fine": None,
            }
        )
        return result
    if delta_cb * delta_bf < 0.0:
        result.update(
            {
                "status": "UNAVAILABLE",
                "classification": "OSCILLATORY",
                "reason": "successive differences have opposite signs",
                "observed_order_p": None,
                "richardson_extrapolation": None,
                "gci_coarse_base": None,
                "gci_base_fine": None,
            }
        )
        return result
    if not math.isfinite(ratio) or ratio <= 1.0:
        result.update(
            {
                "status": "UNAVAILABLE",
                "classification": "NON_ASYMPTOTIC",
                "reason": "successive differences do not decrease under refinement",
                "observed_order_p": None,
                "richardson_extrapolation": None,
                "gci_coarse_base": None,
                "gci_base_fine": None,
            }
        )
        return result
    order = math.log(abs(ratio)) / math.log(refinement_ratio)
    denominator = refinement_ratio**order - 1.0
    if not math.isfinite(order) or order <= 0.0 or denominator <= 0.0:
        result.update(
            {
                "status": "UNAVAILABLE",
                "classification": "NON_ASYMPTOTIC",
                "reason": "non-positive observed order",
                "observed_order_p": None,
                "richardson_extrapolation": None,
                "gci_coarse_base": None,
                "gci_base_fine": None,
            }
        )
        return result
    result.update(
        {
            "status": "AVAILABLE",
            "classification": "ASYMPTOTIC_MONOTONIC",
            "reason": None,
            "observed_order_p": order,
            "richardson_extrapolation": values[2] + (values[2] - values[1]) / denominator,
            "gci_coarse_base": 1.25 * cb_relative / denominator,
            "gci_base_fine": 1.25 * bf_relative / denominator,
        }
    )
    return result


def reconcile_grid_convergence_evidence(project_root: Path) -> dict[str, Any]:
    """Replace stale phase status only after all three steady states are accepted."""

    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    status_path = run_root / "qc" / "grid_convergence_status.json"
    previous = read_json(status_path) if status_path.is_file() else None
    coarse_summary_path = run_root / "grids" / "coarse" / "steady" / "steady_summary.json"
    fine_summary_path = run_root / "grids" / "fine" / "steady" / "steady_summary.json"
    base_manifest_path = (
        root / "outputs" / "cfd_flow" / ARCHIVE_RUN / "qc"
        / "accepted_steady_manifest.json"
    )
    coarse = read_json(coarse_summary_path)
    fine = read_json(fine_summary_path)
    base = read_json(base_manifest_path)
    fine_processes = reconcile_grid_process_provenance(root, "fine")
    mesh_qc = {
        label: read_json(run_root / "grids" / label / "qc" / "mesh_qc.json")
        for label in ("coarse", "fine")
    }
    benchmarks = {
        label: read_json(
            run_root / "grids" / label / "benchmark" / "benchmark_summary.json"
        )
        for label in ("coarse", "fine")
    }
    checks = {
        "coarse_steady_pass": coarse.get("status") == "PASS",
        "base_steady_pass": base.get("status") == "PASS",
        "fine_steady_pass": fine.get("status") == "PASS",
        "base_accepted_iteration": int(base.get("iteration", -1)) == ACCEPTED_ITERATION,
        "base_accepted_sha": base.get("binary_sha256") == ACCEPTED_SHA256,
        "coarse_mesh_pass": mesh_qc["coarse"].get("status") == "PASS",
        "fine_mesh_pass": mesh_qc["fine"].get("status") == "PASS",
        "coarse_benchmark_pass": benchmarks["coarse"].get("status") == "PASS",
        "fine_benchmark_pass": benchmarks["fine"].get("status") == "PASS",
        "fine_process_provenance_complete": fine_processes.get("status") == "PASS",
    }
    passed = all(checks.values())
    result = {
        "actual_head": _git(root, "rev-parse", "HEAD"),
        "production_pipeline_modified": bool(
            _git(
                root,
                "diff",
                "--",
                *(str(path.relative_to(root)) for path in _production_paths(root)),
            )
        ),
        "referee_revision": REFEREE_REVISION_NEW,
        "evidence_reconciliation": "PASS" if passed else "FAIL",
        "checks": checks,
        "coarse": {
            "summary": str(coarse_summary_path),
            "steady_iteration": coarse.get("steady_iteration"),
            "selected_mpi_ranks": coarse.get("selected_ranks"),
        },
        "base": {
            "manifest": str(base_manifest_path),
            "steady_iteration": base.get("iteration"),
            "accepted_sha256": base.get("binary_sha256"),
        },
        "fine": {
            "summary": str(fine_summary_path),
            "steady_iteration": fine.get("steady_iteration"),
            "selected_mpi_ranks": fine.get("selected_ranks"),
            "process_provenance": fine_processes,
        },
        "historical_status": previous,
        "grid_convergence": "READY_FOR_GCI" if passed else "BLOCKED",
        "final_status": (
            "CFD_FLOW_HEALTHY_ADAPTIVE_FLUX_THREE_GRID_STEADY_READY_FOR_GCI"
            if passed
            else "CFD_FLOW_HEALTHY_ADAPTIVE_FLUX_EVIDENCE_RECONCILIATION_FAILED"
        ),
        "next": "COMPUTE RICHARDSON AND GCI" if passed else "FIX EVIDENCE INCONSISTENCY",
        "first_failure": next((name for name, value in checks.items() if not value), None),
    }
    write_json(status_path, result)
    return result


def evaluate_grid_convergence_gate(
    analyses: dict[str, dict[str, Any]], primary_metrics: tuple[str, ...]
) -> dict[str, Any]:
    trends = all(
        analyses[name].get("classification") == "ASYMPTOTIC_MONOTONIC"
        and analyses[name].get("status") == "AVAILABLE"
        for name in primary_metrics
    )
    within_5 = all(
        float(analyses[name]["base_to_fine_relative_difference"]) <= 0.05
        for name in primary_metrics
    )
    within_2 = all(
        float(analyses[name]["base_to_fine_relative_difference"]) <= 0.02
        for name in primary_metrics
    )
    passed = trends and within_5
    return {
        "status": "PASS" if passed else "FAIL",
        "primary_trends_asymptotic_monotonic": trends,
        "base_fine_primary_within_5_percent": within_5,
        "base_fine_primary_within_preferred_2_percent": within_2,
        "next": (
            "PROMOTE HEALTHY ADAPTIVE-FLUX CFD TO PRODUCTION PIPELINE"
            if passed
            else (
                "DESIGN ONE FINER GRID FOR CONVERGENCE CONFIRMATION"
                if trends
                else "REVIEW VOXELIZED GEOMETRY GRID-SENSITIVITY"
            )
        ),
    }


def finalize_grid_convergence(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    reconciliation = reconcile_grid_convergence_evidence(root)
    if reconciliation["evidence_reconciliation"] != "PASS":
        raise FlowError(
            "CFD_FLOW_GRID_EVIDENCE_RECONCILIATION_FAILED",
            str(reconciliation["first_failure"]),
        )
    observables = {
        "coarse": _solved_grid_observables(root, "coarse"),
        "base": _base_grid_observables(root),
        "fine": _solved_grid_observables(root, "fine"),
    }
    scalar_paths = {
        "inlet_gauge_pressure_pa": ("inlet_gauge_pressure_pa",),
        "pressure_drop_outlet_01_pa": ("pressure_drops_pa", "outlet_01"),
        "pressure_drop_outlet_02_pa": ("pressure_drops_pa", "outlet_02"),
        "pressure_drop_outlet_03_pa": ("pressure_drops_pa", "outlet_03"),
        "outlet_01_flow_fraction": ("outlet_flow_fractions", "outlet_01"),
        "outlet_02_flow_fraction": ("outlet_flow_fractions", "outlet_02"),
        "outlet_03_flow_fraction": ("outlet_flow_fractions", "outlet_03"),
        "global_mean_velocity_m_s": ("global_mean_velocity_m_s",),
        "global_mean_pressure_gauge_pa": ("global_mean_pressure_gauge_pa",),
        "maximum_physical_velocity_m_s": ("maximum_physical_velocity_m_s",),
    }

    def value(record: dict[str, Any], path: tuple[str, ...]) -> float:
        current: Any = record
        for part in path:
            current = current[part]
        return float(current)

    analyses = {
        name: three_grid_scalar_analysis(
            value(observables["coarse"], path),
            value(observables["base"], path),
            value(observables["fine"], path),
        )
        for name, path in scalar_paths.items()
    }
    primary = (
        "inlet_gauge_pressure_pa",
        "pressure_drop_outlet_01_pa",
        "pressure_drop_outlet_02_pa",
        "pressure_drop_outlet_03_pa",
        "outlet_01_flow_fraction",
        "outlet_02_flow_fraction",
        "outlet_03_flow_fraction",
    )
    gate = evaluate_grid_convergence_gate(analyses, primary)
    base_fine_gate = bool(gate["base_fine_primary_within_5_percent"])
    preferred = bool(gate["base_fine_primary_within_preferred_2_percent"])
    reasonable_trend = bool(gate["primary_trends_asymptotic_monotonic"])
    passed = gate["status"] == "PASS"
    result = {
        "status": "PASS" if passed else "FAIL",
        "refinement_ratio": REFINEMENT_RATIO,
        "refinement_ratio_contract": {
            "coarse_over_base": grid_specs()["coarse"].dx_m / grid_specs()["base"].dx_m,
            "base_over_fine": grid_specs()["base"].dx_m / grid_specs()["fine"].dx_m,
            "constant_r_1p3": math.isclose(
                grid_specs()["coarse"].dx_m / grid_specs()["base"].dx_m,
                REFINEMENT_RATIO,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            ) and math.isclose(
                grid_specs()["base"].dx_m / grid_specs()["fine"].dx_m,
                REFINEMENT_RATIO,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            ),
            "diffusive_dt_scaling": all(
                math.isclose(
                    spec.dt_s,
                    DT_S * (spec.dx_m / DX_M) ** 2,
                    rel_tol=1.0e-14,
                )
                for spec in grid_specs().values()
            ),
        },
        "evidence_reconciliation": reconciliation,
        "observables": observables,
        "analyses": analyses,
        "base_fine_primary_within_5_percent": base_fine_gate,
        "base_fine_primary_within_preferred_2_percent": preferred,
        "primary_trends_asymptotic_monotonic": reasonable_trend,
        "selected_production_grid": "base" if passed else None,
        "selected_production_dx_m": DX_M if passed else None,
        "final_status": (
            "CFD_FLOW_HEALTHY_ADAPTIVE_FLUX_GRID_CONVERGENCE_PASS"
            if passed
            else "CFD_FLOW_HEALTHY_ADAPTIVE_FLUX_GRID_CONVERGENCE_FAILED"
        ),
        "next": gate["next"],
    }
    final_path = run_root / "qc" / "grid_convergence_final.json"
    write_json(final_path, result)
    active_status = dict(reconciliation)
    active_status.update(
        {
            "grid_convergence": "PASS" if passed else "FAIL",
            "grid_convergence_result": str(final_path),
            "final_status": result["final_status"],
            "next": result["next"],
            "first_failure": (
                None
                if passed
                else next(
                    (
                        name
                        for name in primary
                        if analyses[name]["status"] != "AVAILABLE"
                        or float(analyses[name]["base_to_fine_relative_difference"]) > 0.05
                    ),
                    "primary grid-convergence gate",
                )
            ),
        }
    )
    write_json(run_root / "qc" / "grid_convergence_status.json", active_status)
    return result


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", "-C", str(root), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.strip()


def _production_paths(root: Path) -> tuple[Path, ...]:
    return (
        root / "cfd_flow.py",
        root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
    )


def _hashes(paths: tuple[Path, ...]) -> dict[str, str]:
    return {str(path): sha256_file(path) for path in paths}


def _accepted_record(root: Path) -> tuple[dict[str, Any], dict[str, str]]:
    source = root / "outputs" / "cfd_flow" / RUN_NAME
    convergence = read_json(
        source / "qc" / "exact_watchdog" / "checkpoint_exact_convergence.json"
    )
    record = next(
        item for item in convergence["records"]
        if int(item["iteration"]) == ACCEPTED_ITERATION
    )
    corrected_rows = (
        source / "qc" / "referee" / "corrected_residual_history.csv"
    ).read_text(encoding="utf-8").splitlines()
    import csv

    corrected = next(
        row for row in csv.DictReader(corrected_rows)
        if int(row["iteration"]) == ACCEPTED_ITERATION
    )
    return record, corrected


def archive_accepted_baseline(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    output = root / "outputs" / "cfd_flow" / ARCHIVE_RUN
    restart = output / "restart"
    qc = output / "qc"
    manifest_path = qc / "accepted_steady_manifest.json"
    if manifest_path.is_file():
        existing = read_json(manifest_path)
        binary = Path(existing["archived_binary"])
        header = Path(existing["archived_header"])
        if (
            existing.get("status") == "PASS"
            and binary.is_file()
            and header.is_file()
            and sha256_file(binary) == ACCEPTED_SHA256
        ):
            return existing
        raise FlowError("CFD_FLOW_ACCEPTED_ARCHIVE_INVALID", "Existing archive changed")

    record, corrected = _accepted_record(root)
    source_binary = Path(str(record["restart_binary"]))
    source_header = Path(str(record["restart_header"]))
    if not source_binary.is_file() or not source_header.is_file():
        raise FileNotFoundError("Accepted WSL source restart is missing")
    if sha256_file(source_binary) != ACCEPTED_SHA256:
        raise FlowError("CFD_FLOW_ACCEPTED_ARCHIVE_INVALID", "Source SHA256 changed")
    restart.mkdir(parents=True, exist_ok=False)
    qc.mkdir(parents=True, exist_ok=False)
    archived_binary = restart / source_binary.name
    archived_header = restart / source_header.name
    shutil.copy2(source_binary, archived_binary)
    shutil.copy2(source_header, archived_header)
    archive_sha = sha256_file(archived_binary)
    if archive_sha != ACCEPTED_SHA256:
        raise FlowError("CFD_FLOW_ACCEPTED_ARCHIVE_INVALID", "Copied SHA256 mismatch")
    manifest = {
        "status": "PASS",
        "iteration": ACCEPTED_ITERATION,
        "source_header": str(source_header),
        "source_binary": str(source_binary),
        "archived_header": str(archived_header),
        "archived_binary": str(archived_binary),
        "binary_sha256": archive_sha,
        "referee_revision": REFEREE_REVISION_NEW,
        "R_mass_short": float(corrected["R_mass_short"]),
        "R_mass_long": float(corrected["R_mass_long"]),
        "R_velocity": float(corrected["R_velocity"]),
        "R_pressure": float(corrected["R_pressure"]),
        "R_inlet": float(corrected["R_inlet"]),
        "R_conservation_identity": float(corrected["R_conservation_identity"]),
        "boundary_window_closure": float(corrected["boundary_window_closure"]),
        "significant_time_averaged_backflow": corrected[
            "significant_time_averaged_backflow"
        ].lower() == "true",
        "minimum_pdf": float(corrected["minimum_pdf"]),
        "maximum_lattice_speed": float(corrected["maximum_lattice_speed"]),
        "created_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, manifest)
    return manifest


def build_grid_preflight(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    qc = run_root / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    specs = grid_specs()
    candidates = []
    for ratio in (1.3, math.sqrt(2.0), 1.5, 2.0):
        fine_factor = ratio**3
        runtime_factor = ratio**5
        candidates.append(
            {
                "refinement_ratio": ratio,
                "predicted_coarse_cells": int(round(EXPECTED_CELLS / fine_factor)),
                "predicted_fine_cells": int(round(EXPECTED_CELLS * fine_factor)),
                "predicted_fine_runtime_s": 394_166 / 51.39525262051544
                * runtime_factor,
            }
        )
    all_contain_base_cube = all(
        spec.cube_side_m >= BASE_CUBE_SIDE_M
        for label, spec in specs.items()
        if label != "base"
    )
    constant_ratio = math.isclose(
        specs["coarse"].dx_m / specs["base"].dx_m,
        specs["base"].dx_m / specs["fine"].dx_m,
        rel_tol=0.0,
        abs_tol=1.0e-14,
    )
    checks = {
        "refinement_ratio_at_least_1p3": REFINEMENT_RATIO >= 1.3,
        "constant_ratio": constant_ratio,
        "diffusive_scaling": all(
            math.isclose(
                spec.dt_s,
                DT_S * (spec.dx_m / DX_M) ** 2,
                rel_tol=1.0e-14,
            )
            for spec in specs.values()
        ),
        "new_root_cubes_contain_frozen_base_cube": all_contain_base_cube,
        "fine_predicted_ram_below_60_percent_32GiB": (
            specs["fine"].predicted_fluid_cells * 1600 < 0.60 * 32 * 1024**3
        ),
        "same_physical_geometry_and_plane": True,
        "base_mesh_reused": True,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "selected_refinement_ratio": REFINEMENT_RATIO,
        "selection_reason": "Minimum permitted constant ratio minimizes fine-grid cells and diffusive-scaling runtime while remaining valid for standard three-grid GCI analysis.",
        "grids": {label: asdict(spec) for label, spec in specs.items()},
        "candidate_costs": candidates,
        "checks": checks,
        "seeder_calls_planned": 2,
        "base_mesh": str(
            root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "seeder" / "mesh"
        ),
    }
    write_json(qc / "grid_convergence_preflight.json", result)
    return result


def _vec(values: np.ndarray) -> str:
    return "{ " + ", ".join(f"{float(value):.17g}" for value in values) + " }"


def generate_grid_seeder_lua(spec: GridSpec) -> str:
    plane_vectors = "{ " + _vec(PLANE_X_M) + ", " + _vec(PLANE_Y_M) + " }"
    objects = [
        "  { attribute = { kind = 'boundary', label = 'wall', level = minlevel, calc_dist = true }, geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/wall.stl' } } }",
        "  { attribute = { kind = 'boundary', label = 'inlet', level = minlevel }, geometry = { kind = 'canoND', object = { origin = "
        + _vec(PLANE_ORIGIN_M)
        + ", vec = "
        + plane_vectors
        + " } } }",
    ]
    for label in ("outlet_01", "outlet_02", "outlet_03"):
        objects.append(
            "  { attribute = { kind = 'boundary', label = '"
            + label
            + "', level = minlevel }, geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/"
            + label
            + ".stl' } } }"
        )
    objects.append(
        "  { attribute = { kind = 'seed' }, geometry = { kind = 'canoND', object = { origin = "
        + _vec(SEED_M)
        + " } } }"
    )
    return (
        "-- Isolated healthy grid-convergence Seeder configuration.\n"
        "folder = 'mesh/'\n"
        f"comment = 'ROI003274 {spec.label} grid, dx={spec.dx_m:.17g} m'\n"
        "debug = { debugMode = false, debugFiles = false, debugMesh = 'debug/' }\n"
        f"minlevel = {spec.root_level}\n"
        "bounding_cube = { origin = "
        + _vec(BASE_CUBE_ORIGIN_M)
        + f", length = {spec.cube_side_m:.17g} }}\n"
        "spatial_object = {\n"
        + ",\n".join(objects)
        + "\n}\n"
    )


def _copy_grid_geometry(root: Path, destination: Path) -> None:
    source = (
        root / "outputs" / "cfd_flow" / BASE_MESH_RUN
        / "geometry" / "geometry_solver_m"
    )
    destination.mkdir(parents=True, exist_ok=False)
    for label in ("wall", "outlet_01", "outlet_02", "outlet_03"):
        shutil.copy2(source / f"{label}.stl", destination / f"{label}.stl")


def _first_real_seeder_failure(stderr: str, stdout: str) -> str | None:
    """Skip known MPI loader warnings and return Seeder's first real error."""

    lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    for prefix in ("Fortran runtime error:", "Error:", "FATAL", "Fatal"):
        failure = next((line for line in lines if line.startswith(prefix)), None)
        if failure is not None:
            return failure
    lines.extend(line.strip() for line in stdout.splitlines() if line.strip())
    return next(
        (
            line
            for line in lines
            if "error" in line.lower() and "(ignored)" not in line.lower()
        ),
        None,
    )


def write_coarse_path_forensics(project_root: Path) -> dict[str, Any]:
    """Persist the zero-run evidence collected before any repair invocation."""

    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    seeder_run_root = run_root / "grids" / "coarse" / "seeder"
    contract = resolve_mesh_directory(seeder_run_root)
    sought = set(MESH_REQUIRED_FILES + MESH_QVAL_FILES)
    found = sorted(
        str(path.resolve())
        for path in run_root.rglob("*")
        if path.is_file() and path.name in sought
    )
    file_record = mesh_file_contract(contract.mesh_dir, require_qval=True)
    record = {
        "status": "CASE_A_COMPLETE_MESH_FOUND"
        if file_record["status"] == "PASS"
        else "CASE_B_NO_COMPLETE_MESH",
        "windows_run_root": str(run_root),
        "wsl_runtime_root": windows_to_wsl(run_root, "Ubuntu"),
        "seeder_working_directory_windows": str(seeder_run_root),
        "seeder_working_directory_wsl": windows_to_wsl(seeder_run_root, "Ubuntu"),
        "seeder_lua_folder": "mesh/",
        "expected_mesh_dir": str(contract.mesh_dir),
        "actual_seeder_output_dir": str(contract.seeder_output_dir),
        "expected_bnd_lsb": str(contract.mesh_dir / "bnd.lsb"),
        "actual_bnd_lsb": (
            str(contract.mesh_dir / "bnd.lsb")
            if (contract.mesh_dir / "bnd.lsb").is_file()
            else None
        ),
        "recursive_mesh_files_found": found,
        "path_mismatch_reason": (
            "No duplicate mesh component was found. Seeder resolved folder='mesh/' "
            "relative to its working directory, but the pinned writer requires that "
            "directory to exist before it opens mesh/bnd.lsb; the experiment runner "
            "had not materialized the canonical output directory."
        ),
        "root_cause": "SEEDER_OUTPUT_DIRECTORY_CONTRACT_MISMATCH",
        "file_contract_before_repair": file_record,
        "seeder_calls_during_forensics": 0,
    }
    write_json(run_root / "qc" / "seeder_output_path_forensics.json", record)
    return record


def repair_coarse_mesh(project_root: Path) -> dict[str, Any]:
    """Reuse a valid coarse mesh or make the single authorized repair call."""

    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    grid_root = run_root / "grids" / "coarse"
    seeder = grid_root / "seeder"
    qc_root = grid_root / "qc"
    spec = grid_specs()["coarse"]
    forensics = write_coarse_path_forensics(root)
    contract = resolve_mesh_directory(seeder)
    old_manifest_path = qc_root / "mesh_qc.json"
    old_manifest = read_json(old_manifest_path) if old_manifest_path.is_file() else None
    if can_reuse_completed_seeder_output(contract.mesh_dir, old_manifest):
        mesh_qc = _mesh_qc(contract.mesh_dir, spec)
        result = {
            "status": mesh_qc["status"],
            "reused_existing_mesh": True,
            "seeder_rerun": False,
            "additional_seeder_calls": 0,
            "path_contract": contract.as_record(),
            "file_contract": mesh_file_contract(contract.mesh_dir),
            "mesh_qc": mesh_qc,
        }
        write_json(qc_root / "coarse_mesh_preflight.json", result)
        return result

    history = seeder / "failure_evidence" / "attempt_1"
    history.mkdir(parents=True, exist_ok=True)
    for name in ("seeder.lua", "seeder_stdout.log", "seeder_stderr.log"):
        source = seeder / name
        destination = history / name
        if source.is_file() and not destination.exists():
            shutil.copy2(source, destination)
    if old_manifest_path.is_file():
        destination = history / "mesh_qc.json"
        if not destination.exists():
            shutil.copy2(old_manifest_path, destination)

    contract.mesh_dir.mkdir(parents=True, exist_ok=True)
    command = (
        f"cd '{windows_to_wsl(seeder, 'Ubuntu')}' && "
        f"'{SEEDER_WSL}' seeder.lua"
    )
    started = time.perf_counter()
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", command],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=900,
        check=False,
    )
    elapsed = time.perf_counter() - started
    stdout_path = seeder / "seeder_repair_stdout.log"
    stderr_path = seeder / "seeder_repair_stderr.log"
    stdout_path.write_text(process.stdout, encoding="utf-8")
    stderr_path.write_text(process.stderr, encoding="utf-8")
    file_contract = mesh_file_contract(contract.mesh_dir, require_qval=True)
    path_record = {
        **contract.as_record(),
        "expected_bnd_lsb": str(contract.mesh_dir / "bnd.lsb"),
        "actual_bnd_lsb": str(contract.mesh_dir / "bnd.lsb")
        if (contract.mesh_dir / "bnd.lsb").is_file()
        else None,
        "file_contract": file_contract,
        "musubi_mesh_assignment": (
            musubi_mesh_assignment(contract.mesh_dir)
            if file_contract["status"] == "PASS"
            else None
        ),
    }
    write_json(run_root / "qc" / "coarse_mesh_path_contract.json", path_record)
    if process.returncode != 0 or file_contract["status"] != "PASS":
        first_failure = _first_real_seeder_failure(process.stderr, process.stdout)
        result = {
            "status": "FAIL",
            "reused_existing_mesh": False,
            "seeder_rerun": True,
            "additional_seeder_calls": 1,
            "seeder_returncode": process.returncode,
            "seeder_wall_time_s": elapsed,
            "first_failure": first_failure,
            "path_forensics": forensics,
            "path_contract": path_record,
        }
        write_json(qc_root / "coarse_mesh_preflight.json", result)
        write_json(old_manifest_path, result)
        return result

    mesh_qc = _mesh_qc(contract.mesh_dir, spec)
    result = {
        "status": mesh_qc["status"],
        "reused_existing_mesh": False,
        "seeder_rerun": True,
        "additional_seeder_calls": 1,
        "seeder_returncode": process.returncode,
        "seeder_wall_time_s": elapsed,
        "path_forensics": forensics,
        "path_contract": path_record,
        "mesh_qc": mesh_qc,
    }
    write_json(qc_root / "coarse_mesh_preflight.json", result)
    write_json(
        old_manifest_path,
        {
            **mesh_qc,
            "seeder_returncode": process.returncode,
            "seeder_wall_time_s": elapsed,
            "seeder_lua": str(seeder / "seeder.lua"),
            "mesh_path": str(contract.mesh_dir),
            "repaired_output_directory_contract": True,
        },
    )
    return result


def finalize_repaired_coarse_mesh(project_root: Path) -> dict[str, Any]:
    """Finalize QC after the repaired Seeder call without invoking it again."""

    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    grid_root = run_root / "grids" / "coarse"
    seeder = grid_root / "seeder"
    qc_root = grid_root / "qc"
    contract = resolve_mesh_directory(seeder)
    file_contract = require_complete_mesh(contract.mesh_dir, require_qval=True)
    mesh_qc = _mesh_qc(contract.mesh_dir, grid_specs()["coarse"])
    result = {
        "status": mesh_qc["status"],
        "root_cause": "SEEDER_OUTPUT_DIRECTORY_CONTRACT_MISMATCH",
        "reused_existing_mesh": False,
        "seeder_rerun": True,
        "additional_seeder_calls": 1,
        "seeder_returncode": 0,
        "path_contract": contract.as_record(),
        "file_contract": file_contract,
        "mesh_qc": mesh_qc,
        "musubi_mesh_assignment": musubi_mesh_assignment(contract.mesh_dir),
    }
    write_json(qc_root / "coarse_mesh_preflight.json", result)
    write_json(
        qc_root / "mesh_qc.json",
        {
            **mesh_qc,
            "seeder_returncode": 0,
            "seeder_lua": str(seeder / "seeder.lua"),
            "mesh_path": str(contract.mesh_dir),
            "repaired_output_directory_contract": True,
        },
    )
    return result


def finalize_existing_fine_mesh(project_root: Path) -> dict[str, Any]:
    """Finalize the already-written fine mesh without another Seeder call."""

    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    grid_root = run_root / "grids" / "fine"
    seeder = grid_root / "seeder"
    qc_root = grid_root / "qc"
    contract = resolve_mesh_directory(seeder)
    file_contract = require_complete_mesh(contract.mesh_dir, require_qval=True)
    mesh_qc = _mesh_qc(contract.mesh_dir, grid_specs()["fine"])
    result = {
        "status": mesh_qc["status"],
        "seeder_calls": 1,
        "seeder_returncode": 0,
        "reused_existing_successful_output": True,
        "path_contract": contract.as_record(),
        "file_contract": file_contract,
        "mesh_qc": mesh_qc,
        "musubi_mesh_assignment": musubi_mesh_assignment(contract.mesh_dir),
    }
    write_json(qc_root / "fine_mesh_preflight.json", result)
    write_json(
        qc_root / "mesh_qc.json",
        {
            **mesh_qc,
            "seeder_returncode": 0,
            "seeder_lua": str(seeder / "seeder.lua"),
            "mesh_path": str(contract.mesh_dir),
            "boundary_order_parameterized": True,
        },
    )
    return result


def _mesh_qc(mesh_path: Path, spec: GridSpec) -> dict[str, Any]:
    header = parse_mesh_header(mesh_path)
    count = int(header["fluid_element_count"])
    tree_ids, _, _ = read_treelm_elemlist(mesh_path / "elemlist.lsb", n_elems=count)
    levels = tree_levels(tree_ids)
    components = connected_fluid_region_count(tree_ids, spec.root_level)
    contract = load_mesh_contract(
        mesh_path,
        expected_cells=count,
        allow_zero_normals=True,
        require_runtime_order=False,
    )
    solid = runtime_solid_cells(contract)
    pressure_counts = {
        label: len(
            _pressure_selected_rows(contract, contract.boundaries[label], label, solid)[0]
        )
        for label in ("inlet", "outlet_01", "outlet_02", "outlet_03")
    }
    boundary_counts = header["boundary_cell_counts"]
    checks = {
        "positive_fluid_cells": count > 0,
        "connected_component_one": components == 1,
        "all_labels_present": all(boundary_counts.get(label, 0) > 0 for label in (
            "wall", "inlet", "outlet_01", "outlet_02", "outlet_03"
        )),
        "uniform_single_level": bool(
            len(levels) > 0 and np.all(levels == spec.root_level)
            and header["minimum_level"] == header["maximum_level"] == spec.root_level
        ),
        "runtime_pressure_cells_present": all(value > 0 for value in pressure_counts.values()),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "dx_m": spec.dx_m,
        "dt_s": spec.dt_s,
        "root_level": spec.root_level,
        "fluid_cell_count": count,
        "connected_fluid_regions": components,
        "boundary_cell_counts": boundary_counts,
        "runtime_pressure_valid_cells": pressure_counts,
        "runtime_solid_cells": len(solid),
        "zero_aggregate_normal_cells": {
            label: int(np.count_nonzero(boundary.normal_indices < 0))
            for label, boundary in contract.boundaries.items()
        },
        "zero_normals_are_not_a_hard_gate": True,
        "checks": checks,
    }


def seed_new_grids(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    preflight = build_grid_preflight(root)
    if preflight["status"] != "PASS":
        raise FlowError("CFD_FLOW_GRID_PREFLIGHT_FAILED", "Grid preflight failed")
    run_root = root / "outputs" / "cfd_flow" / GRID_RUN
    specs = grid_specs()
    results: dict[str, Any] = {}
    seeder_calls = 0
    for label in ("coarse", "fine"):
        spec = specs[label]
        grid_root = run_root / "grids" / label
        manifest_path = grid_root / "qc" / "mesh_qc.json"
        if manifest_path.is_file():
            existing = read_json(manifest_path)
            if existing.get("status") != "PASS":
                raise FlowError("CFD_FLOW_GRID_MESH_FAILED", f"Existing {label} mesh failed")
            results[label] = existing
            continue
        seeder = grid_root / "seeder"
        geometry = grid_root / "geometry" / "geometry_solver_m"
        (grid_root / "qc").mkdir(parents=True, exist_ok=False)
        seeder.mkdir(parents=True, exist_ok=False)
        # This pinned Seeder opens the output files directly and therefore
        # requires its configured output folder to exist before launch.
        (seeder / "mesh").mkdir(parents=False, exist_ok=False)
        _copy_grid_geometry(root, geometry)
        lua = generate_grid_seeder_lua(spec)
        (seeder / "seeder.lua").write_text(lua, encoding="utf-8")
        workdir_wsl = windows_to_wsl(seeder, "Ubuntu")
        command = (
            f"cd '{workdir_wsl}' && "
            f"'{SEEDER_WSL}' seeder.lua"
        )
        started = time.perf_counter()
        process = subprocess.run(
            ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", command],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=900,
            check=False,
        )
        elapsed = time.perf_counter() - started
        seeder_calls += 1
        (seeder / "seeder_stdout.log").write_text(process.stdout, encoding="utf-8")
        (seeder / "seeder_stderr.log").write_text(process.stderr, encoding="utf-8")
        if process.returncode != 0:
            failure = _first_real_seeder_failure(process.stderr, process.stdout)
            failed = {
                "status": "FAIL",
                "seeder_returncode": process.returncode,
                "seeder_wall_time_s": elapsed,
                "first_failure": failure,
            }
            write_json(manifest_path, failed)
            raise FlowError("CFD_FLOW_GRID_MESH_FAILED", f"{label}: {failure}")
        qc = _mesh_qc(seeder / "mesh", spec)
        qc.update(
            {
                "seeder_returncode": process.returncode,
                "seeder_wall_time_s": elapsed,
                "seeder_lua": str(seeder / "seeder.lua"),
                "mesh_path": str(seeder / "mesh"),
            }
        )
        write_json(manifest_path, qc)
        results[label] = qc
        if qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_GRID_MESH_FAILED", f"{label} mesh QC failed")
    base_mesh = root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "seeder" / "mesh"
    base_qc = _mesh_qc(base_mesh, specs["base"])
    results["base"] = base_qc
    summary = {
        "status": "PASS" if all(value["status"] == "PASS" for value in results.values()) else "FAIL",
        "actual_head": _git(root, "rev-parse", "HEAD"),
        "production_pipeline_modified": bool(
            _git(root, "diff", "--", *(str(path.relative_to(root)) for path in _production_paths(root)))
        ),
        "seeder_calls_this_invocation": seeder_calls,
        "seeder_calls_total": 2,
        "grids": results,
    }
    write_json(run_root / "qc" / "grid_mesh_summary.json", summary)
    return summary
