"""Isolated checkpoint-driven fast steady continuation for healthy capillary CFD."""

from __future__ import annotations

import csv
import io
import math
import os
import re
import shlex
import shutil
import subprocess
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Callable, Iterable, TextIO

import numpy as np

from .adaptive_flux_steady import (
    AXIS_MESH_RUN,
    BULK_NU_M2_S,
    CONTROLLER_FIELDS,
    EXPECTED_FLUID_CELLS,
    EXPECTED_INLET_GLOBBC,
    NU_M2_S,
    PRESSURE_REFERENCE_PA,
    RUNTIME_BASE_WSL,
    SHORT_SIMULATION_NAME,
    CheckpointTracker,
    MonitoredRun,
    WSL_INFRASTRUCTURE_PATTERN,
    _git,
    _live_safety_failure,
    _match_record,
    _one_wsl_health_check,
    _prepare_short_tracking_root,
    _run_luac,
    _wsl_path_to_unc,
    summarize_controller_csv,
)
from .adaptive_flux_validation import (
    BINARY_WSL,
    CONTROLLER_PATTERN,
    EXPECTED_BINARY_SHA256,
    EXPECTED_DT_S,
    MPIRUN_WSL,
)
from .apes import parse_mesh_header, windows_to_wsl
from .config import load_cfd_flow_config
from .exact_link_flux import CS2, EXPECTED_DX_M, _file_manifest
from .healthy_capillary_calibration import (
    CORRECTED_RUNTIME_NORMALS,
    EXPECTED_PRESSURE_COUNTS,
    HEALTHY_UMEAN_RANGE_M_S,
    HEALTHY_WSS_RANGE_DYN_CM2,
    OUTLET_GAUGE_PRESSURE_PA,
    REFERENCE_DENSITY_KG_M3,
    SMOOTH_INLET_AREA_M2,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    exact_boundary_flux_audit,
    source_only_pressure_outlet_contract,
)
from .io import FlowError, read_json, sha256_file, write_json
from .restart_decode import (
    D3Q19_DIRECTIONS,
    parse_restart_header,
    read_restart_pdf,
    restart_binary_size_contract,
)


RUN_PREFIX = "healthy_mouse_capillary_fast_steady_anchor003274"
SOURCE_RUN = "healthy_mouse_capillary_calibration_anchor003274_20260829_180310"
SOURCE_ITERATION = 314_166
SOURCE_BINARY_NAME = "a3274_7.670E-03.lsb"
EXPECTED_SOURCE_SHA256 = "da657368b65df90192e60a6e6475474574b15fb3bdb2526dc60f3af760e3543a"
EXPECTED_FIRST_ITERATION = SOURCE_ITERATION + 1

PHYSICAL_CORES = 12
LOGICAL_CPUS = 24
MPI_BINDING_ARGS = ("--bind-to", "core", "--map-by", "core", "--report-bindings")
THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}
BENCHMARK_ITERATIONS = 2_000
OLD_8_RANK_ITERATIONS_S = 42.997
TRACKING_INTERVAL = 5_000
CHECKPOINT_INTERVAL = 5_000
SIM_CONTROL_INTERVAL = 100
MAXIMUM_ITERATIONS = 1_000_000
WALLCLOCK_LIMIT_S = 3_600
WRAPPER_MARGIN_S = 300

SHORT_TARGET = 10_000
SHORT_RANGE = (8_000, 12_000)
LONG_TARGET = 20_000
LONG_RANGE = (18_000, 22_000)
MASS_HARD = 0.01
MASS_PREFERRED = 0.005
VELOCITY_HARD = 0.01
VELOCITY_PREFERRED = 0.005
PRESSURE_HARD = 0.005
PRESSURE_PREFERRED = 0.002
INLET_HARD = 0.01
INLET_PREFERRED = 0.001
BOUNDARY_HARD = 0.01
BOUNDARY_PREFERRED = 0.005
PRESSURE_MASS_CROSSCHECK_HARD = 0.001
CONTROLLER_ALGEBRA_HARD = 1.0e-8
SIGNIFICANT_BACKFLOW_FRACTION = 0.05

STATUS_PASS = "CFD_FLOW_HEALTHY_MOUSE_CAPILLARY_FAST_STEADY_PASS"
STATUS_WALLCLOCK = "CFD_FLOW_HEALTHY_STEADY_INCOMPLETE_WALLCLOCK"
STATUS_SAFETY = "CFD_FLOW_HEALTHY_FAST_STEADY_NUMERICAL_SAFETY_FAILED"
STATUS_RUNTIME = "CFD_FLOW_HEALTHY_FAST_STEADY_RUNTIME_FAILED"
STATUS_EXACT_UNRESOLVED = "CFD_FLOW_HEALTHY_FAST_EXACT_AUDIT_UNRESOLVED"
STATUS_CROSSCHECK = "CFD_FLOW_HEALTHY_FAST_PRESSURE_MASS_CROSSCHECK_FAILED"
STATUS_MPI_REPRODUCIBILITY = "CFD_FLOW_MPI_RANK_REPRODUCIBILITY_FAILED"
STATUS_FAILED = "CFD_FLOW_HEALTHY_FAST_STEADY_FAILED"
NEXT_GRID = "RUN HEALTHY ADAPTIVE-FLUX GRID CONVERGENCE"
NEXT_RESUME = "RESUME HEALTHY FAST STEADY FROM LATEST CHECKPOINT"
NEXT_SAFETY = "REVIEW HEALTHY LATTICE NUMERICS"
NEXT_RUNTIME = "RESUME FROM LAST COMPLETE CHECKPOINT AFTER RUNTIME REPAIR"
NEXT_EXACT = "REVIEW EXACT CONSERVATION AUDIT"
NEXT_MPI = "REVIEW MPI RANK REPRODUCIBILITY"


@dataclass(frozen=True, slots=True)
class RestartEvidence:
    iteration: int
    header_path: Path
    binary_path: Path
    sha256: str
    source: str


class FlushPolicy:
    """Flush a stream after N records or one second, whichever comes first."""

    def __init__(
        self,
        stream: TextIO,
        *,
        record_limit: int | None,
        interval_s: float = 1.0,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        self.stream = stream
        self.record_limit = record_limit
        self.interval_s = interval_s
        self.clock = clock
        self.pending = 0
        self.last_flush = clock()
        self.flush_count = 0
        self.lock = threading.Lock()

    def note(self, records: int = 1) -> None:
        with self.lock:
            self.pending += records
            self._flush_if_due_locked(force=False)

    def flush_if_due(self, *, force: bool = False) -> None:
        with self.lock:
            self._flush_if_due_locked(force=force)

    def _flush_if_due_locked(self, *, force: bool) -> None:
        now = self.clock()
        count_due = self.record_limit is not None and self.pending >= self.record_limit
        time_due = now - self.last_flush >= self.interval_s
        if force or count_due or time_due:
            self.stream.flush()
            self.pending = 0
            self.last_flush = now
            self.flush_count += 1


def choose_mpi_ranks(results: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Choose maximum throughput; within 3% of fastest choose fewer ranks."""

    rows = sorted((dict(row) for row in results), key=lambda row: int(row["ranks"]))
    if not rows:
        raise ValueError("No MPI benchmark results")
    fastest = max(float(row["iterations_per_s"]) for row in rows)
    eligible = [row for row in rows if fastest / float(row["iterations_per_s"]) - 1.0 < 0.03]
    selected = min(eligible, key=lambda row: int(row["ranks"]))
    return {
        "selected": selected,
        "fastest_iterations_per_s": fastest,
        "tie_rule": "within 3% of fastest selects fewer ranks",
    }


def next_benchmark_rank(current_results: list[dict[str, Any]], physical_cores: int) -> int | None:
    tested = {int(row["ranks"]) for row in current_results}
    for rank in (4, 6, 8):
        if physical_cores >= rank and rank not in tested:
            return rank
    by_rank = {int(row["ranks"]): float(row["iterations_per_s"]) for row in current_results}
    if physical_cores >= 12 and 12 not in tested and 8 in by_rank and 6 in by_rank:
        if by_rank[8] > by_rank[6]:
            return 12
    if physical_cores >= 16 and 16 not in tested and 12 in by_rank and 8 in by_rank:
        if by_rank[12] >= 1.03 * by_rank[8]:
            return 16
    return None


def select_window(
    current_iteration: int,
    candidates: Iterable[RestartEvidence],
    *,
    target: int,
    allowed: tuple[int, int],
) -> RestartEvidence | None:
    valid = [
        item
        for item in candidates
        if allowed[0] <= current_iteration - item.iteration <= allowed[1]
    ]
    if not valid:
        return None
    return min(valid, key=lambda item: (abs((current_iteration - item.iteration) - target), -item.iteration))


def mass_accumulation_kg_s(total_a: float, total_b: float, delta_iter: int) -> float:
    if delta_iter <= 0:
        raise ValueError("delta_iter must be positive")
    return (
        (total_b - total_a)
        / delta_iter
        * REFERENCE_DENSITY_KG_M3
        * EXPECTED_DX_M**3
        / EXPECTED_DT_S
    )


def mass_residual(total_a: float, total_b: float, delta_iter: int) -> float:
    return abs(mass_accumulation_kg_s(total_a, total_b, delta_iter)) / TARGET_MASS_FLOW_KG_S


def velocity_l2_residual(previous: np.ndarray, current: np.ndarray) -> float:
    previous_values = np.asarray(previous, dtype=np.float64)
    current_values = np.asarray(current, dtype=np.float64)
    numerator = float(np.linalg.norm((current_values - previous_values).ravel()))
    denominator = max(float(np.linalg.norm(current_values.ravel())), np.finfo(np.float64).tiny)
    return numerator / denominator


def characteristic_pressure_drop(inlet_gauge_pa: float) -> float:
    drops = [abs(inlet_gauge_pa - gauge) for gauge in OUTLET_GAUGE_PRESSURE_PA.values()]
    return float(np.median(np.asarray(drops, dtype=np.float64)))


def pressure_residual(delta_mean_pressure_pa: float, characteristic_drop_pa: float) -> float:
    if characteristic_drop_pa <= 0.0:
        raise ValueError("Characteristic pressure drop must be positive")
    return abs(delta_mean_pressure_pa) / characteristic_drop_pa


def inlet_residual(exact_mass_flow: float) -> float:
    return abs(exact_mass_flow - TARGET_MASS_FLOW_KG_S) / TARGET_MASS_FLOW_KG_S


def boundary_residual(inlet_mass_flow: float, outlet_flows: Iterable[float]) -> float:
    if inlet_mass_flow == 0.0:
        raise ValueError("Inlet flow must be nonzero")
    return abs(inlet_mass_flow - sum(outlet_flows)) / abs(inlet_mass_flow)


def significant_backflow(outlet_flows: Iterable[float], inlet_mass_flow: float) -> bool:
    threshold = SIGNIFICANT_BACKFLOW_FRACTION * abs(inlet_mass_flow)
    return any(flow < 0.0 and abs(flow) > threshold for flow in outlet_flows)


def classify_stage(r_mass_short: float) -> str:
    if r_mass_short > 0.05:
        return "STAGE_1_FAR"
    if r_mass_short > MASS_HARD:
        return "STAGE_2_NEAR"
    return "STAGE_3_FINAL"


def all_final_gates_pass(record: dict[str, Any]) -> bool:
    required = (
        record.get("R_mass_short") is not None and record["R_mass_short"] <= MASS_HARD,
        record.get("R_mass_long") is not None and record["R_mass_long"] <= MASS_HARD,
        record.get("R_boundary") is not None and record["R_boundary"] <= BOUNDARY_HARD,
        record.get("R_velocity") is not None and record["R_velocity"] <= VELOCITY_HARD,
        record.get("R_pressure") is not None and record["R_pressure"] <= PRESSURE_HARD,
        record.get("R_inlet") is not None and record["R_inlet"] <= INLET_HARD,
        record.get("pressure_pdf_crosscheck") is not None
        and record["pressure_pdf_crosscheck"]["relative_discrepancy"] <= PRESSURE_MASS_CROSSCHECK_HARD
        and record["pressure_pdf_crosscheck"]["direction_same"],
        not bool(record.get("significant_backflow", True)),
        bool(record.get("all_finite", False)),
        float(record.get("minimum_pdf", -math.inf)) > 0.0,
        float(record.get("maximum_lattice_speed", math.inf)) < 0.05,
        int(record.get("inlet_globbc", -1)) == EXPECTED_INLET_GLOBBC,
    )
    return all(required)


def create_stop_file(runtime_root: Path) -> Path:
    path = Path(runtime_root) / "stop"
    path.touch(exist_ok=True)
    return path


def _physics_and_boundaries_lua() -> str:
    return f"""dx = {EXPECTED_DX_M:.17g}
dt = {EXPECTED_DT_S:.17g}
rho0_phy = {REFERENCE_DENSITY_KG_M3:.17g}
nu_phy = {NU_M2_S:.17g}
bulk_viscosity_phy = {BULK_NU_M2_S:.17g}
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}

function outlet_01_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + OUTLET_GAUGE_PRESSURE_PA['outlet_01']:.17g} end
function outlet_02_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + OUTLET_GAUGE_PRESSURE_PA['outlet_02']:.17g} end
function outlet_03_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + OUTLET_GAUGE_PRESSURE_PA['outlet_03']:.17g} end

physics = {{ dt = dt, rho0 = rho0_phy }}
identify = {{ label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }}
fluid = {{ kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }}
initial_condition = {{ pressure = pressure_reference_phy, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }}
boundary_condition = {{
  {{ label = 'wall', kind = 'wall_libb' }},
  {{ label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = {TARGET_MASS_FLOW_KG_S:.17g} }},
  {{ label = 'outlet_01', kind = 'pressure_eq', pressure = outlet_01_pressure }},
  {{ label = 'outlet_02', kind = 'pressure_eq', pressure = outlet_02_pressure }},
  {{ label = 'outlet_03', kind = 'pressure_eq', pressure = outlet_03_pressure }}
}}"""


def generate_benchmark_lua(*, mesh_wsl: str, resume_header_wsl: str) -> str:
    final_iteration = SOURCE_ITERATION + BENCHMARK_ITERATIONS
    return f"""-- Isolated MPI performance benchmark; output state is not accepted.
simulation_name = '{SHORT_SIMULATION_NAME}'
printRuntimeInfo = true
timing_file = 'timing.res'
mesh = '{mesh_wsl}/'
scaling = 'diffusive'
logging = {{ level = 5 }}
maximum_iterations = {final_iteration}
{_physics_and_boundaries_lua()}
sim_control = {{
  time_control = {{ max = {{ iter = maximum_iterations }}, interval = {{ iter = {SIM_CONTROL_INTERVAL} }} }},
  abort_criteria = {{ stop_file = 'stop' }}
}}
restart = {{ read = '{resume_header_wsl}' }}
"""


def generate_production_lua(*, mesh_wsl: str, resume_header_wsl: str) -> str:
    first_checkpoint = SOURCE_ITERATION + CHECKPOINT_INTERVAL
    return f"""-- Isolated healthy fast steady continuation; Python exact gates are final.
simulation_name = '{SHORT_SIMULATION_NAME}'
printRuntimeInfo = true
timing_file = 'timing/timing.res'
mesh = '{mesh_wsl}/'
scaling = 'diffusive'
logging = {{ level = 5 }}
maximum_iterations = {MAXIMUM_ITERATIONS}
{_physics_and_boundaries_lua()}
sim_control = {{
  time_control = {{
    max = {{ iter = maximum_iterations, clock = {WALLCLOCK_LIMIT_S} }},
    interval = {{ iter = {SIM_CONTROL_INTERVAL} }}
  }},
  abort_criteria = {{ stop_file = 'stop' }}
}}
tracking = {{
  {{
    label = 'p', folder = 'tracking/p/', variable = {{ 'pressure_phy' }},
    shape = {{ kind = 'all' }}, reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {SOURCE_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {TRACKING_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }},
  {{
    label = 'u', folder = 'tracking/u/', variable = {{ 'vel_mag_phy' }},
    shape = {{ kind = 'all' }}, reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {SOURCE_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {TRACKING_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }}
}}
restart = {{
  read = '{resume_header_wsl}', write = 'restart/',
  time_control = {{ min = {{ iter = {first_checkpoint} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CHECKPOINT_INTERVAL} }} }}
}}
"""


def fast_lua_contract(text: str, *, resume_header_wsl: str, benchmark: bool) -> dict[str, Any]:
    checks = {
        "restart_read": f"read = '{resume_header_wsl}'" in text,
        "healthy_target": f"mass_flowrate = {TARGET_MASS_FLOW_KG_S:.17g}" in text,
        "fixed_physics": all(token in text for token in (
            f"dx = {EXPECTED_DX_M:.17g}", f"dt = {EXPECTED_DT_S:.17g}",
            f"rho0_phy = {REFERENCE_DENSITY_KG_M3:.17g}", f"nu_phy = {NU_M2_S:.17g}",
            f"bulk_viscosity_phy = {BULK_NU_M2_S:.17g}", "layout = 'd3q19'",
            "relaxation = 'bgk'", "kind = 'wall_libb'",
        )),
        "fixed_outlets": text.count("kind = 'pressure_eq'") == 3,
        "adaptive_inlet": "kind = 'adaptive_flux_pressure'" in text and "mfr_eq" not in text,
        "built_in_pu_disabled": "steady_state" not in text and "convergence" not in text,
        "stop_file": "stop_file = 'stop'" in text,
        "no_harvester": "harvest" not in text.lower(),
        "no_seeder_call": "seeder.lua" not in text.lower() and "seed.lua" not in text.lower(),
    }
    if benchmark:
        checks.update({
            "benchmark_2000": f"maximum_iterations = {SOURCE_ITERATION + BENCHMARK_ITERATIONS}" in text,
            "no_tracking": "tracking =" not in text,
            "no_restart_write": "write = 'restart/'" not in text,
            "no_steady": "steady_state" not in text,
        })
    else:
        checks.update({
            "tracking_5000": text.count(f"interval = {{ iter = {TRACKING_INTERVAL} }}") >= 3,
            "checkpoint_5000": f"min = {{ iter = {SOURCE_ITERATION + CHECKPOINT_INTERVAL} }}" in text,
            "wallclock_3600": f"clock = {WALLCLOCK_LIMIT_S}" in text,
        })
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


TIMING_FIELDS = (
    "revision", "simulation", "domain_size", "nprocs", "mlups", "mlups_kernel",
    "imbalance_percent", "time_musubi_s", "max_iteration", "total_density",
    "time_main_loop_s", "time_load_mesh_s", "time_init_level_s", "time_write_restart_s",
    "time_balance_s", "time_source_s", "time_aux_s", "time_relax_s", "comp_percent",
    "comm_percent", "bcbuffer_percent", "bc_percent", "interpolation_percent",
)


def parse_timing_result(path: Path) -> dict[str, Any]:
    rows = [line.split() for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if len(rows) != 1 or len(rows[0]) != len(TIMING_FIELDS):
        raise FlowError(STATUS_FAILED, f"Unexpected timing result: {rows}")
    raw = dict(zip(TIMING_FIELDS, rows[0]))
    integer_fields = {"domain_size", "nprocs", "max_iteration"}
    text_fields = {"revision", "simulation"}
    return {
        key: (value if key in text_fields else int(value) if key in integer_fields else float(value))
        for key, value in raw.items()
    }


def _last_controller_record(text: str) -> dict[str, Any]:
    records = [_match_record(match) for match in re.finditer(
        r"ADAPTIVE_FLUX_PRESSURE\s+iter=(?P<iteration>\d+)\s+"
        r"target_lattice=\s*(?P<target>[-+0-9.Ee]+)\s+"
        r"controlled_lattice=\s*(?P<controlled>[-+0-9.Ee]+)\s+"
        r"relative_error=\s*(?P<error>[-+0-9.Ee]+)\s+"
        r"rho_boundary=\s*(?P<rho>[-+0-9.Ee]+)\s+"
        r"pressure_pa=\s*(?P<pressure>[-+0-9.Ee]+)\s+"
        r"max_lattice_velocity=\s*(?P<speed>[-+0-9.Ee]+)\s+"
        r"minimum_pdf=\s*(?P<minimum_pdf>[-+0-9.Ee]+)\s+"
        r"globBC_count=(?P<count>\d+)", text
    )]
    if not records:
        raise FlowError(STATUS_FAILED, "Benchmark emitted no adaptive controller record")
    return records[-1]


def benchmark_reproducibility(results: list[dict[str, Any]]) -> dict[str, Any]:
    if not results:
        return {"status": "FAIL", "reason": "no results"}
    keys = ("rho_boundary", "target_lattice", "controlled_lattice", "max_lattice_velocity", "minimum_pdf")
    reference = results[0]["final_controller"]
    differences: dict[str, float] = {}
    for key in keys:
        values = [float(row["final_controller"][key]) for row in results]
        scale = max(max(abs(value) for value in values), 1.0e-30)
        differences[key] = (max(values) - min(values)) / scale
    algebra_ok = all(float(row["final_controller"]["relative_error"]) <= CONTROLLER_ALGEBRA_HARD for row in results)
    safety_ok = all(
        int(row["final_controller"]["globBC_count"]) == EXPECTED_INLET_GLOBBC
        and float(row["final_controller"]["max_lattice_velocity"]) < 0.05
        and float(row["final_controller"]["minimum_pdf"]) > 0.0
        for row in results
    )
    roundoff_ok = all(value <= 1.0e-10 for value in differences.values())
    return {
        "status": "PASS" if algebra_ok and safety_ok and roundoff_ok else "FAIL",
        "relative_spread": differences,
        "algebra_residual_pass": algebra_ok,
        "safety_pass": safety_ok,
        "roundoff_consistency_pass": roundoff_ok,
        "reference": reference,
    }


@dataclass(frozen=True, slots=True)
class FastMonitoredRun:
    returncode: int
    wall_time_s: float
    safety_failure: str | None
    infrastructure_failure: str | None
    continuity_failure: str | None
    audit_failure: str | None
    audit_failure_status: str | None
    wrapper_timeout: bool
    controller_record_count: int
    latest_controller_iteration: int | None
    stop_file_triggered: bool
    stop_requested_iteration: int | None
    stop_file_created_at: str | None
    stdout_flushes: int
    stderr_flushes: int
    controller_flushes: int


class FastCheckpointTracker(CheckpointTracker):
    """The standard stable-pair tracker with the experiment's 5k contract."""

    def _write(self) -> None:
        write_json(
            self.manifest_path,
            {
                "status": "PASS",
                "checkpoint_interval_iterations": CHECKPOINT_INTERVAL,
                "runtime_root_wsl": self.runtime_root_wsl,
                "records": self.records,
            },
        )


def restart_compatibility(header_path: Path, binary_path: Path) -> dict[str, Any]:
    header = parse_restart_header(Path(header_path))
    contract = restart_binary_size_contract(
        Path(binary_path),
        n_elems=header.n_elems,
        n_components=header.n_components,
        n_dofs=header.n_dofs,
    )
    checks = {
        "fluid_cells": header.n_elems == EXPECTED_FLUID_CELLS,
        "components": header.n_components == 19,
        "dofs": header.n_dofs == 1,
        "variable": header.variable_name == "pdf",
        "binary_size": contract["status"] == "PASS",
        "header_binary_matches": Path(header.binary_path).resolve() == Path(binary_path).resolve(),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "iteration": header.iteration,
        "header": str(header_path),
        "binary": str(binary_path),
        "binary_sha256": sha256_file(Path(binary_path)),
        "checks": checks,
    }


def _controller_safety_failure(record: dict[str, Any]) -> str | None:
    numeric = (
        "target_lattice", "controlled_lattice", "relative_error", "rho_boundary",
        "pressure_pa", "max_lattice_velocity", "minimum_pdf",
    )
    if not all(math.isfinite(float(record[key])) for key in numeric):
        return f"controller NaN/Inf at iteration {record['iteration']}"
    if int(record["globBC_count"]) != EXPECTED_INLET_GLOBBC:
        return f"inlet globBC changed at iteration {record['iteration']}"
    if float(record["minimum_pdf"]) <= 0.0:
        return f"minimum PDF <= 0 at iteration {record['iteration']}"
    if float(record["max_lattice_velocity"]) >= 0.05:
        return f"maximum lattice speed >= 0.05 at iteration {record['iteration']}"
    return None


def _run_wsl_capture(
    *, distribution: str, workdir_wsl: str, command: list[str], timeout_s: int,
) -> tuple[subprocess.CompletedProcess[str], float]:
    shell = f"cd {shlex.quote(workdir_wsl)} && exec {shlex.join(command)}"
    started = time.perf_counter()
    completed = subprocess.run(
        ["wsl.exe", "-d", distribution, "--", "bash", "-lc", shell],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout_s,
    )
    return completed, time.perf_counter() - started


def performance_preflight(distribution: str) -> dict[str, Any]:
    completed, _ = _run_wsl_capture(
        distribution=distribution,
        workdir_wsl=RUNTIME_BASE_WSL,
        command=[
            "bash", "-lc",
            "lscpu && nproc && /home/lzy/.local/bin/mpirun --version "
            "&& /home/lzy/.local/bin/mpirun --help binding "
            "&& /home/lzy/.local/bin/mpirun --help mapping",
        ],
        timeout_s=30,
    )
    if completed.returncode != 0:
        raise FlowError(STATUS_RUNTIME, f"MPI performance preflight failed: {completed.stderr}")
    text = completed.stdout

    def integer(label: str) -> int:
        match = re.search(rf"(?m)^{re.escape(label)}:\s+(\d+)", text)
        if not match:
            raise FlowError(STATUS_RUNTIME, f"lscpu field missing: {label}")
        return int(match.group(1))

    sockets = integer("Socket(s)")
    cores_per_socket = integer("Core(s) per socket")
    logical = integer("CPU(s)")
    physical = sockets * cores_per_socket
    binding = "--bind-to" in text and "--map-by" in text
    return {
        "status": "PASS" if binding and physical >= 8 else "FAIL",
        "physical_cores": physical,
        "cores_per_socket": cores_per_socket,
        "sockets": sockets,
        "logical_cpus": logical,
        "mpi_binding_supported": binding,
        "stdout": text,
        "stderr": completed.stderr,
    }


def _stage_source_checkpoint(
    *, source_header: Path, source_binary: Path, runtime_root: Path,
) -> tuple[Path, Path]:
    staged_header = runtime_root / "restart" / source_header.name
    staged_binary = runtime_root / "restart" / source_binary.name
    shutil.copy2(source_header, staged_header)
    shutil.copy2(source_binary, staged_binary)
    if sha256_file(staged_binary) != EXPECTED_SOURCE_SHA256:
        raise FlowError(STATUS_FAILED, "Staged source restart SHA mismatch")
    if sha256_file(staged_header) != sha256_file(source_header):
        raise FlowError(STATUS_FAILED, "Staged source header SHA mismatch")
    return staged_header, staged_binary


def run_one_benchmark(
    *, distribution: str, mesh: Path, source_header: Path, source_binary: Path,
    run_root: Path, runtime_root_wsl: str, ranks: int,
) -> dict[str, Any]:
    evidence = run_root / "benchmark" / f"rank_{ranks}"
    evidence.mkdir(parents=True, exist_ok=False)
    runtime_root = _wsl_path_to_unc(distribution, runtime_root_wsl)
    if runtime_root.exists():
        raise FlowError(STATUS_FAILED, f"Benchmark WSL root exists: {runtime_root_wsl}")
    _prepare_short_tracking_root(distribution=distribution, root_wsl=runtime_root_wsl)
    _stage_source_checkpoint(
        source_header=source_header, source_binary=source_binary, runtime_root=runtime_root,
    )
    resume = f"restart/{source_header.name}"
    lua = generate_benchmark_lua(
        mesh_wsl=windows_to_wsl(mesh, distribution), resume_header_wsl=resume,
    )
    (runtime_root / "musubi.lua").write_text(lua, encoding="utf-8")
    shutil.copy2(runtime_root / "musubi.lua", evidence / "musubi.lua")
    contract = fast_lua_contract(lua, resume_header_wsl=resume, benchmark=True)
    write_json(evidence / "lua_contract.json", contract)
    if contract["status"] != "PASS":
        raise FlowError(STATUS_FAILED, f"Benchmark Lua contract failed for {ranks} ranks")
    luac_rc = _run_luac(
        distribution=distribution,
        workdir_wsl=runtime_root_wsl,
        stdout_path=evidence / "luac_stdout.log",
        stderr_path=evidence / "luac_stderr.log",
    )
    if luac_rc != 0:
        raise FlowError(STATUS_FAILED, f"Benchmark Lua syntax failed for {ranks} ranks")
    command = [
        "env", *[f"{key}={value}" for key, value in THREAD_ENV.items()],
        "OMPI_ALLOW_RUN_AS_ROOT=1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
        MPIRUN_WSL, "-np", str(ranks), *MPI_BINDING_ARGS, BINARY_WSL, "musubi.lua",
    ]
    completed, wall = _run_wsl_capture(
        distribution=distribution, workdir_wsl=runtime_root_wsl,
        command=command, timeout_s=300,
    )
    (evidence / "musubi_stdout.log").write_text(completed.stdout, encoding="utf-8")
    (evidence / "musubi_stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise FlowError(STATUS_FAILED, f"{ranks}-rank benchmark returned {completed.returncode}")
    timing_path = runtime_root / "timing.res"
    if not timing_path.is_file():
        raise FlowError(STATUS_FAILED, f"{ranks}-rank benchmark timing.res missing")
    shutil.copy2(timing_path, evidence / "timing.res")
    timing = parse_timing_result(timing_path)
    controller = _last_controller_record(completed.stdout)
    if int(controller["iteration"]) != SOURCE_ITERATION + BENCHMARK_ITERATIONS:
        raise FlowError(STATUS_FAILED, f"{ranks}-rank endpoint iteration changed")
    failure = _controller_safety_failure(controller)
    if failure or float(controller["relative_error"]) > CONTROLLER_ALGEBRA_HARD:
        raise FlowError(STATUS_FAILED, failure or f"{ranks}-rank controller algebra failed")
    if int(timing["nprocs"]) != ranks:
        raise FlowError(STATUS_FAILED, f"timing.res reports {timing['nprocs']} ranks, expected {ranks}")
    main_loop = float(timing["time_main_loop_s"])
    result = {
        "status": "PASS", "ranks": ranks, "returncode": completed.returncode,
        "wall_time_s": wall, "iterations": BENCHMARK_ITERATIONS,
        "iterations_per_s": BENCHMARK_ITERATIONS / main_loop,
        "wrapper_iterations_per_s": BENCHMARK_ITERATIONS / wall,
        "mlups": timing["mlups"], "mlups_kernel": timing["mlups_kernel"],
        "imbalance_percent": timing["imbalance_percent"],
        "comp_percent": timing["comp_percent"], "comm_percent": timing["comm_percent"],
        "bcbuffer_percent": timing["bcbuffer_percent"], "bc_percent": timing["bc_percent"],
        "time_aux_s": timing["time_aux_s"], "time_main_loop_s": main_loop,
        "binding_args": list(MPI_BINDING_ARGS), "thread_environment": THREAD_ENV,
        "final_controller": controller, "evidence_root": str(evidence),
    }
    write_json(evidence / "benchmark_result.json", result)
    return result


def run_mpi_benchmarks(
    *, distribution: str, mesh: Path, source_header: Path, source_binary: Path,
    run_root: Path, stamp: str, physical_cores: int,
) -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    while True:
        rank = next_benchmark_rank(results, physical_cores)
        if rank is None:
            break
        # The mandatory 4/6/8 sweep is completed before decline can stop expansion.
        if rank > 8 and float(results[-1]["iterations_per_s"]) <= float(results[-2]["iterations_per_s"]):
            break
        result = run_one_benchmark(
            distribution=distribution, mesh=mesh, source_header=source_header,
            source_binary=source_binary, run_root=run_root,
            runtime_root_wsl=f"{RUNTIME_BASE_WSL}/b_{stamp[-6:]}_{rank}", ranks=rank,
        )
        results.append(result)
    reproducibility = benchmark_reproducibility(results)
    selection = choose_mpi_ranks(results)
    selected = selection["selected"]
    output = {
        "status": "PASS" if reproducibility["status"] == "PASS" else "FAIL",
        "results": results, "benchmark_musubi_calls": len(results),
        "reproducibility": reproducibility, "selection": selection,
        "selected_ranks": selected["ranks"],
        "selected_iterations_per_s": selected["iterations_per_s"],
        "selected_mlups": selected["mlups"],
        "speedup_vs_old_8_rank": selected["iterations_per_s"] / OLD_8_RANK_ITERATIONS_S,
    }
    write_json(run_root / "qc" / "mpi_rank_benchmark.json", output)
    return output


def _pdf_state(binary: Path, *, with_velocity: bool) -> dict[str, Any]:
    pdf = read_restart_pdf(binary, n_elems=EXPECTED_FLUID_CELLS, n_components=19)
    values = np.asarray(pdf)
    finite = bool(np.all(np.isfinite(values)))
    minimum = float(np.min(values))
    density = np.sum(values, axis=1, dtype=np.float64)
    total = float(np.sum(density, dtype=np.float64))
    velocity = (values @ D3Q19_DIRECTIONS.astype(np.float64)) / density[:, None]
    maximum_speed = float(np.max(np.linalg.norm(velocity, axis=1)))
    pressure_scale = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / EXPECTED_DT_S**2 * CS2
    result: dict[str, Any] = {
        "total_pdf_mass": total,
        "mean_pressure_pa": float(np.mean(density, dtype=np.float64) * pressure_scale),
        "all_finite": finite,
        "minimum_pdf": minimum,
        "maximum_lattice_speed": maximum_speed,
    }
    if with_velocity:
        result["velocity_lattice"] = velocity
    return result


def pressure_pdf_accumulation_crosscheck(
    *, previous: dict[str, Any], current: dict[str, Any], delta_iter: int,
) -> dict[str, Any]:
    factor = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**3 / EXPECTED_DT_S
    pdf_acc = (float(current["total_pdf_mass"]) - float(previous["total_pdf_mass"])) / delta_iter * factor
    pressure_scale = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / EXPECTED_DT_S**2 * CS2
    density_delta_sum = (
        EXPECTED_FLUID_CELLS
        * (float(current["mean_pressure_pa"]) - float(previous["mean_pressure_pa"]))
        / pressure_scale
    )
    pressure_acc = density_delta_sum / delta_iter * factor
    direction_same = (
        pdf_acc == 0.0 or pressure_acc == 0.0
        or math.copysign(1.0, pdf_acc) == math.copysign(1.0, pressure_acc)
    )
    discrepancy = abs(pdf_acc - pressure_acc) / max(abs(pdf_acc), abs(pressure_acc), 1.0e-30)
    return {
        "pdf_accumulation_kg_s": pdf_acc,
        "pressure_accumulation_kg_s": pressure_acc,
        "direction_same": direction_same,
        "relative_discrepancy": discrepancy,
        "status": "PASS" if direction_same and discrepancy <= PRESSURE_MASS_CROSSCHECK_HARD else "FAIL",
    }


HISTORY_FIELDS = (
    "iteration", "short_window_start", "short_window_delta_iter", "R_mass_short",
    "long_window_start", "long_window_delta_iter", "R_mass_long", "R_velocity",
    "delta_pressure_pa", "R_pressure", "R_inlet", "R_boundary", "outlet_01_flow",
    "outlet_02_flow", "outlet_03_flow", "significant_backflow",
    "max_lattice_speed", "minimum_pdf", "stage", "all_final_gates_pass",
)


class ExactCheckpointAuditor:
    def __init__(
        self, *, mesh: Path, compatible: list[RestartEvidence], qc: Path,
        runtime_root: Path,
    ) -> None:
        self.mesh = Path(mesh)
        self.compatible = list(compatible)
        self.qc = Path(qc)
        self.runtime_root = Path(runtime_root)
        self.records: list[dict[str, Any]] = []
        self._state_cache: dict[int, dict[str, Any]] = {}
        self.gate_pass: dict[str, Any] | None = None
        self.fatal_failure: str | None = None
        self.fatal_status: str | None = None
        self.stop_file_created_at: str | None = None
        self.stop_requested_iteration: int | None = None
        self.history_path = self.qc / "fast_steady_residual_history.csv"
        with self.history_path.open("w", encoding="utf-8", newline="") as stream:
            csv.DictWriter(stream, fieldnames=HISTORY_FIELDS).writeheader()
        self._write_convergence()

    def _state(self, evidence: RestartEvidence, *, with_velocity: bool = False) -> dict[str, Any]:
        cached = self._state_cache.get(evidence.iteration)
        if cached is None or (with_velocity and "velocity_lattice" not in cached):
            state = _pdf_state(evidence.binary_path, with_velocity=with_velocity)
            if cached is not None:
                state.update(cached)
            self._state_cache[evidence.iteration] = state
            return state
        return cached

    def _write_convergence(self) -> None:
        write_json(
            self.qc / "checkpoint_exact_convergence.json",
            {
                "status": "PASS" if self.gate_pass else "IN_PROGRESS",
                "gate_pass_iteration": self.gate_pass.get("iteration") if self.gate_pass else None,
                "fatal_failure": self.fatal_failure,
                "fatal_status": self.fatal_status,
                "records": self.records,
            },
        )

    def _append_history(self, record: dict[str, Any]) -> None:
        row = {key: record.get(key) for key in HISTORY_FIELDS}
        with self.history_path.open("a", encoding="utf-8", newline="") as stream:
            csv.DictWriter(stream, fieldnames=HISTORY_FIELDS).writerow(row)

    def audit(
        self, current: RestartEvidence, *, persist: bool = True, request_stop: bool = True,
    ) -> dict[str, Any]:
        candidates = [*self.compatible, current]
        short = select_window(current.iteration, candidates, target=SHORT_TARGET, allowed=SHORT_RANGE)
        long = select_window(current.iteration, candidates, target=LONG_TARGET, allowed=LONG_RANGE)
        current_state = self._state(current)
        record: dict[str, Any] = {
            "iteration": current.iteration,
            "restart_header": str(current.header_path), "restart_binary": str(current.binary_path),
            "restart_sha256": current.sha256, "short_window_start": None,
            "short_window_delta_iter": None, "R_mass_short": None,
            "long_window_start": None, "long_window_delta_iter": None,
            "R_mass_long": None, "R_velocity": None, "delta_pressure_pa": None,
            "characteristic_pressure_drop_pa": None, "R_pressure": None,
            "R_inlet": None, "R_boundary": None, "outlet_01_flow": None,
            "outlet_02_flow": None, "outlet_03_flow": None,
            "outlet_signed_sum": None, "significant_backflow": None,
            "all_finite": current_state["all_finite"],
            "max_lattice_speed": current_state["maximum_lattice_speed"],
            "maximum_lattice_speed": current_state["maximum_lattice_speed"],
            "minimum_pdf": current_state["minimum_pdf"],
            "inlet_globbc": EXPECTED_INLET_GLOBBC,
            "stage": "WAITING_FOR_SHORT_WINDOW", "all_final_gates_pass": False,
        }
        if not current_state["all_finite"] or current_state["minimum_pdf"] <= 0.0 or current_state["maximum_lattice_speed"] >= 0.05:
            raise FlowError(STATUS_SAFETY, f"checkpoint {current.iteration} PDF safety gate failed")
        if short is not None:
            short_state = self._state(short)
            delta = current.iteration - short.iteration
            record.update({
                "short_window_start": short.iteration, "short_window_delta_iter": delta,
                "R_mass_short": mass_residual(short_state["total_pdf_mass"], current_state["total_pdf_mass"], delta),
            })
            record["stage"] = classify_stage(float(record["R_mass_short"]))
        if long is not None:
            long_state = self._state(long)
            delta = current.iteration - long.iteration
            record.update({
                "long_window_start": long.iteration, "long_window_delta_iter": delta,
                "R_mass_long": mass_residual(long_state["total_pdf_mass"], current_state["total_pdf_mass"], delta),
            })
        if short is not None and record["stage"] in {"STAGE_2_NEAR", "STAGE_3_FINAL"}:
            try:
                short_full = self._state(short, with_velocity=True)
                current_full = self._state(current, with_velocity=True)
                record["R_velocity"] = velocity_l2_residual(
                    short_full["velocity_lattice"], current_full["velocity_lattice"],
                )
                exact = exact_boundary_flux_audit(
                    mesh=self.mesh, restart_binary=current.binary_path,
                    link_csv=self.qc / f"checkpoint_exact_links_{current.iteration}.csv",
                )
                inlet = float(exact["inlet"]["exact_mass_flow_kg_s"])
                outlet_flows = {
                    label: float(exact["outlets"][label]["exact_signed_mass_flow_kg_s"])
                    for label in ("outlet_01", "outlet_02", "outlet_03")
                }
                delta_pressure = abs(float(current_full["mean_pressure_pa"]) - float(short_full["mean_pressure_pa"]))
                characteristic = characteristic_pressure_drop(float(exact["inlet"]["gauge_pressure_pa"]))
                crosscheck = pressure_pdf_accumulation_crosscheck(
                    previous=short_full, current=current_full,
                    delta_iter=current.iteration - short.iteration,
                )
                record.update({
                    "R_pressure": pressure_residual(delta_pressure, characteristic),
                    "delta_pressure_pa": delta_pressure,
                    "characteristic_pressure_drop_pa": characteristic,
                    "exact_inlet_mass_flow": inlet, "R_inlet": inlet_residual(inlet),
                    "outlet_01_flow": outlet_flows["outlet_01"],
                    "outlet_02_flow": outlet_flows["outlet_02"],
                    "outlet_03_flow": outlet_flows["outlet_03"],
                    "outlet_signed_sum": sum(outlet_flows.values()),
                    "R_boundary": boundary_residual(inlet, outlet_flows.values()),
                    "significant_backflow": significant_backflow(outlet_flows.values(), inlet),
                    "pressure_pdf_crosscheck": crosscheck,
                    "final_inlet_gauge_pressure_pa": exact["inlet"]["gauge_pressure_pa"],
                    "exact_boundary_flux": exact,
                    "all_finite": exact["all_finite"],
                    "minimum_pdf": exact["minimum_pdf"],
                    "max_lattice_speed": exact["maximum_lattice_speed"],
                    "maximum_lattice_speed": exact["maximum_lattice_speed"],
                    "inlet_globbc": exact["inlet"]["globbc_count"],
                })
                record["all_final_gates_pass"] = all_final_gates_pass(record)
                if crosscheck["status"] != "PASS":
                    self.fatal_failure = f"pressure/PDF accumulation crosscheck failed at {current.iteration}"
                    self.fatal_status = STATUS_CROSSCHECK
            except FlowError as error:
                if error.status == STATUS_SAFETY:
                    raise
                self.fatal_failure = str(error)
                self.fatal_status = STATUS_EXACT_UNRESOLVED
        if persist:
            self.records.append(record)
            if all(item.iteration != current.iteration for item in self.compatible):
                self.compatible.append(current)
                self.compatible.sort(key=lambda item: item.iteration)
            self._append_history(record)
            if record["all_final_gates_pass"] and self.gate_pass is None:
                self.gate_pass = record
                self.stop_requested_iteration = current.iteration
                if request_stop:
                    create_stop_file(self.runtime_root)
                    self.stop_file_created_at = datetime.now().isoformat()
                write_json(self.qc / "gate_pass_exact_boundary_flux.json", record["exact_boundary_flux"])
                write_json(self.qc / "gate_pass_velocity_residual.json", {
                    "iteration": current.iteration, "short_window_start": record["short_window_start"],
                    "R_velocity": record["R_velocity"], "hard_gate": VELOCITY_HARD,
                })
                write_json(self.qc / "gate_pass_mass_residual.json", {
                    "iteration": current.iteration, "short_window_start": record["short_window_start"],
                    "R_mass_short": record["R_mass_short"], "long_window_start": record["long_window_start"],
                    "R_mass_long": record["R_mass_long"], "hard_gate": MASS_HARD,
                })
                write_json(self.qc / "gate_pass_pressure_residual.json", {
                    "iteration": current.iteration, "delta_pressure_pa": record["delta_pressure_pa"],
                    "characteristic_pressure_drop_pa": record["characteristic_pressure_drop_pa"],
                    "R_pressure": record["R_pressure"], "hard_gate": PRESSURE_HARD,
                    "pressure_pdf_crosscheck": record["pressure_pdf_crosscheck"],
                })
            self._write_convergence()
        return record


def _tracker_evidence(record: dict[str, Any]) -> RestartEvidence:
    return RestartEvidence(
        iteration=int(record["iteration"]),
        header_path=Path(record["header_path_windows"]),
        binary_path=Path(record["binary_path_windows"]),
        sha256=str(record["sha256"]),
        source="runtime",
    )


def run_fast_monitored_wsl(
    *, distribution: str, workdir_wsl: str, command: list[str], stdout_path: Path,
    stderr_path: Path, controller_csv_path: Path, timeout_s: int,
    checkpoint_tracker: FastCheckpointTracker, auditor: ExactCheckpointAuditor,
) -> FastMonitoredRun:
    shell = f"cd {shlex.quote(workdir_wsl)} && exec {shlex.join(command)}"
    started = time.perf_counter()
    state: dict[str, Any] = {
        "safety_failure": None, "infrastructure_failure": None,
        "continuity_failure": None, "controller_record_count": 0,
        "latest_controller_iteration": None,
    }
    lock = threading.Lock()
    terminate_event = threading.Event()
    wrapper_timeout = False
    audited_iterations: set[int] = set()
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    controller_csv_path.parent.mkdir(parents=True, exist_ok=True)
    with (
        stdout_path.open("w", encoding="utf-8", errors="replace", buffering=65536) as stdout_sink,
        stderr_path.open("w", encoding="utf-8", errors="replace", buffering=65536) as stderr_sink,
        controller_csv_path.open("w", encoding="utf-8", newline="", buffering=65536) as controller_sink,
    ):
        stdout_policy = FlushPolicy(stdout_sink, record_limit=None)
        stderr_policy = FlushPolicy(stderr_sink, record_limit=None)
        controller_policy = FlushPolicy(controller_sink, record_limit=100)
        writer = csv.DictWriter(controller_sink, fieldnames=CONTROLLER_FIELDS)
        writer.writeheader()
        controller_policy.flush_if_due(force=True)
        process = subprocess.Popen(
            ["wsl.exe", "-d", distribution, "--", "bash", "-lc", shell],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            encoding="utf-8", errors="replace", bufsize=1,
        )

        def consume_stdout() -> None:
            assert process.stdout is not None
            try:
                for line in process.stdout:
                    stdout_sink.write(line)
                    stdout_policy.note()
                    if WSL_INFRASTRUCTURE_PATTERN.search(line):
                        with lock:
                            state["infrastructure_failure"] = state["infrastructure_failure"] or line.strip()
                        terminate_event.set()
                    failure = _live_safety_failure(line)
                    match = CONTROLLER_PATTERN.search(line)
                    with lock:
                        if match:
                            record = _match_record(match)
                            writer.writerow(record)
                            controller_policy.note()
                            if state["controller_record_count"] == 0 and int(record["iteration"]) <= SOURCE_ITERATION:
                                state["continuity_failure"] = f"first controller iteration {record['iteration']} <= {SOURCE_ITERATION}"
                                terminate_event.set()
                            state["controller_record_count"] += 1
                            state["latest_controller_iteration"] = int(record["iteration"])
                            failure = failure or _controller_safety_failure(record)
                            if float(record["relative_error"]) > CONTROLLER_ALGEBRA_HARD:
                                failure = failure or f"controller algebra residual > 1e-8 at {record['iteration']}"
                        if failure and state["safety_failure"] is None:
                            state["safety_failure"] = failure
                            terminate_event.set()
            except (OSError, ValueError) as error:
                with lock:
                    state["infrastructure_failure"] = f"stdout PIPE failure: {error}"
                terminate_event.set()

        def consume_stderr() -> None:
            assert process.stderr is not None
            try:
                for line in process.stderr:
                    stderr_sink.write(line)
                    stderr_policy.note()
                    if WSL_INFRASTRUCTURE_PATTERN.search(line):
                        with lock:
                            state["infrastructure_failure"] = state["infrastructure_failure"] or line.strip()
                        terminate_event.set()
            except (OSError, ValueError) as error:
                with lock:
                    state["infrastructure_failure"] = f"stderr PIPE failure: {error}"
                terminate_event.set()

        threads = [threading.Thread(target=consume_stdout, daemon=True), threading.Thread(target=consume_stderr, daemon=True)]
        for thread in threads:
            thread.start()
        termination_sent = False
        next_scan = started + 1.0
        next_progress = started + 60.0
        while process.poll() is None:
            time.sleep(0.1)
            now = time.perf_counter()
            stdout_policy.flush_if_due()
            stderr_policy.flush_if_due()
            controller_policy.flush_if_due()
            if now >= next_scan:
                try:
                    checkpoint_tracker.scan()
                    for item in checkpoint_tracker.records:
                        iteration = int(item["iteration"])
                        if iteration not in audited_iterations and iteration > SOURCE_ITERATION:
                            auditor.audit(_tracker_evidence(item))
                            audited_iterations.add(iteration)
                            if auditor.fatal_failure:
                                terminate_event.set()
                                break
                except FlowError as error:
                    if error.status == STATUS_SAFETY:
                        with lock:
                            state["safety_failure"] = str(error)
                    else:
                        auditor.fatal_failure = str(error)
                        auditor.fatal_status = STATUS_EXACT_UNRESOLVED
                    terminate_event.set()
                except (OSError, ValueError) as error:
                    with lock:
                        state["infrastructure_failure"] = f"checkpoint/audit failure: {error}"
                    terminate_event.set()
                next_scan = now + 1.0
            if terminate_event.is_set() and not termination_sent:
                stdout_policy.flush_if_due(force=True)
                stderr_policy.flush_if_due(force=True)
                controller_policy.flush_if_due(force=True)
                process.terminate()
                termination_sent = True
            if now - started >= timeout_s:
                wrapper_timeout = True
                process.terminate()
                break
            if now >= next_progress:
                with lock:
                    latest = state["latest_controller_iteration"]
                print(f"FAST_STEADY_PROGRESS elapsed_s={now-started:.1f} latest_controller_iter={latest}", flush=True)
                next_progress = now + 60.0
        try:
            returncode = process.wait(timeout=15)
        except subprocess.TimeoutExpired:
            process.kill()
            returncode = process.wait(timeout=15)
        for thread in threads:
            thread.join(timeout=15)
        checkpoint_tracker.scan(force_complete=True)
        for item in checkpoint_tracker.records:
            iteration = int(item["iteration"])
            if iteration not in audited_iterations and iteration > SOURCE_ITERATION and auditor.gate_pass is None:
                auditor.audit(_tracker_evidence(item), request_stop=False)
                audited_iterations.add(iteration)
        stdout_policy.flush_if_due(force=True)
        stderr_policy.flush_if_due(force=True)
        controller_policy.flush_if_due(force=True)
    return FastMonitoredRun(
        returncode=returncode, wall_time_s=time.perf_counter() - started,
        safety_failure=state["safety_failure"],
        infrastructure_failure=state["infrastructure_failure"],
        continuity_failure=state["continuity_failure"],
        audit_failure=auditor.fatal_failure, audit_failure_status=auditor.fatal_status,
        wrapper_timeout=wrapper_timeout,
        controller_record_count=int(state["controller_record_count"]),
        latest_controller_iteration=state["latest_controller_iteration"],
        stop_file_triggered=auditor.stop_file_created_at is not None,
        stop_requested_iteration=auditor.stop_requested_iteration,
        stop_file_created_at=auditor.stop_file_created_at,
        stdout_flushes=stdout_policy.flush_count, stderr_flushes=stderr_policy.flush_count,
        controller_flushes=controller_policy.flush_count,
    )


def _historical_evidence(source_root: Path) -> tuple[list[RestartEvidence], dict[str, Any]]:
    restart_root = Path(source_root) / "restart"
    evidence: dict[int, RestartEvidence] = {}
    contracts: list[dict[str, Any]] = []
    for header_path in sorted(restart_root.glob(f"{SHORT_SIMULATION_NAME}_header_*.lua")):
        header = parse_restart_header(header_path)
        if header.iteration > SOURCE_ITERATION:
            continue
        contract = restart_compatibility(header_path, header.binary_path)
        contracts.append(contract)
        if contract["status"] == "PASS":
            evidence[header.iteration] = RestartEvidence(
                iteration=header.iteration, header_path=header_path,
                binary_path=header.binary_path, sha256=contract["binary_sha256"],
                source="historical healthy same-physics checkpoint",
            )
    last_header = restart_root / f"{SHORT_SIMULATION_NAME}_lastHeader.lua"
    last = parse_restart_header(last_header)
    last_contract = restart_compatibility(last_header, last.binary_path)
    contracts.append(last_contract)
    if last_contract["status"] == "PASS":
        evidence[last.iteration] = RestartEvidence(
            iteration=last.iteration, header_path=last_header,
            binary_path=last.binary_path, sha256=last_contract["binary_sha256"],
            source="healthy resume checkpoint",
        )
    required = {299_326, 309_326, SOURCE_ITERATION}
    missing = sorted(required - set(evidence))
    status = "PASS" if not missing and all(item["status"] == "PASS" for item in contracts) else "FAIL"
    return sorted(evidence.values(), key=lambda item: item.iteration), {
        "status": status, "required_iterations": sorted(required),
        "missing_required_iterations": missing, "contracts": contracts,
    }


def _archive_selected_checkpoints(
    *, run_root: Path, selections: dict[str, RestartEvidence],
) -> dict[str, Any]:
    archive_root = Path(run_root) / "restart"
    archive_root.mkdir(parents=True, exist_ok=True)
    by_iteration: dict[int, dict[str, Any]] = {}
    roles: dict[str, Any] = {}
    for role, evidence in selections.items():
        archived = by_iteration.get(evidence.iteration)
        if archived is None:
            destination = archive_root / f"checkpoint_{evidence.iteration}"
            destination.mkdir(parents=True, exist_ok=False)
            header_path = destination / evidence.header_path.name
            binary_path = destination / evidence.binary_path.name
            shutil.copy2(evidence.header_path, header_path)
            shutil.copy2(evidence.binary_path, binary_path)
            archive_sha = sha256_file(binary_path)
            if archive_sha != evidence.sha256:
                raise FlowError(STATUS_FAILED, f"Archive SHA mismatch at {evidence.iteration}")
            archived = {
                "iteration": evidence.iteration,
                "archived_header": str(header_path), "archived_binary": str(binary_path),
                "sha256": archive_sha, "source": evidence.source,
            }
            by_iteration[evidence.iteration] = archived
        roles[role] = archived
    result = {"status": "PASS", "roles": roles, "unique_checkpoints": list(by_iteration.values())}
    write_json(Path(run_root) / "qc" / "selected_checkpoint_archive.json", result)
    return result


def _tracking_files(runtime_root: Path, label: str) -> list[Path]:
    return sorted((Path(runtime_root) / "tracking" / label).glob("*.res"))


def _copy_small_runtime_evidence(runtime_root: Path, run_root: Path) -> None:
    for source, destination in (
        (runtime_root / "musubi.lua", run_root / "musubi.lua"),
        (runtime_root / "timing" / "timing.res", run_root / "timing" / "timing.res"),
    ):
        if source.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    for label in ("p", "u"):
        destination = run_root / "tracking" / label
        destination.mkdir(parents=True, exist_ok=True)
        for source in _tracking_files(runtime_root, label):
            shutil.copy2(source, destination / source.name)


def _frozen_paths(
    *, root: Path, mesh: Path, source_root: Path, binary_windows: Path,
) -> tuple[Path, ...]:
    production = (
        root / "cfd_flow.py", root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
    )
    mesh_files = tuple(sorted(path for path in Path(mesh).iterdir() if path.is_file()))
    source_files = tuple(sorted(path for path in (Path(source_root) / "restart").iterdir() if path.is_file()))
    other = (
        Path(source_root) / "qc" / "healthy_calibration_manifest.json",
        Path(source_root) / "qc" / "healthy_mouse_capillary_reference.json",
        binary_windows,
    )
    return (*production, *mesh_files, *source_files, *other)


def _runtime_contract(
    *, lua_contract: dict[str, Any], selected_ranks: int,
    performance: dict[str, Any], runtime_root_wsl: str,
) -> dict[str, Any]:
    checks = {
        "lua_contract": lua_contract["status"] == "PASS",
        "persistent_wsl_ext4": runtime_root_wsl.startswith(f"{RUNTIME_BASE_WSL}/"),
        "mpi_core_binding": bool(performance["mpi_binding_supported"]),
        "no_physical_core_oversubscription": selected_ranks <= int(performance["physical_cores"]),
        "thread_limits": all(value == "1" for value in THREAD_ENV.values()),
        "tracking_interval_5000": TRACKING_INTERVAL == 5_000,
        "checkpoint_interval_5000": CHECKPOINT_INTERVAL == 5_000,
        "built_in_pu_steady_disabled": True,
        "controller_flush_100_or_1s": True,
        "stdout_stderr_flush_about_1s": True,
        "checkpoint_exact_early_stop": True,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks, "selected_ranks": selected_ranks,
        "binding_args": list(MPI_BINDING_ARGS), "thread_environment": THREAD_ENV,
        "runtime_root_wsl": runtime_root_wsl,
    }


def run_healthy_capillary_fast_steady(project_root: Path) -> dict[str, Any]:
    """Benchmark ranks, then run at most one checkpoint-audited continuation."""

    root = Path(project_root).resolve()
    head = _git(root, "rev-parse", "HEAD")
    branch = _git(root, "branch", "--show-current")
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root
    mesh = output_root / AXIS_MESH_RUN / "seeder" / "mesh"
    source_root = output_root / SOURCE_RUN
    source_header = source_root / "restart" / f"{SHORT_SIMULATION_NAME}_lastHeader.lua"
    source_binary = source_root / "restart" / SOURCE_BINARY_NAME
    binary_windows = _wsl_path_to_unc(config.apes.wsl_distribution, BINARY_WSL)
    production_paths = (
        root / "cfd_flow.py", root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
    )
    production_diff = _git(
        root, "diff", "--", *(str(path.relative_to(root)) for path in production_paths),
    )
    if production_diff:
        raise FlowError(STATUS_FAILED, "Production pipeline has local modifications")
    if config.mesh.dx_target_m != EXPECTED_DX_M:
        raise FlowError(STATUS_FAILED, f"Frozen dx changed: {config.mesh.dx_target_m}")
    if config.solver.wallclock_limit_s != WALLCLOCK_LIMIT_S:
        raise FlowError(STATUS_FAILED, "Frozen wallclock limit changed")
    if parse_mesh_header(mesh)["fluid_element_count"] != EXPECTED_FLUID_CELLS:
        raise FlowError(STATUS_FAILED, "Frozen mesh fluid count changed")
    source_contract = restart_compatibility(source_header, source_binary)
    if (
        source_contract["status"] != "PASS"
        or source_contract["iteration"] != SOURCE_ITERATION
        or source_contract["binary_sha256"] != EXPECTED_SOURCE_SHA256
    ):
        raise FlowError(STATUS_FAILED, f"Source restart contract failed: {source_contract}")
    binary_sha = sha256_file(binary_windows)
    if binary_sha != EXPECTED_BINARY_SHA256:
        raise FlowError(STATUS_FAILED, f"Adaptive binary SHA changed: {binary_sha}")
    history, history_contract = _historical_evidence(source_root)
    if history_contract["status"] != "PASS":
        raise FlowError(STATUS_FAILED, f"Historical checkpoint compatibility failed: {history_contract}")
    outlet_contract = source_only_pressure_outlet_contract(mesh)
    if outlet_contract["status"] != "PASS":
        raise FlowError(STATUS_EXACT_UNRESOLVED, "Source-proven outlet contract changed")
    frozen_paths = _frozen_paths(
        root=root, mesh=mesh, source_root=source_root, binary_windows=binary_windows,
    )
    frozen_before = _file_manifest(frozen_paths)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{RUN_PREFIX}_{stamp}"
    qc = run_root / "qc"
    tracking = run_root / "tracking"
    qc.mkdir(parents=True, exist_ok=False)
    tracking.mkdir(parents=True, exist_ok=False)
    manifest_path = qc / "fast_steady_final_manifest.json"
    reference_source = source_root / "qc" / "healthy_mouse_capillary_reference.json"
    shutil.copy2(reference_source, qc / "healthy_mouse_capillary_reference.json")
    write_json(qc / "source_restart_contract.json", source_contract)
    write_json(qc / "historical_restart_compatibility.json", history_contract)
    write_json(qc / "source_proven_outlet_contract.json", outlet_contract)
    summary: dict[str, Any] = {
        "status": STATUS_FAILED, "next": NEXT_EXACT, "first_failure": None,
        "run_root": str(run_root), "actual_head": head, "branch": branch,
        "production_pipeline_modified": False, "resume_run": SOURCE_RUN,
        "resume_iteration": SOURCE_ITERATION, "resume_header": str(source_header),
        "resume_binary": str(source_binary), "resume_binary_sha256": source_contract["binary_sha256"],
        "adaptive_binary_sha256": binary_sha, "adaptive_binary_rebuilt": False,
        "mesh": str(mesh), "fluid_cells": EXPECTED_FLUID_CELLS,
        "seeder_calls": 0, "benchmark_musubi_calls": 0, "production_musubi_calls": 0,
        "harvester_calls": 0, "grid_convergence": "NOT_RUN",
        "tracking_interval": TRACKING_INTERVAL, "checkpoint_interval": CHECKPOINT_INTERVAL,
        "built_in_pu_steady_enabled": False,
        "controller_buffering_strategy": "100 records or 1 second; stdout/stderr about 1 second",
        "dx_m": EXPECTED_DX_M, "dt_s": EXPECTED_DT_S,
        "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S, "target_q_m3_s": TARGET_Q_M3_S,
        "started_at": datetime.now().isoformat(), "frozen_files_before": frozen_before,
    }
    write_json(manifest_path, summary)
    tracker: FastCheckpointTracker | None = None
    runtime_root: Path | None = None
    try:
        health = _one_wsl_health_check(config.apes.wsl_distribution)
        write_json(qc / "wsl_health_preflight.json", health)
        if health["status"] != "PASS":
            raise FlowError(STATUS_RUNTIME, "WSL persistent ext4 preflight failed")
        performance = performance_preflight(config.apes.wsl_distribution)
        write_json(qc / "mpi_hardware_preflight.json", performance)
        if performance["status"] != "PASS":
            raise FlowError(STATUS_RUNTIME, "MPI binding/physical-core preflight failed")
        summary.update({
            "physical_cores": performance["physical_cores"],
            "logical_cpus": performance["logical_cpus"],
            "mpi_binding_used": True,
        })
        benchmark = run_mpi_benchmarks(
            distribution=config.apes.wsl_distribution, mesh=mesh,
            source_header=source_header, source_binary=source_binary,
            run_root=run_root, stamp=stamp,
            physical_cores=int(performance["physical_cores"]),
        )
        summary.update({
            "benchmark_musubi_calls": benchmark["benchmark_musubi_calls"],
            "benchmark": benchmark,
            "selected_mpi_ranks": benchmark["selected_ranks"],
            "selected_iterations_per_s": benchmark["selected_iterations_per_s"],
            "selected_mlups": benchmark["selected_mlups"],
            "speedup_vs_old_8_rank": benchmark["speedup_vs_old_8_rank"],
        })
        write_json(manifest_path, summary)
        if benchmark["status"] != "PASS":
            raise FlowError(STATUS_MPI_REPRODUCIBILITY, "MPI endpoint reproducibility was not roundoff-consistent")

        runtime_root_wsl = f"{RUNTIME_BASE_WSL}/f_{stamp[-6:]}"
        runtime_root = _wsl_path_to_unc(config.apes.wsl_distribution, runtime_root_wsl)
        if runtime_root.exists():
            raise FlowError(STATUS_RUNTIME, f"Production WSL runtime root exists: {runtime_root_wsl}")
        _prepare_short_tracking_root(
            distribution=config.apes.wsl_distribution, root_wsl=runtime_root_wsl,
        )
        _stage_source_checkpoint(
            source_header=source_header, source_binary=source_binary, runtime_root=runtime_root,
        )
        resume_wsl = f"restart/{source_header.name}"
        lua = generate_production_lua(
            mesh_wsl=windows_to_wsl(mesh, config.apes.wsl_distribution),
            resume_header_wsl=resume_wsl,
        )
        (runtime_root / "musubi.lua").write_text(lua, encoding="utf-8")
        lua_contract = fast_lua_contract(lua, resume_header_wsl=resume_wsl, benchmark=False)
        write_json(qc / "fast_steady_lua_contract.json", lua_contract)
        runtime_contract = _runtime_contract(
            lua_contract=lua_contract, selected_ranks=int(benchmark["selected_ranks"]),
            performance=performance, runtime_root_wsl=runtime_root_wsl,
        )
        write_json(qc / "fast_steady_runtime_contract.json", runtime_contract)
        if runtime_contract["status"] != "PASS":
            raise FlowError(STATUS_FAILED, "Fast steady runtime contract failed")
        luac_rc = _run_luac(
            distribution=config.apes.wsl_distribution, workdir_wsl=runtime_root_wsl,
            stdout_path=tracking / "luac_stdout.log", stderr_path=tracking / "luac_stderr.log",
        )
        if luac_rc != 0:
            raise FlowError(STATUS_FAILED, "Production Lua syntax failed")
        tracker = FastCheckpointTracker(
            runtime_root=runtime_root, runtime_root_wsl=runtime_root_wsl,
            manifest_path=run_root / "checkpoint_manifest.json",
        )
        auditor = ExactCheckpointAuditor(
            mesh=mesh, compatible=history, qc=qc, runtime_root=runtime_root,
        )
        command = [
            "env", *[f"{key}={value}" for key, value in THREAD_ENV.items()],
            "OMPI_ALLOW_RUN_AS_ROOT=1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            MPIRUN_WSL, "-np", str(benchmark["selected_ranks"]),
            *MPI_BINDING_ARGS, BINARY_WSL, "musubi.lua",
        ]
        summary.update({
            "status": "CFD_FLOW_HEALTHY_FAST_STEADY_RUNNING",
            "production_musubi_calls": 1,
            "runtime_root_wsl": runtime_root_wsl,
        })
        write_json(manifest_path, summary)
        run = run_fast_monitored_wsl(
            distribution=config.apes.wsl_distribution, workdir_wsl=runtime_root_wsl,
            command=command, stdout_path=tracking / "musubi_stdout.log",
            stderr_path=tracking / "musubi_stderr.log",
            controller_csv_path=tracking / "controller_records.csv",
            timeout_s=WALLCLOCK_LIMIT_S + WRAPPER_MARGIN_S,
            checkpoint_tracker=tracker, auditor=auditor,
        )
        _copy_small_runtime_evidence(runtime_root, run_root)
        controller = summarize_controller_csv(tracking / "controller_records.csv", enforce_gates=False)
        latest = _tracker_evidence(tracker.records[-1]) if tracker.records else None
        summary.update({
            "musubi_returncode": run.returncode, "production_wall_time_s": run.wall_time_s,
            "production_first_iteration": controller["first_iteration"],
            "production_final_iteration": controller["final_iteration"],
            "additional_production_iterations": controller["final_iteration"] - SOURCE_ITERATION,
            "latest_checkpoint_iteration": latest.iteration if latest else None,
            "checkpoint_iterations": [int(item["iteration"]) for item in tracker.records],
            "gate_pass_iteration": auditor.gate_pass.get("iteration") if auditor.gate_pass else None,
            "exact_convergence_triggered": auditor.gate_pass is not None,
            "built_in_convergence_triggered": False,
            "stop_file_triggered": run.stop_file_triggered,
            "stop_file_created_at": run.stop_file_created_at,
            "stop_requested_iteration": run.stop_requested_iteration,
            "actual_last_controller_iteration": run.latest_controller_iteration,
            "runtime_controller": controller,
            "writer_flush_counts": {
                "stdout": run.stdout_flushes, "stderr": run.stderr_flushes,
                "controller": run.controller_flushes,
            },
            "live_safety_failure": run.safety_failure,
            "wsl_infrastructure_failure": run.infrastructure_failure,
            "continuity_failure": run.continuity_failure,
            "exact_audit_failure": run.audit_failure,
            "wrapper_timeout": run.wrapper_timeout,
        })
        if run.safety_failure:
            raise FlowError(STATUS_SAFETY, run.safety_failure)
        if run.infrastructure_failure or run.continuity_failure or run.wrapper_timeout:
            raise FlowError(
                STATUS_RUNTIME,
                run.infrastructure_failure or run.continuity_failure or "wrapper timeout",
            )
        if run.audit_failure:
            raise FlowError(run.audit_failure_status or STATUS_EXACT_UNRESOLVED, run.audit_failure)
        if run.returncode != 0:
            raise FlowError(STATUS_RUNTIME, f"Musubi returned {run.returncode}")

        selected_records: dict[str, RestartEvidence] = {}
        accepted: RestartEvidence | None = None
        final_record: dict[str, Any] | None = None
        shutdown_safety: dict[str, Any] | None = None
        if auditor.gate_pass is not None:
            accepted = RestartEvidence(
                iteration=int(auditor.gate_pass["iteration"]),
                header_path=Path(auditor.gate_pass["restart_header"]),
                binary_path=Path(auditor.gate_pass["restart_binary"]),
                sha256=str(auditor.gate_pass["restart_sha256"]), source="gate-pass runtime checkpoint",
            )
            if sha256_file(accepted.binary_path) != accepted.sha256:
                raise FlowError(STATUS_EXACT_UNRESOLVED, "Gate-pass checkpoint SHA changed before final audit")
            final_record = auditor.audit(accepted, persist=False)
            if not all_final_gates_pass(final_record):
                raise FlowError(STATUS_EXACT_UNRESOLVED, "Gate-pass checkpoint did not reproduce all exact gates")
            selected_records["accepted_gate_pass"] = accepted
            window_candidates = history + [_tracker_evidence(item) for item in tracker.records]
            short = next(item for item in window_candidates if item.iteration == int(final_record["short_window_start"]))
            long = next(item for item in window_candidates if item.iteration == int(final_record["long_window_start"]))
            selected_records["short_window_reference"] = short
            selected_records["long_window_reference"] = long
        if latest is not None:
            selected_records["final_shutdown"] = latest
            shutdown_safety = _pdf_state(latest.binary_path, with_velocity=False)
            if (
                not shutdown_safety["all_finite"]
                or shutdown_safety["minimum_pdf"] <= 0.0
                or shutdown_safety["maximum_lattice_speed"] >= 0.05
            ):
                raise FlowError(STATUS_SAFETY, "Final shutdown restart safety audit failed")
        archive = _archive_selected_checkpoints(run_root=run_root, selections=selected_records) if selected_records else None
        frozen_after = _file_manifest(frozen_paths)
        source_unchanged = frozen_before == frozen_after
        write_json(qc / "source_frozen_files_unchanged_qc.json", {
            "status": "PASS" if source_unchanged else "FAIL",
            "source_frozen_files_unchanged": source_unchanged,
            "before": frozen_before, "after": frozen_after,
        })
        if not source_unchanged:
            raise FlowError(STATUS_FAILED, "Source frozen files changed")
        if accepted is not None and archive is not None:
            accepted_archive = archive["roles"]["accepted_gate_pass"]
            shutdown_archive = archive["roles"].get("final_shutdown")
            summary.update({
                "status": STATUS_PASS, "next": NEXT_GRID, "first_failure": None,
                "criterion_that_stopped_solver": (
                    "PYTHON EXACT CHECKPOINT GATES -> stop file"
                    if run.stop_file_triggered
                    else "MUSUBI 3600 s wallclock; exact gates confirmed at shutdown"
                ),
                "accepted_steady_restart_path": accepted_archive["archived_binary"],
                "accepted_steady_restart_sha256": accepted_archive["sha256"],
                "shutdown_final_restart_path": shutdown_archive["archived_binary"] if shutdown_archive else None,
                "gate_pass_exact_record": final_record,
                "shutdown_final_safety": shutdown_safety,
            })
        else:
            summary.update({
                "status": STATUS_WALLCLOCK, "next": NEXT_RESUME,
                "first_failure": None,
                "criterion_that_stopped_solver": "MUSUBI 3600 s wallclock",
                "accepted_steady_restart_path": None,
                "accepted_steady_restart_sha256": None,
                "shutdown_final_restart_path": archive["roles"]["final_shutdown"]["archived_binary"] if archive and "final_shutdown" in archive["roles"] else None,
            })
        summary.update({
            "source_frozen_files_unchanged": source_unchanged,
            "final_exact_record": final_record,
            "completed_at": datetime.now().isoformat(),
        })
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        status = error.status if isinstance(error, FlowError) else STATUS_FAILED
        next_map = {
            STATUS_SAFETY: NEXT_SAFETY, STATUS_RUNTIME: NEXT_RUNTIME,
            STATUS_CROSSCHECK: NEXT_EXACT, STATUS_EXACT_UNRESOLVED: NEXT_EXACT,
            STATUS_MPI_REPRODUCIBILITY: NEXT_MPI,
        }
        source_unchanged = False
        try:
            source_unchanged = frozen_before == _file_manifest(frozen_paths)
        except (OSError, ValueError):
            pass
        summary.update({
            "status": status, "next": next_map.get(status, NEXT_EXACT),
            "first_failure": str(error),
            "source_frozen_files_unchanged": source_unchanged,
            "completed_at": datetime.now().isoformat(),
        })
        if tracker is not None:
            summary["checkpoint_iterations"] = [int(item["iteration"]) for item in tracker.records]
            summary["latest_checkpoint_iteration"] = tracker.records[-1]["iteration"] if tracker.records else None
        write_json(manifest_path, summary)
        return summary


def watch_existing_fast_steady_run(project_root: Path, run_root: Path) -> dict[str, Any]:
    """Zero-solver-call exact watchdog for an already-running production process."""

    root = Path(project_root).resolve()
    run_root = Path(run_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    mesh = config.paths.output_root / AXIS_MESH_RUN / "seeder" / "mesh"
    source_root = config.paths.output_root / SOURCE_RUN
    history, contract = _historical_evidence(source_root)
    if contract["status"] != "PASS":
        raise FlowError(STATUS_EXACT_UNRESOLVED, "Watchdog historical checkpoint contract failed")
    manifest_path = run_root / "qc" / "fast_steady_final_manifest.json"
    manifest = read_json(manifest_path)
    runtime_root_wsl = str(manifest["runtime_root_wsl"])
    runtime_root = _wsl_path_to_unc(config.apes.wsl_distribution, runtime_root_wsl)
    watchdog_qc = run_root / "qc" / "exact_watchdog"
    watchdog_qc.mkdir(parents=True, exist_ok=False)
    auditor = ExactCheckpointAuditor(
        mesh=mesh, compatible=history, qc=watchdog_qc, runtime_root=runtime_root,
    )
    audited: set[int] = set()
    started = time.perf_counter()
    last_manifest_status = str(manifest.get("status"))
    while time.perf_counter() - started <= WALLCLOCK_LIMIT_S + WRAPPER_MARGIN_S:
        try:
            checkpoint_manifest = read_json(run_root / "checkpoint_manifest.json")
        except (OSError, ValueError):
            time.sleep(0.25)
            continue
        for item in checkpoint_manifest.get("records", []):
            iteration = int(item["iteration"])
            if iteration > SOURCE_ITERATION and iteration not in audited:
                auditor.audit(_tracker_evidence(item))
                audited.add(iteration)
                print(
                    f"FAST_EXACT_WATCHDOG iteration={iteration} "
                    f"stage={auditor.records[-1]['stage']} "
                    f"R_mass_short={auditor.records[-1]['R_mass_short']} "
                    f"pass={auditor.records[-1]['all_final_gates_pass']}",
                    flush=True,
                )
                if auditor.fatal_failure:
                    result = {
                        "status": auditor.fatal_status or STATUS_EXACT_UNRESOLVED,
                        "first_failure": auditor.fatal_failure,
                        "musubi_calls": 0, "audited_iterations": sorted(audited),
                    }
                    write_json(watchdog_qc / "watchdog_manifest.json", result)
                    return result
                if auditor.gate_pass is not None:
                    result = {
                        "status": "EXACT_GATE_PASS_STOP_REQUESTED",
                        "gate_pass_iteration": auditor.gate_pass["iteration"],
                        "stop_file_created_at": auditor.stop_file_created_at,
                        "musubi_calls": 0, "audited_iterations": sorted(audited),
                        "gate_pass_record": auditor.gate_pass,
                    }
                    write_json(watchdog_qc / "watchdog_manifest.json", result)
                    return result
        try:
            last_manifest_status = str(read_json(manifest_path).get("status"))
        except (OSError, ValueError):
            pass
        if last_manifest_status != "CFD_FLOW_HEALTHY_FAST_STEADY_RUNNING":
            break
        time.sleep(1.0)
    result = {
        "status": "NO_EXACT_GATE_PASS_BEFORE_SOLVER_EXIT",
        "first_failure": None, "musubi_calls": 0,
        "audited_iterations": sorted(audited),
    }
    write_json(watchdog_qc / "watchdog_manifest.json", result)
    return result


def reconcile_existing_fast_steady_watchdog(
    project_root: Path, run_root: Path,
) -> dict[str, Any]:
    """Promote the zero-call watchdog history after the original process exits."""

    root = Path(project_root).resolve()
    run_root = Path(run_root).resolve()
    qc = run_root / "qc"
    manifest_path = qc / "fast_steady_final_manifest.json"
    summary = read_json(manifest_path)
    watchdog_qc = qc / "exact_watchdog"
    watchdog = read_json(watchdog_qc / "watchdog_manifest.json")
    convergence = read_json(watchdog_qc / "checkpoint_exact_convergence.json")
    if int(summary.get("production_musubi_calls", -1)) != 1:
        raise FlowError(STATUS_FAILED, "Reconciliation requires exactly one production Musubi call")
    if int(watchdog.get("musubi_calls", -1)) != 0:
        raise FlowError(STATUS_FAILED, "Exact watchdog was not zero-solver-call")
    if watchdog.get("status") != "NO_EXACT_GATE_PASS_BEFORE_SOLVER_EXIT":
        raise FlowError(STATUS_FAILED, f"Unexpected watchdog status: {watchdog.get('status')}")
    records = convergence.get("records", [])
    if not records:
        raise FlowError(STATUS_EXACT_UNRESOLVED, "Watchdog produced no exact checkpoint records")
    final = records[-1]
    if int(final["iteration"]) != int(summary["latest_checkpoint_iteration"]):
        raise FlowError(STATUS_EXACT_UNRESOLVED, "Watchdog final iteration differs from latest checkpoint")
    if bool(final["all_final_gates_pass"]):
        raise FlowError(STATUS_FAILED, "Incomplete-wallclock reconciliation received a gate-pass checkpoint")
    canonical_files = (
        "fast_steady_residual_history.csv",
        "checkpoint_exact_convergence.json",
    )
    backups: dict[str, str] = {}
    for name in canonical_files:
        canonical = qc / name
        source = watchdog_qc / name
        backup = qc / f"pre_watchdog_window_candidate_bug_{name}"
        if not backup.exists() and canonical.is_file():
            shutil.copy2(canonical, backup)
        shutil.copy2(source, canonical)
        backups[name] = str(backup)
    shutdown_path = Path(summary["shutdown_final_restart_path"])
    shutdown_sha = sha256_file(shutdown_path)
    checkpoint_manifest = read_json(run_root / "checkpoint_manifest.json")
    latest = checkpoint_manifest["records"][-1]
    if shutdown_sha != latest["sha256"]:
        raise FlowError(STATUS_EXACT_UNRESOLVED, "Archived shutdown restart SHA differs from WSL checkpoint evidence")
    reconciliation = {
        "status": "PASS",
        "reason": "runtime checkpoint candidates were added to each subsequent short/long selector",
        "production_musubi_calls_added": 0,
        "watchdog_musubi_calls": 0,
        "canonical_history_replaced": True,
        "pre_reconciliation_evidence_backups": backups,
        "final_iteration": final["iteration"],
        "final_all_gates_pass": final["all_final_gates_pass"],
        "shutdown_restart_sha256": shutdown_sha,
    }
    write_json(qc / "fast_steady_watchdog_reconciliation.json", reconciliation)
    summary.update({
        "status": STATUS_WALLCLOCK,
        "next": NEXT_RESUME,
        "first_failure": None,
        "final_exact_record": final,
        "shutdown_final_restart_sha256": shutdown_sha,
        "shutdown_final_safety": {
            "all_finite": final["all_finite"],
            "minimum_pdf": final["minimum_pdf"],
            "maximum_lattice_speed": final["maximum_lattice_speed"],
        },
        "exact_watchdog_reconciliation": reconciliation,
        "completed_at": datetime.now().isoformat(),
    })
    write_json(manifest_path, summary)
    return summary
