"""One isolated adaptive-flux Musubi steady baseline on the frozen mesh."""

from __future__ import annotations

import csv
import math
import re
import shlex
import shutil
import subprocess
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

from .adaptive_flux_validation import (
    BINARY_WSL,
    CONTROLLER_PATTERN,
    EXPECTED_BINARY_SHA256,
    EXPECTED_DT_S,
    MPIRUN_WSL,
)
from .apes import parse_mesh_header, windows_to_wsl
from .config import load_cfd_flow_config
from .diagnostics import parse_official_steady_termination
from .exact_link_flux import (
    EXPECTED_DX_M,
    REFERENCE_DENSITY_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    _file_manifest,
)
from .io import FlowError, sha256_file, write_json
from .restart_decode import parse_restart_header, restart_binary_size_contract


STEADY_PREFIX = "mcclure_adaptive_flux_steady_anchor003274"
AXIS_MESH_RUN = "axis_aligned_ideal_plane_inlet_preflight_anchor003274_20260829_120444"
STEADY_PENDING_AUDIT = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_STEADY_PASS_PENDING_EXACT_AUDIT"
STEADY_FAILED = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_STEADY_FAILED"
NEXT_EXACT_AUDIT = "RUN INDEPENDENT FINAL-RESTART PDF-LINK AUDIT"
NEXT_REVIEW_STEADY = "REVIEW ADAPTIVE-FLUX STEADY FAILURE"

EXPECTED_FLUID_CELLS = 221_309
EXPECTED_INLET_GLOBBC = 287
MPI_RANKS = 8
PRESSURE_REFERENCE_PA = 23622.32012800001
NU_M2_S = 3.27e-6
BULK_NU_M2_S = 2.18e-6
PRESSURE_THRESHOLD_PA = 0.14590517589648699
VELOCITY_THRESHOLD_M_S = 9.8385580075361729e-08
CONVERGENCE_INTERVAL = 100
CONVERGENCE_NVALS = 100
WRAPPER_MARGIN_S = 300
CHECKPOINT_INTERVAL = 20_000
EXPECTED_RESTART_BYTES = EXPECTED_FLUID_CELLS * 19 * 8
RUNTIME_BASE_WSL = "/home/lzy/u3da"
SHORT_SIMULATION_NAME = "a3274"
SYNTHETIC_RECORDS = 30_000
WSL_PREFLIGHT_FAILED = "CFD_FLOW_WSL_RUNTIME_IO_PREFLIGHT_FAILED"
STREAM_PREFLIGHT_FAILED = "CFD_FLOW_WSL_STREAMING_MONITOR_PREFLIGHT_FAILED"
WSL_RUNTIME_UNSTABLE = "CFD_FLOW_WSL_RUNTIME_UNSTABLE"
NEXT_REPAIR_WSL = "REPAIR WSL BEFORE CFD RETRY"
NEXT_REVIEW_WSL = "REVIEW WSL HOST STABILITY OR RESUME FROM LAST CHECKPOINT"

CONTROLLER_FIELDS = (
    "iteration",
    "target_lattice",
    "controlled_lattice",
    "relative_error",
    "rho_boundary",
    "pressure_pa",
    "max_lattice_velocity",
    "minimum_pdf",
    "globBC_count",
)


@dataclass(frozen=True, slots=True)
class MonitoredRun:
    returncode: int
    wall_time_s: float
    safety_failure: str | None
    wrapper_timeout: bool
    infrastructure_failure: str | None
    continuity_failure: str | None
    controller_record_count: int
    latest_controller_iteration: int | None


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", "-C", str(root), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.strip()


def generate_adaptive_steady_lua(
    *,
    mesh_wsl: str,
    restart_wsl: str,
    pressure_tracking_wsl: str,
    velocity_tracking_wsl: str,
    maximum_iterations: int,
    wallclock_limit_s: int,
) -> str:
    """Preserve the accepted project steady criterion and change only the inlet BC."""

    return f"""-- Isolated adaptive-flux steady baseline; fresh initialization.
simulation_name = '{SHORT_SIMULATION_NAME}'
printRuntimeInfo = true
timing_file = 'timing/timing.res'
mesh = '{mesh_wsl}/'
scaling = 'diffusive'
logging = {{ level = 5 }}

dx = {EXPECTED_DX_M:.17g}
dt = {EXPECTED_DT_S:.17g}
rho0_phy = {REFERENCE_DENSITY_KG_M3:.17g}
nu_phy = {NU_M2_S:.17g}
bulk_viscosity_phy = {BULK_NU_M2_S:.17g}
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}
maximum_iterations = {int(maximum_iterations)}

function outlet_01_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + 14.544978101274268:.17g} end
function outlet_02_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + 132.20454922317552:.17g} end
function outlet_03_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA - 13.700626673311461:.17g} end

sim_control = {{
  time_control = {{
    max = {{ iter = maximum_iterations, clock = {int(wallclock_limit_s)} }},
    interval = {{ iter = {CONVERGENCE_INTERVAL} }}
  }},
  abort_criteria = {{
    stop_file = 'stop',
    steady_state = true,
    convergence = {{
      variable = {{ 'pressure_phy', 'vel_mag_phy' }},
      shape = {{ kind = 'all' }},
      reduction = {{ 'average', 'average' }},
      time_control = {{ min = {{ iter = 0 }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
      norm = 'average',
      nvals = {CONVERGENCE_NVALS},
      absolute = true,
      condition = {{
        {{ threshold = {PRESSURE_THRESHOLD_PA:.17g}, operator = '<=' }},
        {{ threshold = {VELOCITY_THRESHOLD_M_S:.17g}, operator = '<=' }}
      }}
    }}
  }}
}}

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
}}

tracking = {{
  {{
    label = 'p', folder = '{pressure_tracking_wsl}/',
    variable = {{ 'pressure_phy' }}, shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = 0 }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }},
  {{
    label = 'u', folder = '{velocity_tracking_wsl}/',
    variable = {{ 'vel_mag_phy' }}, shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = 0 }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }}
}}

restart = {{
  write = '{restart_wsl}/',
  time_control = {{
    min = {{ iter = {CHECKPOINT_INTERVAL} }},
    max = {{ iter = maximum_iterations }},
    interval = {{ iter = {CHECKPOINT_INTERVAL} }}
  }}
}}
"""


def steady_lua_contract(text: str) -> dict[str, Any]:
    checks = {
        "fresh_initialization": "restart = {" in text
        and "write = 'restart/'" in text
        and "read =" not in text,
        "adaptive_inlet": "kind = 'adaptive_flux_pressure'" in text,
        "target_mass_flow": f"mass_flowrate = {TARGET_MASS_FLOW_KG_S:.17g}" in text,
        "pressure_outlets": text.count("kind = 'pressure_eq'") == 3,
        "wall_libb": "kind = 'wall_libb'" in text,
        "d3q19_bgk": "layout = 'd3q19'" in text and "relaxation = 'bgk'" in text,
        "criterion_variables": "variable = { 'pressure_phy', 'vel_mag_phy' }" in text,
        "criterion_average": "norm = 'average'" in text and "absolute = true" in text,
        "criterion_nvals": f"nvals = {CONVERGENCE_NVALS}" in text,
        "criterion_interval": f"interval = {{ iter = {CONVERGENCE_INTERVAL} }}" in text,
        "pressure_threshold": f"threshold = {PRESSURE_THRESHOLD_PA:.17g}" in text,
        "velocity_threshold": f"threshold = {VELOCITY_THRESHOLD_M_S:.17g}" in text,
        "restart_write": "write = 'restart/'" in text,
        "restart_interval": f"interval = {{ iter = {CHECKPOINT_INTERVAL} }}" in text,
        "no_full_field_export": "format = 'vtk'" not in text and "velocity_phy'" not in text,
        "no_harvester": "harvest" not in text.lower(),
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _match_record(match: re.Match[str]) -> dict[str, float | int]:
    return {
        "iteration": int(match.group("iteration")),
        "target_lattice": float(match.group("target")),
        "controlled_lattice": float(match.group("controlled")),
        "relative_error": float(match.group("error")),
        "rho_boundary": float(match.group("rho")),
        "pressure_pa": float(match.group("pressure")),
        "max_lattice_velocity": float(match.group("speed")),
        "minimum_pdf": float(match.group("minimum_pdf")),
        "globBC_count": int(match.group("count")),
    }


def summarize_runtime_controller(paths: Iterable[Path]) -> dict[str, Any]:
    records: list[dict[str, float | int]] = []
    for path in paths:
        with Path(path).open("r", encoding="utf-8", errors="replace") as stream:
            for line in stream:
                match = CONTROLLER_PATTERN.search(line)
                if match:
                    records.append(_match_record(match))
    if not records:
        raise FlowError(STEADY_FAILED, "No adaptive controller records were emitted")
    numeric_fields = (
        "target_lattice",
        "controlled_lattice",
        "relative_error",
        "rho_boundary",
        "pressure_pa",
        "max_lattice_velocity",
        "minimum_pdf",
    )
    if not all(
        math.isfinite(float(record[field]))
        for record in records
        for field in numeric_fields
    ):
        raise FlowError(STEADY_FAILED, "Adaptive controller emitted NaN/Inf")
    counts = sorted({int(record["globBC_count"]) for record in records})
    if counts != [EXPECTED_INLET_GLOBBC]:
        raise FlowError(STEADY_FAILED, f"Adaptive inlet globBC changed: {counts}")
    maximum_speed = max(float(record["max_lattice_velocity"]) for record in records)
    minimum_pdf = min(float(record["minimum_pdf"]) for record in records)
    if maximum_speed >= 0.05:
        raise FlowError(STEADY_FAILED, f"Maximum lattice speed {maximum_speed} >= 0.05")
    if minimum_pdf <= 0.0:
        raise FlowError(STEADY_FAILED, f"Minimum PDF {minimum_pdf} <= 0")
    final = records[-1]
    return {
        "status": "PASS",
        "record_count": len(records),
        "first_iteration": int(records[0]["iteration"]),
        "final_iteration": int(final["iteration"]),
        "globbc_counts": counts,
        "maximum_relative_error": max(float(item["relative_error"]) for item in records),
        "rho_boundary_range": [
            min(float(item["rho_boundary"]) for item in records),
            max(float(item["rho_boundary"]) for item in records),
        ],
        "pressure_range_pa": [
            min(float(item["pressure_pa"]) for item in records),
            max(float(item["pressure_pa"]) for item in records),
        ],
        "maximum_lattice_speed": maximum_speed,
        "minimum_pdf": minimum_pdf,
        "final_record": final,
        "all_finite": True,
    }


def _live_safety_failure(text: str) -> str | None:
    for line in text.splitlines():
        if "ADAPTIVE_FLUX_PRESSURE" not in line:
            continue
        if re.search(r"(?i)(?<![A-Za-z])(?:nan|[-+]?inf(?:inity)?)(?![A-Za-z])", line):
            return "Adaptive controller emitted NaN/Inf"
        match = CONTROLLER_PATTERN.search(line)
        if not match:
            continue
        record = _match_record(match)
        if float(record["max_lattice_velocity"]) >= 0.05:
            return f"Adaptive controller speed gate failed at iteration {record['iteration']}"
        if float(record["minimum_pdf"]) <= 0.0:
            return f"Adaptive controller PDF gate failed at iteration {record['iteration']}"
    return None


class CheckpointTracker:
    """Record only restart pairs that are complete and stable on WSL ext4."""

    def __init__(
        self,
        *,
        runtime_root: Path,
        runtime_root_wsl: str,
        manifest_path: Path,
    ) -> None:
        self.runtime_root = Path(runtime_root)
        self.runtime_root_wsl = runtime_root_wsl.rstrip("/")
        self.manifest_path = Path(manifest_path)
        self.records: list[dict[str, Any]] = []
        self._recorded_headers: set[str] = set()
        self._signatures: dict[str, tuple[int, int]] = {}
        self._write()

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

    def scan(self, *, force_complete: bool = False) -> None:
        restart_root = self.runtime_root / "restart"
        for header_path in sorted(
            restart_root.glob(f"{SHORT_SIMULATION_NAME}_header_*.lua")
        ):
            header_key = str(header_path)
            if header_key in self._recorded_headers:
                continue
            try:
                header = parse_restart_header(header_path)
            except FlowError:
                # A header name can become visible before its final bytes are flushed.
                continue
            binary_path = header.binary_path
            if (
                header.n_elems != EXPECTED_FLUID_CELLS
                or header.n_components != 19
                or header.n_dofs != 1
                or not binary_path.is_file()
            ):
                continue
            stat = binary_path.stat()
            if stat.st_size != EXPECTED_RESTART_BYTES:
                self._signatures[header_key] = (stat.st_size, stat.st_mtime_ns)
                continue
            signature = (stat.st_size, stat.st_mtime_ns)
            previously_seen = self._signatures.get(header_key)
            self._signatures[header_key] = signature
            if not force_complete and previously_seen != signature:
                continue
            header_relative = header_path.relative_to(self.runtime_root).as_posix()
            binary_relative = binary_path.relative_to(self.runtime_root).as_posix()
            self.records.append(
                {
                    "iteration": header.iteration,
                    "header_path_wsl": f"{self.runtime_root_wsl}/{header_relative}",
                    "binary_path_wsl": f"{self.runtime_root_wsl}/{binary_relative}",
                    "header_path_windows": str(header_path),
                    "binary_path_windows": str(binary_path),
                    "file_size_bytes": stat.st_size,
                    "sha256": sha256_file(binary_path),
                    "complete": True,
                }
            )
            self.records.sort(key=lambda item: int(item["iteration"]))
            self._recorded_headers.add(header_key)
            self._write()


WSL_INFRASTRUCTURE_PATTERN = re.compile(
    r"(?i)(getpwnam|getpwuid|UtilInitGroups|CreateInstance|E_FAIL|"
    r"Input/output error|I/O error|OSError\s*22|WSL instance stopped)"
)


def _run_monitored_wsl(
    *,
    distribution: str,
    workdir_wsl: str,
    command: list[str],
    stdout_path: Path,
    stderr_path: Path,
    controller_csv_path: Path,
    timeout_s: int,
    checkpoint_tracker: CheckpointTracker | None = None,
    minimum_controller_iteration_exclusive: int | None = None,
) -> MonitoredRun:
    """Stream WSL stdout/stderr directly into durable Windows evidence."""

    shell = f"cd {shlex.quote(workdir_wsl)} && exec {shlex.join(command)}"
    started = time.perf_counter()
    state: dict[str, Any] = {
        "safety_failure": None,
        "infrastructure_failure": None,
        "continuity_failure": None,
        "controller_record_count": 0,
        "latest_controller_iteration": None,
    }
    state_lock = threading.Lock()
    stop_event = threading.Event()
    wrapper_timeout = False
    next_progress = started + 60.0
    next_checkpoint_scan = started + 10.0
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    controller_csv_path.parent.mkdir(parents=True, exist_ok=True)
    with (
        stdout_path.open("w", encoding="utf-8", errors="replace", buffering=1) as stdout_sink,
        stderr_path.open("w", encoding="utf-8", errors="replace", buffering=1) as stderr_sink,
        controller_csv_path.open("w", encoding="utf-8", newline="", buffering=1) as controller_sink,
    ):
        writer = csv.DictWriter(controller_sink, fieldnames=CONTROLLER_FIELDS)
        writer.writeheader()
        controller_sink.flush()
        process = subprocess.Popen(
            ["wsl.exe", "-d", distribution, "--", "bash", "-lc", shell],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )

        def consume_stdout() -> None:
            assert process.stdout is not None
            try:
                for line in process.stdout:
                    stdout_sink.write(line)
                    stdout_sink.flush()
                    if WSL_INFRASTRUCTURE_PATTERN.search(line):
                        with state_lock:
                            if state["infrastructure_failure"] is None:
                                state["infrastructure_failure"] = line.strip()
                        stop_event.set()
                    failure = _live_safety_failure(line)
                    match = CONTROLLER_PATTERN.search(line)
                    with state_lock:
                        if failure and state["safety_failure"] is None:
                            state["safety_failure"] = failure
                            stop_event.set()
                        if match:
                            record = _match_record(match)
                            writer.writerow(record)
                            controller_sink.flush()
                            if (
                                state["controller_record_count"] == 0
                                and minimum_controller_iteration_exclusive is not None
                                and int(record["iteration"])
                                <= minimum_controller_iteration_exclusive
                            ):
                                state["continuity_failure"] = (
                                    "First controller iteration "
                                    f"{record['iteration']} <= resume iteration "
                                    f"{minimum_controller_iteration_exclusive}"
                                )
                                stop_event.set()
                            state["controller_record_count"] += 1
                            state["latest_controller_iteration"] = int(
                                record["iteration"]
                            )
            except (OSError, ValueError) as error:
                with state_lock:
                    state["infrastructure_failure"] = f"stdout PIPE failure: {error}"
                stop_event.set()

        def consume_stderr() -> None:
            assert process.stderr is not None
            try:
                for line in process.stderr:
                    stderr_sink.write(line)
                    stderr_sink.flush()
                    if WSL_INFRASTRUCTURE_PATTERN.search(line):
                        with state_lock:
                            if state["infrastructure_failure"] is None:
                                state["infrastructure_failure"] = line.strip()
                        stop_event.set()
            except (OSError, ValueError) as error:
                with state_lock:
                    state["infrastructure_failure"] = f"stderr PIPE failure: {error}"
                stop_event.set()

        stdout_thread = threading.Thread(target=consume_stdout, daemon=True)
        stderr_thread = threading.Thread(target=consume_stderr, daemon=True)
        stdout_thread.start()
        stderr_thread.start()
        termination_sent = False
        while process.poll() is None:
            time.sleep(0.1)
            elapsed = time.perf_counter() - started
            if checkpoint_tracker is not None and time.perf_counter() >= next_checkpoint_scan:
                try:
                    checkpoint_tracker.scan()
                except (OSError, ValueError, FlowError) as error:
                    with state_lock:
                        state["infrastructure_failure"] = (
                            f"persistent restart evidence failure: {error}"
                        )
                    stop_event.set()
                next_checkpoint_scan += 10.0
            if stop_event.is_set() and not termination_sent:
                try:
                    process.terminate()
                except OSError:
                    pass
                termination_sent = True
            if elapsed >= timeout_s:
                wrapper_timeout = True
                try:
                    process.terminate()
                except OSError:
                    pass
                break
            if time.perf_counter() >= next_progress:
                with state_lock:
                    latest_iteration = state["latest_controller_iteration"]
                    record_count = state["controller_record_count"]
                print(
                    "ADAPTIVE_STEADY_PROGRESS "
                    f"elapsed_s={elapsed:.1f} "
                    f"latest_controller_iter={latest_iteration} "
                    f"controller_records={record_count}",
                    flush=True,
                )
                next_progress += 60.0
        try:
            returncode = process.wait(timeout=15)
        except subprocess.TimeoutExpired:
            process.kill()
            returncode = process.wait(timeout=15)
        stdout_thread.join(timeout=15)
        stderr_thread.join(timeout=15)
        if stdout_thread.is_alive() or stderr_thread.is_alive():
            with state_lock:
                state["infrastructure_failure"] = "WSL PIPE reader did not terminate"
        if checkpoint_tracker is not None:
            try:
                checkpoint_tracker.scan(force_complete=True)
            except (OSError, ValueError, FlowError) as error:
                with state_lock:
                    if state["infrastructure_failure"] is None:
                        state["infrastructure_failure"] = (
                            f"final persistent restart evidence failure: {error}"
                        )
    return MonitoredRun(
        returncode=returncode,
        wall_time_s=time.perf_counter() - started,
        safety_failure=state["safety_failure"],
        wrapper_timeout=wrapper_timeout,
        infrastructure_failure=state["infrastructure_failure"],
        continuity_failure=state["continuity_failure"],
        controller_record_count=int(state["controller_record_count"]),
        latest_controller_iteration=state["latest_controller_iteration"],
    )


def _run_luac(
    *, distribution: str, workdir_wsl: str, stdout_path: Path, stderr_path: Path
) -> int:
    completed = subprocess.run(
        [
            "wsl.exe",
            "-d",
            distribution,
            "--",
            "bash",
            "-lc",
            f"cd {shlex.quote(workdir_wsl)} && /home/lzy/.local/bin/luac -p musubi.lua",
        ],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
    )
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    return completed.returncode


def _one_wsl_health_check(distribution: str) -> dict[str, Any]:
    script = (
        'id && getent passwd "$(id -u)" && pwd && df -h /home/lzy && '
        f"mkdir -p {RUNTIME_BASE_WSL} && touch {RUNTIME_BASE_WSL}/io_test && "
        f"cat {RUNTIME_BASE_WSL}/io_test && "
        f"findmnt -n -o SOURCE,FSTYPE -T {RUNTIME_BASE_WSL}"
    )
    completed = subprocess.run(
        ["wsl.exe", "-d", distribution, "--", "bash", "-lc", script],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
    )
    combined = f"{completed.stdout}\n{completed.stderr}"
    forbidden = WSL_INFRASTRUCTURE_PATTERN.search(combined)
    filesystem_lines = [
        line.strip()
        for line in completed.stdout.splitlines()
        if "ext4" in line.lower()
    ]
    passed = completed.returncode == 0 and forbidden is None and bool(filesystem_lines)
    return {
        "status": "PASS" if passed else "FAIL",
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "forbidden_error": forbidden.group(0) if forbidden else None,
        "runtime_filesystem": filesystem_lines[-1] if filesystem_lines else None,
        "runtime_root_persistence_verified": passed,
    }


def _wsl_health_preflight(distribution: str) -> dict[str, Any]:
    first = _one_wsl_health_check(distribution)
    result: dict[str, Any] = {
        "status": first["status"],
        "restart_performed": False,
        "attempts": [first],
    }
    if first["status"] == "PASS":
        return result
    shutdown = subprocess.run(
        ["wsl.exe", "--shutdown"],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
    )
    result["restart_performed"] = True
    result["shutdown_returncode"] = shutdown.returncode
    second = _one_wsl_health_check(distribution)
    result["attempts"].append(second)
    result["status"] = second["status"]
    return result


def _run_synthetic_stream_preflight(
    *, distribution: str, workdir_wsl: str, evidence_root: Path
) -> dict[str, Any]:
    evidence_root.mkdir(parents=True, exist_ok=False)
    line = (
        "ADAPTIVE_FLUX_PRESSURE iter={} "
        "target_lattice=2.3478724595760924E-003 "
        "controlled_lattice=2.3478724595760924E-003 "
        "relative_error=1.0000000000000000E-014 "
        "rho_boundary=1.0001000000000000E+000 "
        "pressure_pa=2.3624680000000000E+004 "
        "max_lattice_velocity=1.0000000000000000E-005 "
        "minimum_pdf=2.0000000000000000E-002 "
        f"globBC_count={EXPECTED_INLET_GLOBBC}"
    )
    script = (
        f"line={line!r}\n"
        f"for iteration in range(1, {SYNTHETIC_RECORDS + 1}):\n"
        "    print(line.format(iteration))\n"
    )
    run = _run_monitored_wsl(
        distribution=distribution,
        workdir_wsl=workdir_wsl,
        command=["python3", "-u", "-c", script],
        stdout_path=evidence_root / "synthetic_stdout.log",
        stderr_path=evidence_root / "synthetic_stderr.log",
        controller_csv_path=evidence_root / "synthetic_controller_records.csv",
        timeout_s=120,
    )
    passed = (
        run.returncode == 0
        and run.controller_record_count == SYNTHETIC_RECORDS
        and run.safety_failure is None
        and run.infrastructure_failure is None
        and not run.wrapper_timeout
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "expected_records": SYNTHETIC_RECORDS,
        "received_records": run.controller_record_count,
        "returncode": run.returncode,
        "safety_failure": run.safety_failure,
        "infrastructure_failure": run.infrastructure_failure,
        "wrapper_timeout": run.wrapper_timeout,
        "stdout": str(evidence_root / "synthetic_stdout.log"),
        "stderr": str(evidence_root / "synthetic_stderr.log"),
        "controller_csv": str(evidence_root / "synthetic_controller_records.csv"),
    }


def summarize_controller_csv(path: Path, *, enforce_gates: bool = True) -> dict[str, Any]:
    count = 0
    first: dict[str, float | int] | None = None
    final: dict[str, float | int] | None = None
    globbc_counts: set[int] = set()
    maximum_relative_error = -math.inf
    minimum_rho = math.inf
    maximum_rho = -math.inf
    minimum_pressure = math.inf
    maximum_pressure = -math.inf
    maximum_speed = -math.inf
    minimum_pdf = math.inf
    all_finite = True
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as stream:
        for raw in csv.DictReader(stream):
            record: dict[str, float | int] = {
                "iteration": int(raw["iteration"]),
                "target_lattice": float(raw["target_lattice"]),
                "controlled_lattice": float(raw["controlled_lattice"]),
                "relative_error": float(raw["relative_error"]),
                "rho_boundary": float(raw["rho_boundary"]),
                "pressure_pa": float(raw["pressure_pa"]),
                "max_lattice_velocity": float(raw["max_lattice_velocity"]),
                "minimum_pdf": float(raw["minimum_pdf"]),
                "globBC_count": int(raw["globBC_count"]),
            }
            numeric = [float(record[field]) for field in CONTROLLER_FIELDS[1:-1]]
            all_finite = all_finite and all(math.isfinite(value) for value in numeric)
            count += 1
            first = record if first is None else first
            final = record
            globbc_counts.add(int(record["globBC_count"]))
            maximum_relative_error = max(
                maximum_relative_error, abs(float(record["relative_error"]))
            )
            minimum_rho = min(minimum_rho, float(record["rho_boundary"]))
            maximum_rho = max(maximum_rho, float(record["rho_boundary"]))
            minimum_pressure = min(minimum_pressure, float(record["pressure_pa"]))
            maximum_pressure = max(maximum_pressure, float(record["pressure_pa"]))
            maximum_speed = max(maximum_speed, float(record["max_lattice_velocity"]))
            minimum_pdf = min(minimum_pdf, float(record["minimum_pdf"]))
    if first is None or final is None:
        raise FlowError(STEADY_FAILED, "No adaptive controller records were persisted")
    if enforce_gates:
        if not all_finite:
            raise FlowError(STEADY_FAILED, "Adaptive controller emitted NaN/Inf")
        if sorted(globbc_counts) != [EXPECTED_INLET_GLOBBC]:
            raise FlowError(
                STEADY_FAILED, f"Adaptive inlet globBC changed: {globbc_counts}"
            )
        if maximum_speed >= 0.05:
            raise FlowError(STEADY_FAILED, f"Maximum lattice speed {maximum_speed} >= 0.05")
        if minimum_pdf <= 0.0:
            raise FlowError(STEADY_FAILED, f"Minimum PDF {minimum_pdf} <= 0")
    return {
        "status": "PASS" if all_finite else "FAIL",
        "record_count": count,
        "first_iteration": int(first["iteration"]),
        "final_iteration": int(final["iteration"]),
        "globbc_counts": sorted(globbc_counts),
        "maximum_relative_error": maximum_relative_error,
        "rho_boundary_range": [minimum_rho, maximum_rho],
        "pressure_range_pa": [minimum_pressure, maximum_pressure],
        "maximum_lattice_speed": maximum_speed,
        "minimum_pdf": minimum_pdf,
        "final_record": final,
        "all_finite": all_finite,
    }


def _prepare_short_tracking_root(*, distribution: str, root_wsl: str) -> None:
    """Create the complete short runtime tree before any solver call."""

    completed = subprocess.run(
        [
            "wsl.exe",
            "-d",
            distribution,
            "--",
            "mkdir",
            "-p",
            f"{root_wsl}/tracking/p",
            f"{root_wsl}/tracking/u",
            f"{root_wsl}/restart",
            f"{root_wsl}/timing",
        ],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
    )
    if completed.returncode != 0:
        raise FlowError(STEADY_FAILED, f"Could not create short tracking root: {completed.stderr}")


def _wsl_path_to_unc(distribution: str, path_wsl: str) -> Path:
    tail = path_wsl.strip("/").replace("/", "\\")
    return Path(f"\\\\wsl.localhost\\{distribution}\\{tail}")


def _runtime_path_length(root_wsl: str, path: Path, runtime_root: Path) -> int:
    relative = path.relative_to(runtime_root).as_posix()
    return len(f"{root_wsl}/{relative}")


def _maximum_runtime_output_path_length(root_wsl: str, runtime_root: Path) -> int:
    paths = [runtime_root, *runtime_root.rglob("*")]
    return max(
        len(root_wsl)
        if path == runtime_root
        else _runtime_path_length(root_wsl, path, runtime_root)
        for path in paths
    )


def _predicted_maximum_runtime_path_length(root_wsl: str) -> int:
    candidates = (
        f"{root_wsl}/musubi.lua",
        f"{root_wsl}/tracking/p/{SHORT_SIMULATION_NAME}_p_p00000.res",
        f"{root_wsl}/tracking/u/{SHORT_SIMULATION_NAME}_u_p00000.res",
        f"{root_wsl}/restart/{SHORT_SIMULATION_NAME}_header_1.000E+00.lua",
        f"{root_wsl}/restart/{SHORT_SIMULATION_NAME}_1.000E+00.lsb",
        f"{root_wsl}/timing/timing.res",
    )
    return max(len(path) for path in candidates)


def _archive_runtime(runtime_root: Path, archive_root: Path) -> None:
    """Copy short runtime evidence without removing or altering its source."""

    for item in runtime_root.iterdir():
        destination = archive_root / item.name
        if item.is_dir():
            shutil.copytree(item, destination, dirs_exist_ok=True)
        elif item.is_file():
            shutil.copy2(item, destination)


def run_adaptive_flux_steady(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    head = _git(root, "rev-parse", "HEAD")
    branch = _git(root, "branch", "--show-current")
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    binary_windows = Path(
        r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees\musubi_mcclure_adaptive_flux_20260829_1300\build\musubi_adaptive_flux"
    )
    binary_sha = sha256_file(binary_windows)
    if binary_sha != EXPECTED_BINARY_SHA256:
        raise FlowError(STEADY_FAILED, f"Adaptive binary SHA changed: {binary_sha}")
    if config.mesh.dx_target_m != EXPECTED_DX_M:
        raise FlowError(STEADY_FAILED, f"dx changed: {config.mesh.dx_target_m}")
    if config.solver.maximum_iterations != 1_000_000:
        raise FlowError(STEADY_FAILED, "Current formal maximum_iterations changed")
    if config.solver.wallclock_limit_s != 3600:
        raise FlowError(STEADY_FAILED, "Current formal wallclock_limit_s changed")
    if config.solver.convergence_interval_iterations != CONVERGENCE_INTERVAL:
        raise FlowError(STEADY_FAILED, "Current formal convergence interval changed")
    if config.solver.convergence_nvals != CONVERGENCE_NVALS:
        raise FlowError(STEADY_FAILED, "Current formal convergence nvals changed")

    output_root = config.paths.output_root
    mesh = output_root / AXIS_MESH_RUN / "seeder" / "mesh"
    mesh_qc = parse_mesh_header(mesh)
    if mesh_qc["fluid_element_count"] != EXPECTED_FLUID_CELLS:
        raise FlowError(STEADY_FAILED, "Frozen adaptive mesh fluid count changed")
    inlet_rim = output_root / AXIS_MESH_RUN / "qc" / "inlet_rim_audit.json"
    rim_text = inlet_rim.read_text(encoding="utf-8")
    inlet_match = re.search(r'"total_inlet_d3q19_cells"\s*:\s*(\d+)', rim_text)
    if not inlet_match or int(inlet_match.group(1)) != EXPECTED_INLET_GLOBBC:
        raise FlowError(STEADY_FAILED, "Frozen adaptive inlet globBC evidence changed")

    frozen_paths = (
        root / "cfd_flow.py",
        root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
        *sorted(path for path in mesh.iterdir() if path.is_file()),
        inlet_rim,
        binary_windows,
    )
    frozen_before = _file_manifest(frozen_paths)
    local_frozen_before = _file_manifest(frozen_paths[:-1])
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{STEADY_PREFIX}_{stamp}"
    short_runtime_root_wsl = f"{RUNTIME_BASE_WSL}/a_{stamp[-6:]}"
    short_runtime_root = _wsl_path_to_unc(
        config.apes.wsl_distribution, short_runtime_root_wsl
    )
    restart_dir = run_root / "restart"
    tracking_dir = run_root / "tracking"
    qc_dir = run_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=False)
    tracking_dir.mkdir(parents=True, exist_ok=False)

    summary_path = qc_dir / "adaptive_flux_steady_manifest.json"
    summary: dict[str, Any] = {
        "status": STEADY_FAILED,
        "next": NEXT_REVIEW_STEADY,
        "run_root": str(run_root),
        "actual_head": head,
        "branch": branch,
        "production_pipeline_modified": False,
        "adaptive_binary_sha256": binary_sha,
        "adaptive_binary_rebuilt": False,
        "mesh_path": str(mesh),
        "fluid_cell_count": EXPECTED_FLUID_CELLS,
        "inlet_globbc_count": EXPECTED_INLET_GLOBBC,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "mpi_ranks": MPI_RANKS,
        "harvester_calls": 0,
        "grid_convergence": "NOT_RUN",
        "dx_m": EXPECTED_DX_M,
        "dt_s": EXPECTED_DT_S,
        "target_q_m3_s": TARGET_Q_M3_S,
        "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
        "fresh_initialization": True,
        "maximum_iterations": config.solver.maximum_iterations,
        "wallclock_limit_s": config.solver.wallclock_limit_s,
        "short_runtime_root_wsl": short_runtime_root_wsl,
        "runtime_root_persistence_verified": False,
        "live_wsl_file_polling_removed": True,
        "stdout_pipe_monitoring_used": True,
        "synthetic_stream_expected_records": SYNTHETIC_RECORDS,
        "synthetic_stream_received_records": 0,
        "synthetic_stream_status": "NOT_RUN",
        "checkpoint_interval_iterations": CHECKPOINT_INTERVAL,
        "checkpoint_manifest": str(run_root / "checkpoint_manifest.json"),
        "frozen_files_before": frozen_before,
        "started_at": datetime.now().isoformat(),
    }
    write_json(summary_path, summary)
    checkpoint_tracker: CheckpointTracker | None = None

    try:
        health = _wsl_health_preflight(config.apes.wsl_distribution)
        write_json(qc_dir / "wsl_health_preflight.json", health)
        summary["wsl_health_preflight"] = health["status"]
        summary["wsl_restart_performed"] = health["restart_performed"]
        summary["runtime_root_persistence_verified"] = bool(
            health["attempts"][-1].get("runtime_root_persistence_verified")
        )
        write_json(summary_path, summary)
        if health["status"] != "PASS":
            raise FlowError(WSL_PREFLIGHT_FAILED, NEXT_REPAIR_WSL)
        if short_runtime_root.exists():
            raise FlowError(
                STEADY_FAILED,
                f"Short runtime root already exists: {short_runtime_root_wsl}",
            )
        _prepare_short_tracking_root(
            distribution=config.apes.wsl_distribution,
            root_wsl=short_runtime_root_wsl,
        )
        lua = generate_adaptive_steady_lua(
            mesh_wsl=windows_to_wsl(mesh, config.apes.wsl_distribution),
            restart_wsl="restart",
            pressure_tracking_wsl="tracking/p",
            velocity_tracking_wsl="tracking/u",
            maximum_iterations=config.solver.maximum_iterations,
            wallclock_limit_s=config.solver.wallclock_limit_s,
        )
        lua_path = short_runtime_root / "musubi.lua"
        lua_path.write_text(lua, encoding="utf-8")
        contract = steady_lua_contract(lua)
        write_json(qc_dir / "adaptive_flux_steady_lua_contract.json", contract)
        if contract["status"] != "PASS":
            raise FlowError(STEADY_FAILED, f"Steady Lua contract failed: {contract}")
        luac_returncode = _run_luac(
            distribution=config.apes.wsl_distribution,
            workdir_wsl=short_runtime_root_wsl,
            stdout_path=short_runtime_root / "tracking" / "luac_stdout.log",
            stderr_path=short_runtime_root / "tracking" / "luac_stderr.log",
        )
        if luac_returncode != 0:
            raise FlowError(STEADY_FAILED, "Adaptive steady Lua syntax failed")
        required_mesh = tuple(mesh / name for name in ("header.lua", "elemlist.lsb", "bnd.lua", "bnd.lsb"))
        preflight_checks = {
            "short_runtime_root_exists": short_runtime_root.is_dir(),
            "tracking_exists": (short_runtime_root / "tracking").is_dir(),
            "restart_exists": (short_runtime_root / "restart").is_dir(),
            "timing_exists": (short_runtime_root / "timing").is_dir(),
            "mesh_required_files_exist": all(path.is_file() for path in required_mesh),
            "rendered_musubi_lua_exists": lua_path.is_file(),
            "lua_syntax_pass": luac_returncode == 0,
            "binary_sha_pass": binary_sha == EXPECTED_BINARY_SHA256,
            "cwd_is_persistent_wsl_ext4_root": short_runtime_root_wsl.startswith(
                f"{RUNTIME_BASE_WSL}/"
            ),
            "tracking_paths_are_short_relative": "folder = 'tracking/p/'" in lua
            and "folder = 'tracking/u/'" in lua,
            "restart_path_is_short_relative": "write = 'restart/'" in lua,
            "restart_interval_is_20000": (
                f"interval = {{ iter = {CHECKPOINT_INTERVAL} }}" in lua
            ),
            "timing_path_is_short_relative": "timing_file = 'timing/timing.res'" in lua,
            "only_mesh_uses_project_absolute_path": str(run_root).replace("\\", "/") not in lua,
        }
        preflight = {
            "status": "PASS" if all(preflight_checks.values()) else "FAIL",
            "short_runtime_root_wsl": short_runtime_root_wsl,
            "short_runtime_root_windows": str(short_runtime_root),
            "checks": preflight_checks,
            "predicted_runtime_path_length": _predicted_maximum_runtime_path_length(
                short_runtime_root_wsl
            ),
        }
        preflight_checks["maximum_path_length_at_most_80"] = (
            preflight["predicted_runtime_path_length"] <= 80
        )
        preflight["status"] = "PASS" if all(preflight_checks.values()) else "FAIL"
        write_json(qc_dir / "short_runtime_static_preflight.json", preflight)
        if preflight["status"] != "PASS":
            raise FlowError(STEADY_FAILED, f"Short runtime static preflight failed: {preflight}")

        synthetic = _run_synthetic_stream_preflight(
            distribution=config.apes.wsl_distribution,
            workdir_wsl=short_runtime_root_wsl,
            evidence_root=qc_dir / "streaming_preflight",
        )
        write_json(qc_dir / "streaming_preflight.json", synthetic)
        summary["synthetic_stream_received_records"] = synthetic["received_records"]
        summary["synthetic_stream_status"] = synthetic["status"]
        write_json(summary_path, summary)
        if synthetic["status"] != "PASS":
            raise FlowError(
                STREAM_PREFLIGHT_FAILED,
                "Synthetic WSL stdout PIPE / Windows writer contract failed",
            )
        post_stream_health = _one_wsl_health_check(config.apes.wsl_distribution)
        write_json(qc_dir / "wsl_health_after_streaming_preflight.json", post_stream_health)
        if post_stream_health["status"] != "PASS":
            raise FlowError(
                STREAM_PREFLIGHT_FAILED,
                "WSL health failed after the synthetic streaming preflight",
            )

        checkpoint_tracker = CheckpointTracker(
            runtime_root=short_runtime_root,
            runtime_root_wsl=short_runtime_root_wsl,
            manifest_path=run_root / "checkpoint_manifest.json",
        )

        summary["musubi_calls"] = 1
        summary["status"] = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_STEADY_RUNNING"
        write_json(summary_path, summary)
        run = _run_monitored_wsl(
            distribution=config.apes.wsl_distribution,
            workdir_wsl=short_runtime_root_wsl,
            command=[
                "env",
                "OMPI_ALLOW_RUN_AS_ROOT=1",
                "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
                MPIRUN_WSL,
                "-np",
                str(MPI_RANKS),
                BINARY_WSL,
                "musubi.lua",
            ],
            stdout_path=tracking_dir / "musubi_stdout.log",
            stderr_path=tracking_dir / "musubi_stderr.log",
            controller_csv_path=tracking_dir / "controller_records.csv",
            timeout_s=config.solver.wallclock_limit_s + WRAPPER_MARGIN_S,
            checkpoint_tracker=checkpoint_tracker,
        )
        maximum_runtime_path_length = _maximum_runtime_output_path_length(
            short_runtime_root_wsl, short_runtime_root
        )
        summary.update(
            {
                "musubi_returncode": run.returncode,
                "musubi_wall_time_s": run.wall_time_s,
                "live_safety_failure": run.safety_failure,
                "wrapper_timeout": run.wrapper_timeout,
                "wsl_infrastructure_failure": run.infrastructure_failure,
                "controller_record_count": run.controller_record_count,
                "last_controller_iteration": run.latest_controller_iteration,
                "checkpoint_iterations": [
                    int(item["iteration"]) for item in checkpoint_tracker.records
                ],
                "maximum_runtime_output_path_length": maximum_runtime_path_length,
                "short_runtime_evidence_archived": False,
            }
        )
        write_json(summary_path, summary)
        if maximum_runtime_path_length > 80:
            raise FlowError(
                STEADY_FAILED,
                f"Musubi-visible path length {maximum_runtime_path_length} > 80",
            )
        if run.infrastructure_failure:
            raise FlowError(WSL_RUNTIME_UNSTABLE, run.infrastructure_failure)
        if run.safety_failure:
            raise FlowError(STEADY_FAILED, run.safety_failure)
        if run.wrapper_timeout:
            raise FlowError(STEADY_FAILED, "Adaptive steady wrapper timeout")
        if run.returncode != 0:
            failure_text = (tracking_dir / "musubi_stderr.log").read_text(
                encoding="utf-8", errors="replace"
            )
            if WSL_INFRASTRUCTURE_PATTERN.search(failure_text):
                raise FlowError(WSL_RUNTIME_UNSTABLE, failure_text.strip())
            if "hvs_ascii_module.f90" in failure_text and "End of record" in failure_text:
                raise FlowError(
                    STEADY_FAILED,
                    "Musubi hvs_ascii fixed-record overflow in persistent short runtime",
                )
            raise FlowError(STEADY_FAILED, f"Adaptive Musubi returned {run.returncode}")

        _archive_runtime(short_runtime_root, run_root)
        summary["short_runtime_evidence_archived"] = True
        for checkpoint in checkpoint_tracker.records:
            header_relative = PurePosixPath(checkpoint["header_path_wsl"]).relative_to(
                PurePosixPath(short_runtime_root_wsl)
            )
            binary_relative = PurePosixPath(checkpoint["binary_path_wsl"]).relative_to(
                PurePosixPath(short_runtime_root_wsl)
            )
            archived_header = run_root.joinpath(*header_relative.parts)
            archived_binary = run_root.joinpath(*binary_relative.parts)
            archive_sha = sha256_file(archived_binary)
            checkpoint.update(
                {
                    "archived_header_path": str(archived_header),
                    "archived_binary_path": str(archived_binary),
                    "archived_binary_sha256": archive_sha,
                    "archive_sha256_match": archive_sha == checkpoint["sha256"],
                }
            )
            if not checkpoint["archive_sha256_match"]:
                raise FlowError(
                    STEADY_FAILED,
                    f"Archived checkpoint SHA mismatch at {checkpoint['iteration']}",
                )
        checkpoint_tracker._write()
        stdout_path = tracking_dir / "musubi_stdout.log"
        stderr_path = tracking_dir / "musubi_stderr.log"
        runtime = summarize_controller_csv(tracking_dir / "controller_records.csv")
        write_json(qc_dir / "adaptive_flux_runtime_controller_qc.json", runtime)
        combined_log = stdout_path.read_text(encoding="utf-8", errors="replace") + "\n" + stderr_path.read_text(
            encoding="utf-8", errors="replace"
        )
        steady = parse_official_steady_termination(combined_log)
        write_json(qc_dir / "adaptive_flux_steady_termination_qc.json", steady)
        if not steady["official_steady_termination"]:
            raise FlowError(STEADY_FAILED, "Existing steady criterion was not reached")
        total_iterations = int(steady["confirmation_iteration"])
        if runtime["final_iteration"] != total_iterations:
            raise FlowError(
                STEADY_FAILED,
                f"Final controller iteration {runtime['final_iteration']} != steady iteration {total_iterations}",
            )

        restart_header = restart_dir / f"{SHORT_SIMULATION_NAME}_lastHeader.lua"
        if not restart_header.is_file():
            raise FlowError(STEADY_FAILED, "Final restart lastHeader is missing")
        header = parse_restart_header(restart_header)
        archived_restart_binary = restart_dir / header.binary_path.name
        if not archived_restart_binary.is_file():
            raise FlowError(STEADY_FAILED, "Archived final restart binary is missing")
        if (
            header.iteration != total_iterations
            or header.n_elems != EXPECTED_FLUID_CELLS
            or header.n_components != 19
            or header.n_dofs != 1
        ):
            raise FlowError(STEADY_FAILED, f"Final restart header contract failed: {header}")
        binary_contract = restart_binary_size_contract(
            archived_restart_binary,
            n_elems=header.n_elems,
            n_components=header.n_components,
            n_dofs=header.n_dofs,
        )
        if binary_contract["status"] != "PASS":
            raise FlowError(STEADY_FAILED, "Final restart binary size contract failed")
        write_json(
            qc_dir / "adaptive_flux_final_restart_qc.json",
            {
                "status": "PASS",
                "last_header": str(restart_header),
                "last_header_sha256": sha256_file(restart_header),
                "runtime_binary": str(header.binary_path),
                "binary": str(archived_restart_binary),
                "binary_sha256": sha256_file(archived_restart_binary),
                "iteration": header.iteration,
                "n_elems": header.n_elems,
                "n_components": header.n_components,
                "binary_size_contract": binary_contract,
            },
        )

        frozen_after = _file_manifest(frozen_paths)
        unchanged = frozen_before == frozen_after
        write_json(
            qc_dir / "source_frozen_files_unchanged_qc.json",
            {
                "status": "PASS" if unchanged else "FAIL",
                "source_frozen_files_unchanged": unchanged,
                "before": frozen_before,
                "after": frozen_after,
            },
        )
        if not unchanged:
            raise FlowError(STEADY_FAILED, "Source frozen files changed during steady run")
        summary.update(
            {
                "status": STEADY_PENDING_AUDIT,
                "next": NEXT_EXACT_AUDIT,
                "total_steady_iterations": total_iterations,
                "steady_criterion": "PROJECT_CHARACTERISTIC_SCALE_STEADY_0P1_PERCENT",
                "steady_criterion_status": "PASS",
                "runtime_controller": runtime,
                "final_restart_header": str(restart_header),
                "final_restart_binary": str(archived_restart_binary),
                "runtime_restart_header": str(
                    short_runtime_root / "restart" / restart_header.name
                ),
                "runtime_restart_binary": str(
                    short_runtime_root / "restart" / archived_restart_binary.name
                ),
                "checkpoint_iterations": [
                    int(item["iteration"]) for item in checkpoint_tracker.records
                ],
                "latest_checkpoint_path": (
                    checkpoint_tracker.records[-1]["binary_path_wsl"]
                    if checkpoint_tracker.records
                    else None
                ),
                "short_runtime_root_preserved_through_audit": True,
                "source_frozen_files_unchanged": True,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(summary_path, summary)
        return summary
    except Exception as error:
        error_status = error.status if isinstance(error, FlowError) else STEADY_FAILED
        if error_status == WSL_PREFLIGHT_FAILED:
            next_step = NEXT_REPAIR_WSL
        elif error_status == WSL_RUNTIME_UNSTABLE:
            next_step = NEXT_REVIEW_WSL
        else:
            next_step = NEXT_REVIEW_STEADY
        controller_csv = tracking_dir / "controller_records.csv"
        if controller_csv.is_file():
            try:
                summary["runtime_controller"] = summarize_controller_csv(
                    controller_csv, enforce_gates=False
                )
            except Exception as controller_error:
                summary["runtime_controller_evidence_failure"] = str(controller_error)
        if checkpoint_tracker is not None:
            summary["checkpoint_iterations"] = [
                int(item["iteration"]) for item in checkpoint_tracker.records
            ]
            summary["latest_checkpoint_path"] = (
                checkpoint_tracker.records[-1]["binary_path_wsl"]
                if checkpoint_tracker.records
                else None
            )
        local_frozen_after = _file_manifest(frozen_paths[:-1])
        local_unchanged = local_frozen_before == local_frozen_after
        frozen_after: dict[str, Any] | None = None
        full_unchanged = False
        if error_status != WSL_RUNTIME_UNSTABLE:
            try:
                frozen_after = _file_manifest(frozen_paths)
                full_unchanged = frozen_before == frozen_after
            except (OSError, ValueError) as manifest_error:
                summary["full_frozen_manifest_failure"] = str(manifest_error)
        summary.update(
            {
                "status": error_status,
                "next": next_step,
                "first_failure": str(error),
                "source_frozen_files_unchanged": (
                    full_unchanged if frozen_after is not None else local_unchanged
                ),
                "source_frozen_local_files_unchanged": local_unchanged,
                "adaptive_binary_post_run_sha_verified": frozen_after is not None,
                "frozen_files_after": frozen_after,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(summary_path, summary)
        return summary
