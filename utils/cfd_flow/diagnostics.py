"""Low-cost read-only diagnostics for the frozen Musubi restart."""

from __future__ import annotations

import csv
import json
import math
import re
import shutil
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from .apes import (
    compute_lattice_scaling,
    inspect_apes_environment,
    load_boundary_conditions,
    run_wsl_tool,
    windows_to_wsl,
)
from .config import load_cfd_flow_config
from .geometry import load_frozen_surface_partition
from .io import FlowError, load_flow_inputs, sha256_file, write_json
from .qc import numerical_port_fluxes, validate_and_convert_flow_vtu


DIAGNOSTIC_STATUS = "CFD_FLOW_DIAGNOSTIC_COMPLETE"
PROBE_STATUS = "CFD_FLOW_CONVERGENCE_PROBE_COMPLETE"
RESTART_INVALID = "CFD_FLOW_RESTART_INVALID"
HARVEST_FAILED = "CFD_FLOW_DIAGNOSTIC_HARVEST_FAILED"
SOURCE_RUN_NAME = "musubi_solver_recovery_anchor003274_20260828_164912"
FROZEN_SEEDER_NAME = "musubi_recovery_anchor003274_20260828_162530"
START_ITERATION = 162_464
ADDITIONAL_ITERATIONS = 5_000
END_ITERATION = START_ITERATION + ADDITIONAL_ITERATIONS
DIAGNOSTIC_WALLCLOCK_S = 300
WRAPPER_HARD_TIMEOUT_S = 330
ASCII_TMP_ROOT_WSL = "/tmp/u3d"
ASCII_PRESSURE_FOLDER_WSL = "/tmp/u3d/p"
ASCII_VELOCITY_FOLDER_WSL = "/tmp/u3d/u"
ASCII_PRESSURE_LABEL = "p"
ASCII_VELOCITY_LABEL = "u"
PINNED_TREELM_LABEL_LEN = 80
ASCII_PROJECT_LENGTH_GATE = 60
OLD_ASCII_FILENAME_LENGTH = 193
PROBE_REVISION = "SHORT_ASCII_PATH_V2"
OFFICIAL_EQUIVALENT_STATUS = "CFD_FLOW_OFFICIAL_EQUIVALENT_CONVERGENCE_PROBE_COMPLETE"
OFFICIAL_EQUIVALENT_REVISION = "OFFLINE_OFFICIAL_EQUIVALENT_NVALS100_V1"
PREVIOUS_PROBE_NAME = "musubi_convergence_probe_anchor003274_20260828_210236"
OFFICIAL_EQUIVALENT_START_ITERATION = 167_464
OFFICIAL_EQUIVALENT_ADDITIONAL_ITERATIONS = 6_000
OFFICIAL_EQUIVALENT_END_ITERATION = (
    OFFICIAL_EQUIVALENT_START_ITERATION + OFFICIAL_EQUIVALENT_ADDITIONAL_ITERATIONS
)
OFFICIAL_EQUIVALENT_SOLVER_BUDGET_S = 180
OFFICIAL_EQUIVALENT_WRAPPER_TIMEOUT_S = 240
OFFICIAL_EQUIVALENT_NVALS = 100
PROJECT_STEADY_CONFIRMED_STATUS = "CFD_FLOW_PROJECT_STEADY_0P1PCT_CONFIRMED"
PROJECT_STEADY_NOT_REACHED_STATUS = "CFD_FLOW_PROJECT_STEADY_CRITERION_NOT_REACHED"
PROJECT_STEADY_REVISION = "PROJECT_CHARACTERISTIC_SCALE_STEADY_0P1_PERCENT_V1"
PROJECT_STEADY_POLICY = "PROJECT_CHARACTERISTIC_SCALE_STEADY_0P1_PERCENT"
PROJECT_STEADY_SOURCE_RUN_NAME = (
    "musubi_official_equivalent_probe_anchor003274_20260828_223213"
)
PROJECT_STEADY_START_ITERATION = 173_464
PROJECT_STEADY_ADDITIONAL_ITERATIONS = 30_000
PROJECT_STEADY_END_ITERATION = (
    PROJECT_STEADY_START_ITERATION + PROJECT_STEADY_ADDITIONAL_ITERATIONS
)
PROJECT_STEADY_RELATIVE_FRACTION = 1.0e-3
PROJECT_STEADY_EXPECTED_BUDGET_S = 750
PROJECT_STEADY_WRAPPER_TIMEOUT_S = 900
PRESSURE_THRESHOLD_PA = 1.0e-3
VELOCITY_THRESHOLD_M_S = 1.0e-9
PRESSURE_REFERENCE_PA = 23622.32012800001
TARGET_Q_M3_S = 7.693508475538942e-16


def _float_pattern() -> str:
    return r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def parse_total_density_log(text: str) -> list[dict[str, float | int]]:
    """Extract the compact runtime blocks written by Musubi."""

    number = _float_pattern()
    pattern = re.compile(
        rf"iterations:\s*(\d+).*?"
        rf"simTime\s*:\s*({number}).*?"
        rf"wallClock\s*:\s*({number}).*?"
        rf"field\s*\|\s*total density\s*\n\s*1\s*\|\s*({number})",
        re.DOTALL,
    )
    by_iteration: dict[int, dict[str, float | int]] = {}
    for match in pattern.finditer(text):
        iteration = int(match.group(1))
        by_iteration[iteration] = {
            "iteration": iteration,
            "physical_time_s": float(match.group(2)),
            "wall_clock_s": float(match.group(3)),
            "total_density": float(match.group(4)),
        }
    records = [by_iteration[key] for key in sorted(by_iteration)]
    if len(records) < 2:
        raise FlowError("CFD_FLOW_DIAGNOSTIC_INPUT_INVALID", "No total-density series in Musubi log")
    return records


def write_total_density_csv(path: Path, records: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("iteration", "physical_time_s", "wall_clock_s", "total_density"),
        )
        writer.writeheader()
        writer.writerows(records)


def summarize_total_density(records: list[dict[str, float | int]]) -> dict[str, Any]:
    iterations = np.asarray([item["iteration"] for item in records], dtype=float)
    density = np.asarray([item["total_density"] for item in records], dtype=float)
    final_iteration = int(iterations[-1])
    final_density = float(density[-1])

    landmarks: dict[str, dict[str, float | int]] = {}
    for target in (0, 40_000, 80_000, 120_000, 150_000, START_ITERATION):
        index = int(np.argmin(np.abs(iterations - target)))
        landmarks[str(target)] = {
            "requested_iteration": target,
            "actual_iteration": int(iterations[index]),
            "total_density": float(density[index]),
        }

    windows: dict[str, dict[str, float]] = {}
    for width in (100, 1_000, 5_000, 10_000):
        start_iteration = final_iteration - width
        start_density = float(np.interp(start_iteration, iterations, density))
        delta = final_density - start_density
        windows[str(width)] = {
            "start_density_interpolated": start_density,
            "final_density": final_density,
            "delta": delta,
            "absolute_delta": abs(delta),
            "relative_delta": delta / start_density,
            "per_1000_iteration_slope": delta * 1000.0 / width,
        }
    return {
        "status": "PASS",
        "interpretation": "DIAGNOSTIC_ONLY_NOT_A_STEADY_HARD_GATE",
        "sample_count": len(records),
        "first_iteration": int(iterations[0]),
        "final_iteration": final_iteration,
        "landmarks": landmarks,
        "windows": windows,
    }


def _extract_lua_number(text: str, name: str) -> float:
    match = re.search(rf"(?m)^\s*{re.escape(name)}\s*=\s*({_float_pattern()})\s*,?", text)
    if match is None:
        raise FlowError(RESTART_INVALID, f"Missing {name} in restart header")
    return float(match.group(1))


def _directory_manifest(directory: Path) -> dict[str, dict[str, Any]]:
    return {
        path.relative_to(directory).as_posix(): {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def validate_restart(
    source_run: Path,
    frozen_seeder_run: Path,
) -> dict[str, Any]:
    """Validate the terminal restart without changing any historical output."""

    musubi_dir = source_run / "musubi"
    restart_dir = musubi_dir / "restart"
    last_header = restart_dir / "roi003274_steady_lbm_lastHeader.lua"
    timestamp_header = restart_dir / "roi003274_steady_lbm_header_3.966E-03.lua"
    for path in (last_header, timestamp_header):
        if not path.is_file():
            raise FlowError(RESTART_INVALID, f"Missing restart header: {path}")
    text = last_header.read_text(encoding="utf-8", errors="strict")
    binary_match = re.search(r"['\"]([^'\"]+\.lsb)['\"]", text)
    mesh_match = re.search(r"(?m)^\s*mesh\s*=\s*['\"]([^'\"]+)['\"]", text)
    if binary_match is None or mesh_match is None:
        raise FlowError(RESTART_INVALID, "Restart header lacks binary_name or mesh")
    binary = (musubi_dir / binary_match.group(1)).resolve()
    mesh = (musubi_dir / mesh_match.group(1)).resolve()
    expected_mesh = (frozen_seeder_run / "seeder" / "mesh").resolve()
    iteration = int(_extract_lua_number(text, "iter"))
    sim_time = _extract_lua_number(text, "sim")
    elements = int(_extract_lua_number(text, "nElems"))
    checks = {
        "header_references_existing_lsb": binary.is_file(),
        "timestamp_header_matches_last_header": (
            sha256_file(timestamp_header) == sha256_file(last_header)
        ),
        "n_elems_exact": elements == 221_109,
        "iteration_exact": iteration == START_ITERATION,
        "physical_time_exact": math.isclose(
            sim_time,
            0.00396640625,
            rel_tol=0.0,
            abs_tol=1.0e-14,
        ),
        "mesh_path_is_frozen_seeder": mesh == expected_mesh and mesh.is_dir(),
        "solver_config_exists": (musubi_dir / "musubi.lua").is_file(),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "iteration": iteration,
        "physical_time_s": sim_time,
        "n_elems": elements,
        "binary": str(binary),
        "mesh": str(mesh),
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in (binary, timestamp_header, last_header)
            if path.is_file()
        },
    }
    if result["status"] != "PASS":
        raise FlowError(RESTART_INVALID, f"Restart contract failed: {checks}")
    return result


def validate_continuation_restart(
    restart_header: Path,
    *,
    solver_workdir: Path,
    frozen_seeder_run: Path,
    expected_iteration: int,
) -> dict[str, Any]:
    """Validate a diagnostic restart while preserving its absolute WSL binary reference."""

    header = Path(restart_header).resolve()
    if not header.is_file():
        raise FlowError(RESTART_INVALID, f"Missing restart header: {header}")
    text = header.read_text(encoding="utf-8", errors="strict")
    binary_match = re.search(r"['\"]([^'\"]+\.lsb)['\"]", text)
    mesh_match = re.search(r"(?m)^\s*mesh\s*=\s*['\"]([^'\"]+)['\"]", text)
    if binary_match is None or mesh_match is None:
        raise FlowError(RESTART_INVALID, "Continuation restart lacks binary_name or mesh")
    binary = header.parent / Path(binary_match.group(1)).name
    mesh = (Path(solver_workdir).resolve() / mesh_match.group(1)).resolve()
    expected_mesh = (Path(frozen_seeder_run).resolve() / "seeder" / "mesh").resolve()
    timestamp_headers = sorted(header.parent.glob("*_header_*.lua"))
    iteration = int(_extract_lua_number(text, "iter"))
    sim_time = _extract_lua_number(text, "sim")
    elements = int(_extract_lua_number(text, "nElems"))
    checks = {
        "header_references_existing_lsb": binary.is_file(),
        "timestamp_header_matches_last_header": (
            len(timestamp_headers) == 1
            and sha256_file(timestamp_headers[0]) == sha256_file(header)
        ),
        "n_elems_exact": elements == 221_109,
        "iteration_exact": iteration == expected_iteration,
        "physical_time_matches_iteration": math.isclose(
            sim_time,
            expected_iteration * 2.44140625e-8,
            rel_tol=0.0,
            abs_tol=1.0e-14,
        ),
        "mesh_path_is_frozen_seeder": mesh == expected_mesh and mesh.is_dir(),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "iteration": iteration,
        "physical_time_s": sim_time,
        "n_elems": elements,
        "last_header": str(header),
        "last_header_sha256": sha256_file(header),
        "binary": str(binary),
        "binary_sha256": sha256_file(binary) if binary.is_file() else None,
        "mesh": str(mesh),
    }
    if result["status"] != "PASS":
        raise FlowError(RESTART_INVALID, f"Continuation restart contract failed: {checks}")
    return result


def generate_diagnostic_harvester_lua(output_folder_wsl: str) -> str:
    return f"""-- Diagnostic harvest of the non-converged iteration-162464 restart.
require 'musubi'
restart = {{ read = 'restart/roi003274_steady_lbm_lastHeader.lua' }}
tracking = {{
  label = 'current_snapshot_162464',
  folder = '{output_folder_wsl.rstrip('/')}/',
  variable = {{ 'pressure_phy', 'velocity_phy', 'vel_mag_phy' }},
  shape = {{ kind = 'all' }},
  output = {{ format = 'vtk' }}
}}
"""


def _replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
            f"Expected one occurrence of {old!r}, found {text.count(old)}",
        )
    return text.replace(old, new, 1)


def generate_diagnostic_musubi_lua(
    production_lua: str,
    *,
    pressure_folder_wsl: str,
    velocity_folder_wsl: str,
    restart_folder_wsl: str,
    timing_file_wsl: str,
    start_iteration: int = START_ITERATION,
    end_iteration: int = END_ITERATION,
    restart_read: str = "restart/roi003274_steady_lbm_lastHeader.lua",
    pressure_convergence_threshold_pa: float = PRESSURE_THRESHOLD_PA,
    velocity_convergence_threshold_m_s: float = VELOCITY_THRESHOLD_M_S,
) -> str:
    """Change only restart, stop ceiling, tracking, and output destinations."""

    text = production_lua
    text = _replace_once(
        text,
        "-- Generated production configuration; official Musubi syntax.",
        "-- Diagnostic restart continuation; frozen production physics and BC.",
    )
    text = _replace_once(text, "timing_file = 'mus_timing.res'", f"timing_file = '{timing_file_wsl}'")
    text = _replace_once(text, "maximum_iterations = 1000000", f"maximum_iterations = {end_iteration}")
    text = _replace_once(
        text,
        "{ threshold = 0.001, operator = '<=' }",
        f"{{ threshold = {pressure_convergence_threshold_pa:.17g}, operator = '<=' }}",
    )
    text = _replace_once(
        text,
        "{ threshold = 1.0000000000000001e-09, operator = '<=' }",
        f"{{ threshold = {velocity_convergence_threshold_m_s:.17g}, operator = '<=' }}",
    )
    text = _replace_once(
        text,
        "max = { iter = maximum_iterations, clock = 3600 }",
        "max = { iter = maximum_iterations }",
    )
    text = _replace_once(
        text,
        "time_control = { min = { iter = 0 }, max = { iter = maximum_iterations }, interval = { iter = 100 } },",
        f"time_control = {{ min = {{ iter = {start_iteration} }}, "
        "max = { iter = maximum_iterations }, interval = { iter = 100 } },",
    )
    pressure_folder = pressure_folder_wsl.rstrip("/") + "/"
    velocity_folder = velocity_folder_wsl.rstrip("/") + "/"
    restart_folder = restart_folder_wsl.rstrip("/") + "/"
    diagnostic_output = f"""tracking = {{
  {{
    label = '{ASCII_PRESSURE_LABEL}',
    folder = '{pressure_folder}',
    variable = {{ 'pressure_phy' }},
    shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {start_iteration} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = 100 }} }},
    output = {{ format = 'ascii' }}
  }},
  {{
    label = '{ASCII_VELOCITY_LABEL}',
    folder = '{velocity_folder}',
    variable = {{ 'vel_mag_phy' }},
    shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {start_iteration} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = 100 }} }},
    output = {{ format = 'ascii' }}
  }}
}}

restart = {{
  read = '{restart_read}',
  write = '{restart_folder}'
}}
"""
    text = _replace_once(text, "restart = { write = 'restart/' }\n", diagnostic_output)
    return text


def ascii_path_preflight(simulation_name: str = "roi003274_steady_lbm") -> dict[str, Any]:
    """Model TreElm's full ASCII result path, including its rank suffix."""

    candidates: list[str] = []
    required_examples: list[str] = []
    for folder, label in (
        (ASCII_PRESSURE_FOLDER_WSL, ASCII_PRESSURE_LABEL),
        (ASCII_VELOCITY_FOLDER_WSL, ASCII_VELOCITY_LABEL),
    ):
        for rank in (0, 7):
            required_examples.append(f"{folder}/{label}_p{rank:05d}.res")
            candidates.append(f"{folder}/{simulation_name}_{label}_p{rank:05d}.res")
    lengths = {path: len(path) for path in (*required_examples, *candidates)}
    maximum = max(lengths.values())
    checks = {
        "pinned_label_len_exact": PINNED_TREELM_LABEL_LEN == 80,
        "all_below_pinned_label_len": maximum < PINNED_TREELM_LABEL_LEN,
        "project_safety_margin": maximum <= ASCII_PROJECT_LENGTH_GATE,
        "required_rank_0_and_7_examples_checked": len(required_examples) == 4,
        "full_path_and_rank_suffix_checked": all(path.endswith(("p00000.res", "p00007.res")) for path in lengths),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "pinned_label_len": PINNED_TREELM_LABEL_LEN,
        "old_predicted_filename_length": OLD_ASCII_FILENAME_LENGTH,
        "max_ascii_filename_length": maximum,
        "predicted_filename_lengths": lengths,
    }


def diagnostic_lua_contract(
    production_lua: str,
    diagnostic_lua: str,
    *,
    start_iteration: int = START_ITERATION,
    end_iteration: int = END_ITERATION,
    restart_read: str = "restart/roi003274_steady_lbm_lastHeader.lua",
    pressure_convergence_threshold_pa: float = PRESSURE_THRESHOLD_PA,
    velocity_convergence_threshold_m_s: float = VELOCITY_THRESHOLD_M_S,
) -> dict[str, Any]:
    frozen_fragments = (
        "mesh = '../../musubi_recovery_anchor003274_20260828_162530/seeder/mesh/'",
        "dx = 1.9999999999999999e-07",
        "dt = 2.4414062499999991e-08",
        "rho0_phy = 1056",
        "nu_phy = 3.27e-06",
        "bulk_viscosity_phy = (2.0 / 3.0) * nu_phy",
        "kind = 'fluid', layout = 'd3q19', relaxation = 'bgk'",
        "label = 'wall', kind = 'wall_libb'",
        "label = 'inlet', kind = 'mfr_eq', mass_flowrate = 8.124344950169123e-13",
        "label = 'outlet_01', kind = 'pressure_eq'",
        "label = 'outlet_02', kind = 'pressure_eq'",
        "label = 'outlet_03', kind = 'pressure_eq'",
        "function outlet_01_pressure(x, y, z, t) return 23636.865106101286 end",
        "function outlet_02_pressure(x, y, z, t) return 23754.524677223188 end",
        "function outlet_03_pressure(x, y, z, t) return 23608.619501326699 end",
    )
    checks = {
        "production_fragments_present": all(item in production_lua for item in frozen_fragments),
        "frozen_physics_bc_preserved": all(item in diagnostic_lua for item in frozen_fragments),
        "restart_read_present": f"read = '{restart_read}'" in diagnostic_lua,
        "iteration_ceiling_exact": f"maximum_iterations = {end_iteration}" in diagnostic_lua,
        "tracking_start_exact": diagnostic_lua.count(f"min = {{ iter = {start_iteration} }}") >= 2,
        "internal_wallclock_ceiling_absent": "clock =" not in diagnostic_lua,
        "pinned_convergence_norm_average": "norm = 'average'" in diagnostic_lua,
        "pinned_convergence_nvals_100": "nvals = 100" in diagnostic_lua,
        "pinned_convergence_absolute": "absolute = true" in diagnostic_lua,
        "pinned_convergence_spatial_average": (
            "reduction = { 'average', 'average' }" in diagnostic_lua
        ),
        "pressure_convergence_threshold_exact": (
            f"{{ threshold = {pressure_convergence_threshold_pa:.17g}, operator = '<=' }}"
            in diagnostic_lua
        ),
        "velocity_convergence_threshold_exact": (
            f"{{ threshold = {velocity_convergence_threshold_m_s:.17g}, operator = '<=' }}"
            in diagnostic_lua
        ),
        "pressure_tracking_ascii": (
            f"label = '{ASCII_PRESSURE_LABEL}'" in diagnostic_lua
            and "variable = { 'pressure_phy' }" in diagnostic_lua
        ),
        "velocity_tracking_ascii": (
            f"label = '{ASCII_VELOCITY_LABEL}'" in diagnostic_lua
            and "variable = { 'vel_mag_phy' }" in diagnostic_lua
        ),
        "tracking_interval_100": diagnostic_lua.count("interval = { iter = 100 }") >= 2,
        "no_tracking_vtk": "output = { format = 'vtk' }" not in diagnostic_lua,
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _find_harvested_vtk(directory: Path) -> Path:
    parallel = sorted(directory.rglob("*.pvtu"))
    serial = sorted(directory.rglob("*.vtu"))
    candidates = parallel or serial
    if not candidates:
        raise FlowError(HARVEST_FAILED, "mus_harvesting produced no VTU/PVTU")
    return candidates[-1]


def _mass_classification(error: float) -> str:
    if error <= 0.01:
        return "ALREADY_GOOD"
    if error <= 0.05:
        return "CLOSE"
    if error <= 0.20:
        return "STILL_DEVELOPING"
    return "FAR_FROM_STEADY_OR_BC_REVIEW"


def current_snapshot_qc(
    raw_vtk: Path,
    output_vtu: Path,
    *,
    partition: Any,
    scaling: Any,
    target_q_m3_s: float,
) -> tuple[pv.UnstructuredGrid, dict[str, Any], dict[str, Any]]:
    grid, field_qc = validate_and_convert_flow_vtu(
        raw_vtk,
        output_vtu,
        pressure_reference_pa=PRESSURE_REFERENCE_PA,
    )
    velocity_max = field_qc["velocity_m_s"]["max"]
    mach = velocity_max * scaling.dt_s / scaling.dx_m / math.sqrt(1.0 / 3.0)
    field_qc.update(
        {
            "snapshot_status": "NON_CONVERGED_DIAGNOSTIC_SNAPSHOT",
            "actual_lattice_mach": mach,
            "actual_lattice_mach_below_0_05": mach < 0.05,
            "proteus_final_compatibility_claimed": False,
        }
    )
    fluxes, measured_pressure = numerical_port_fluxes(grid, partition, scaling.dx_m)
    q_in = float(fluxes["inlet"])
    q_out = tuple(float(fluxes[f"outlet_{index:02d}"]) for index in range(1, 4))
    denominator = abs(q_in)
    mass_error = math.inf if denominator == 0.0 else abs(abs(q_in) - sum(abs(q) for q in q_out)) / denominator
    signed_outward = {"inlet": -q_in, **{f"outlet_{i:02d}": q for i, q in enumerate(q_out, 1)}}
    inlet_error = math.inf if target_q_m3_s == 0.0 else abs(abs(q_in) - target_q_m3_s) / target_q_m3_s
    flow_qc = {
        "status": "PASS",
        "snapshot_status": "NON_CONVERGED_DIAGNOSTIC_SNAPSHOT",
        "q_in_m3_s": q_in,
        "q_out_01_m3_s": q_out[0],
        "q_out_02_m3_s": q_out[1],
        "q_out_03_m3_s": q_out[2],
        "signed_outward_flux_m3_s": signed_outward,
        "flow_direction": {
            "inlet_inward": q_in > 0.0,
            "outlet_01_outward": q_out[0] > 0.0,
            "outlet_02_outward": q_out[1] > 0.0,
            "outlet_03_outward": q_out[2] > 0.0,
        },
        "flow_direction_pass": q_in > 0.0 and all(value > 0.0 for value in q_out),
        "mass_conservation_relative_error": mass_error,
        "mass_conservation_classification": _mass_classification(mass_error),
        "mass_error_is_diagnostic_not_hard_fail": True,
        "inlet_target_m3_s": target_q_m3_s,
        "inlet_relative_error": inlet_error,
        "outlet_measured_gauge_pressure_pa": measured_pressure,
        "outlet_target_gauge_pressure_pa": {
            "outlet_01": 14.544978101274268,
            "outlet_02": 132.20454922317552,
            "outlet_03": -13.700626673311461,
        },
    }
    return grid, field_qc, flow_qc


def create_velocity_snapshot(grid: pv.UnstructuredGrid, path: Path) -> None:
    centers = np.asarray(grid.cell_centers().points, dtype=float) * 1.0e6
    speed = np.linalg.norm(np.asarray(grid.cell_data["velocity_phy"], dtype=float), axis=1)
    extents = np.ptp(centers, axis=0)
    axes = np.argsort(extents)[-2:]
    labels = ("x", "y", "z")
    figure, axis = plt.subplots(figsize=(8.0, 6.2), constrained_layout=True)
    image = axis.hexbin(
        centers[:, axes[0]],
        centers[:, axes[1]],
        C=speed,
        gridsize=180,
        reduce_C_function=np.max,
        mincnt=1,
        cmap="viridis",
    )
    axis.set_xlabel(f"{labels[axes[0]]} (µm)")
    axis.set_ylabel(f"{labels[axes[1]]} (µm)")
    axis.set_aspect("equal", adjustable="box")
    axis.set_title("Current non-converged velocity snapshot (maximum projection)")
    figure.colorbar(image, ax=axis, label="Velocity magnitude (m/s)")
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _numeric_tracking_rows(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = [float(value) for value in re.findall(_float_pattern(), line)]
        if len(values) >= 2:
            rows.append((values[0], values[-1]))
    return rows


def _load_tracking_series(directory: Path, series_name: str, dt_s: float) -> dict[int, float]:
    files = sorted(directory.glob("*.res"))
    if not files:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_TRACKING_INVALID",
            f"No tracking file for {series_name}",
        )
    values: dict[int, float] = {}
    for path in files:
        for time_value, result in _numeric_tracking_rows(path):
            iteration = int(round(time_value / dt_s)) if abs(time_value) < 1_000.0 else int(round(time_value))
            values[iteration] = result
    if len(values) < 2:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_TRACKING_INVALID",
            f"Too few samples for {series_name}",
        )
    return values


def parse_tracking(directory: Path, dt_s: float) -> list[dict[str, float | int | None]]:
    pressure = _load_tracking_series(directory / ASCII_PRESSURE_LABEL, "pressure", dt_s)
    velocity = _load_tracking_series(directory / ASCII_VELOCITY_LABEL, "velocity", dt_s)
    common = sorted(set(pressure).intersection(velocity))
    if len(common) < 2:
        raise FlowError("CFD_FLOW_DIAGNOSTIC_TRACKING_INVALID", "Pressure/velocity samples do not align")
    records: list[dict[str, float | int | None]] = []
    previous_pressure: float | None = None
    previous_velocity: float | None = None
    for iteration in common:
        pressure_value = pressure[iteration]
        velocity_value = velocity[iteration]
        records.append(
            {
                "iteration": iteration,
                "physical_time_s": iteration * dt_s,
                "avg_pressure_pa": pressure_value,
                "avg_velocity_m_s": velocity_value,
                "delta_pressure_vs_previous": (
                    None if previous_pressure is None else pressure_value - previous_pressure
                ),
                "delta_velocity_vs_previous": (
                    None if previous_velocity is None else velocity_value - previous_velocity
                ),
            }
        )
        previous_pressure = pressure_value
        previous_velocity = velocity_value
    return records


def write_tracking_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fields = (
        "iteration",
        "physical_time_s",
        "avg_pressure_pa",
        "avg_velocity_m_s",
        "delta_pressure_vs_previous",
        "delta_velocity_vs_previous",
    )
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def combine_continuous_tracking(
    previous: list[dict[str, Any]],
    current: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Deduplicate one shared restart endpoint and enforce 100-iteration continuity."""

    if not previous or not current:
        raise FlowError("CFD_FLOW_TRACKING_CONTINUITY_FAILED", "Empty tracking series")
    shared_iteration = int(previous[-1]["iteration"])
    if shared_iteration != int(current[0]["iteration"]):
        raise FlowError(
            "CFD_FLOW_TRACKING_CONTINUITY_FAILED",
            "Previous final and current initial iterations do not match",
        )
    pressure_difference = abs(
        float(previous[-1]["avg_pressure_pa"]) - float(current[0]["avg_pressure_pa"])
    )
    velocity_difference = abs(
        float(previous[-1]["avg_velocity_m_s"])
        - float(current[0]["avg_velocity_m_s"])
    )
    pressure_tolerance = max(
        1.0e-9,
        abs(float(previous[-1]["avg_pressure_pa"])) * 1.0e-12,
    )
    velocity_tolerance = max(
        1.0e-15,
        abs(float(previous[-1]["avg_velocity_m_s"])) * 1.0e-9,
    )
    endpoint_checks = {
        "pressure_shared_endpoint_consistent": pressure_difference <= pressure_tolerance,
        "velocity_shared_endpoint_consistent": velocity_difference <= velocity_tolerance,
    }
    if not all(endpoint_checks.values()):
        raise FlowError(
            "CFD_FLOW_TRACKING_CONTINUITY_FAILED",
            f"Shared endpoint values are discontinuous: {endpoint_checks}",
        )

    raw = [*previous, *current[1:]]
    iterations = [int(item["iteration"]) for item in raw]
    interval_checks = [later - earlier for earlier, later in zip(iterations, iterations[1:])]
    if not interval_checks or any(interval != 100 for interval in interval_checks):
        raise FlowError(
            "CFD_FLOW_TRACKING_CONTINUITY_FAILED",
            "Combined samples are not strictly increasing at 100-iteration intervals",
        )
    combined: list[dict[str, Any]] = []
    previous_pressure: float | None = None
    previous_velocity: float | None = None
    for item in raw:
        pressure = float(item["avg_pressure_pa"])
        velocity = float(item["avg_velocity_m_s"])
        combined.append(
            {
                "iteration": int(item["iteration"]),
                "physical_time_s": float(item["physical_time_s"]),
                "avg_pressure_pa": pressure,
                "avg_velocity_m_s": velocity,
                "pressure_adjacent_delta": (
                    None if previous_pressure is None else pressure - previous_pressure
                ),
                "velocity_adjacent_delta": (
                    None if previous_velocity is None else velocity - previous_velocity
                ),
            }
        )
        previous_pressure = pressure
        previous_velocity = velocity
    qc = {
        "status": "PASS",
        "sampling_interval_iterations": 100,
        "samples_strictly_increasing": True,
        "shared_endpoint_iteration": shared_iteration,
        "shared_endpoint_deduplicated": True,
        "pressure_shared_endpoint_absolute_difference": pressure_difference,
        "pressure_shared_endpoint_tolerance": pressure_tolerance,
        "velocity_shared_endpoint_absolute_difference": velocity_difference,
        "velocity_shared_endpoint_tolerance": velocity_tolerance,
        "checks": endpoint_checks,
        "previous_sample_count": len(previous),
        "current_sample_count": len(current),
        "combined_unique_sample_count": len(combined),
    }
    return combined, qc


def write_combined_tracking_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fields = (
        "iteration",
        "physical_time_s",
        "avg_pressure_pa",
        "avg_velocity_m_s",
        "pressure_adjacent_delta",
        "velocity_adjacent_delta",
    )
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def official_equivalent_residuals(
    records: list[dict[str, Any]],
    *,
    nvals: int = OFFICIAL_EQUIVALENT_NVALS,
) -> list[dict[str, float | int]]:
    """Exactly emulate norm=average using previous values, excluding current."""

    if nvals <= 0:
        raise ValueError("nvals must be positive")
    pressure = np.asarray([item["avg_pressure_pa"] for item in records], dtype=float)
    velocity = np.asarray([item["avg_velocity_m_s"] for item in records], dtype=float)
    residuals: list[dict[str, float | int]] = []
    for index in range(nvals, len(records)):
        pressure_history_mean = float(np.mean(pressure[index - nvals : index]))
        velocity_history_mean = float(np.mean(velocity[index - nvals : index]))
        pressure_residual = abs(float(pressure[index]) - pressure_history_mean)
        velocity_residual = abs(float(velocity[index]) - velocity_history_mean)
        residuals.append(
            {
                "iteration": int(records[index]["iteration"]),
                "pressure_current": float(pressure[index]),
                "pressure_history_mean": pressure_history_mean,
                "pressure_residual_pa": pressure_residual,
                "pressure_ratio_to_threshold": pressure_residual / PRESSURE_THRESHOLD_PA,
                "velocity_current": float(velocity[index]),
                "velocity_history_mean": velocity_history_mean,
                "velocity_residual_m_s": velocity_residual,
                "velocity_ratio_to_threshold": velocity_residual / VELOCITY_THRESHOLD_M_S,
            }
        )
    return residuals


def write_official_equivalent_residuals_csv(
    path: Path,
    residuals: list[dict[str, Any]],
) -> None:
    fields = (
        "iteration",
        "pressure_current",
        "pressure_history_mean",
        "pressure_residual_pa",
        "pressure_ratio_to_threshold",
        "velocity_current",
        "velocity_history_mean",
        "velocity_residual_m_s",
        "velocity_ratio_to_threshold",
    )
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(residuals)


def _official_residual_trend(values: np.ndarray) -> dict[str, Any]:
    if len(values) < 5:
        raise FlowError(
            "CFD_FLOW_OFFICIAL_EQUIVALENT_RESIDUAL_INVALID",
            "Fewer than five official-equivalent residual points",
        )
    first_five = float(np.mean(values[:5]))
    last_five = float(np.mean(values[-5:]))
    ratio = math.inf if first_five == 0.0 and last_five > 0.0 else (
        0.0 if first_five == 0.0 else last_five / first_five
    )
    if ratio < 0.8:
        interpretation = "CLEARLY_DECAYING"
    elif ratio < 0.95:
        interpretation = "DECAYING"
    elif ratio <= 1.05:
        interpretation = "PLATEAU"
    else:
        interpretation = "GROWING"
    return {
        "first_residual": float(values[0]),
        "last_residual": float(values[-1]),
        "mean_first_5": first_five,
        "mean_last_5": last_five,
        "trend_ratio_last5_over_first5": ratio,
        "trend_interpretation": interpretation,
    }


def summarize_official_equivalent_residuals(
    records: list[dict[str, Any]],
    residuals: list[dict[str, Any]],
) -> dict[str, Any]:
    if not residuals:
        raise FlowError(
            "CFD_FLOW_OFFICIAL_EQUIVALENT_RESIDUAL_INVALID",
            "No complete 100-history residual point",
        )
    pressure_values = np.asarray(
        [item["pressure_residual_pa"] for item in residuals],
        dtype=float,
    )
    velocity_values = np.asarray(
        [item["velocity_residual_m_s"] for item in residuals],
        dtype=float,
    )
    final = residuals[-1]
    final_pressure = float(records[-1]["avg_pressure_pa"])
    final_velocity = float(records[-1]["avg_velocity_m_s"])
    return {
        "status": "PASS",
        "metric_name": "OFFLINE_OFFICIAL_EQUIVALENT_RESIDUAL",
        "equivalence_basis": "PINNED_TREELM_NORM_AVERAGE_NVALS_100_ABSOLUTE",
        "not_musubi_official_steady_termination": True,
        "musubi_internal_convergence_history_persists_across_restart": False,
        "spatial_reduction": "DOMAIN_AVERAGE",
        "sampling_interval_iterations": 100,
        "nvals": OFFICIAL_EQUIVALENT_NVALS,
        "residual_point_count": len(residuals),
        "pressure": {
            "final_average": final_pressure,
            "final_residual": float(final["pressure_residual_pa"]),
            "threshold": PRESSURE_THRESHOLD_PA,
            "ratio_to_threshold": float(final["pressure_ratio_to_threshold"]),
            "relative_to_final_average": (
                float(final["pressure_residual_pa"]) / abs(final_pressure)
            ),
            **_official_residual_trend(pressure_values),
        },
        "velocity": {
            "final_average": final_velocity,
            "final_residual": float(final["velocity_residual_m_s"]),
            "threshold": VELOCITY_THRESHOLD_M_S,
            "ratio_to_threshold": float(final["velocity_ratio_to_threshold"]),
            "relative_to_final_average": (
                float(final["velocity_residual_m_s"]) / abs(final_velocity)
            ),
            **_official_residual_trend(velocity_values),
        },
    }


def classify_official_equivalent(
    summary: dict[str, Any],
    *,
    numerical_runtime_error: bool,
) -> tuple[str, str]:
    """Classify only from official-equivalent residuals, never adjacent deltas."""

    pressure = summary["pressure"]
    velocity = summary["velocity"]
    pressure_ratio = float(pressure["ratio_to_threshold"])
    velocity_ratio = float(velocity["ratio_to_threshold"])
    pressure_trend = str(pressure["trend_interpretation"])
    velocity_trend = str(velocity["trend_interpretation"])
    if numerical_runtime_error or "GROWING" in (pressure_trend, velocity_trend):
        return "NUMERICAL_PROBLEM", "REVIEW_NUMERICAL_OR_BOUNDARY_CONDITION"
    if pressure_ratio <= 1.0 and velocity_ratio <= 1.0:
        return "OFFLINE_EQUIVALENT_STEADY_PASS", "RUN_SINGLE_FORMAL_STEADY_CONFIRMATION"
    if pressure_ratio <= 10.0 and velocity_ratio <= 10.0:
        return "OFFLINE_EQUIVALENT_NEAR_STEADY", "CONTINUE_FROM_RESTART_SHORT"
    decaying = {"DECAYING", "CLEARLY_DECAYING"}
    if pressure_trend in decaying and velocity_trend in decaying:
        return "STILL_CONVERGING", "CONTINUE_FROM_RESTART_LONGER"
    if "PLATEAU" in (pressure_trend, velocity_trend):
        return "CONVERGENCE_CRITERION_REVIEW_NEEDED", "REVIEW_CONVERGENCE_CRITERION"
    raise FlowError(
        "CFD_FLOW_OFFICIAL_EQUIVALENT_CLASSIFICATION_INVALID",
        "Residuals do not match an authorized convergence class",
    )


def build_project_criterion_review(
    *,
    outlet_gauge_pressures_pa: list[float] | tuple[float, ...],
    inlet_target_mean_velocity_m_s: float,
    current_pressure_residual_pa: float,
    current_velocity_residual_m_s: float,
) -> dict[str, Any]:
    """Build the single pre-registered 0.1% characteristic-scale criterion."""

    pressures = np.asarray(outlet_gauge_pressures_pa, dtype=float)
    if len(pressures) != 3 or not np.all(np.isfinite(pressures)):
        raise ValueError("Exactly three finite outlet gauge pressures are required")
    pressure_scale = float(np.max(pressures) - np.min(pressures))
    velocity_scale = float(inlet_target_mean_velocity_m_s)
    if pressure_scale <= 0.0 or velocity_scale <= 0.0:
        raise ValueError("Characteristic scales must be positive")
    pressure_threshold = PROJECT_STEADY_RELATIVE_FRACTION * pressure_scale
    velocity_threshold = PROJECT_STEADY_RELATIVE_FRACTION * velocity_scale
    return {
        "status": "PASS",
        "criterion_policy_name": PROJECT_STEADY_POLICY,
        "criterion_origin": "PROJECT_ENGINEERING_TOLERANCE",
        "criterion_is_literature_measured": False,
        "criterion_is_PROTEUS_requirement": False,
        "legacy_thresholds": {
            "status": "LEGACY_HEURISTIC_ABSOLUTE_THRESHOLDS",
            "origin": "EARLY_PROJECT_DIAGNOSTIC_ASSUMPTION",
            "pressure_pa": PRESSURE_THRESHOLD_PA,
            "velocity_m_s": VELOCITY_THRESHOLD_M_S,
            "not_proteus_requirement": True,
            "not_musubi_requirement": True,
            "not_experimentally_measured_tolerance": True,
        },
        "steady_relative_fraction": PROJECT_STEADY_RELATIVE_FRACTION,
        "pressure_characteristic_scale": {
            "definition": "OUTLET_GAUGE_PRESSURE_SPAN",
            "source_values_pa": [float(value) for value in pressures],
            "common_pressure_reference_excluded": True,
            "value_pa": pressure_scale,
        },
        "velocity_characteristic_scale": {
            "definition": "FIXED_INLET_TARGET_MEAN_VELOCITY",
            "not_current_domain_mean": True,
            "value_m_s": velocity_scale,
        },
        "derived_thresholds": {
            "pressure_pa": pressure_threshold,
            "velocity_m_s": velocity_threshold,
        },
        "current_offline_equivalent_residuals": {
            "pressure_pa": current_pressure_residual_pa,
            "velocity_m_s": current_velocity_residual_m_s,
            "pressure_ratio_to_project_threshold": (
                current_pressure_residual_pa / pressure_threshold
            ),
            "velocity_ratio_to_project_threshold": (
                current_velocity_residual_m_s / velocity_threshold
            ),
            "pre_run_status": "NOT_YET_PASS",
        },
        "musubi_convergence_contract": {
            "absolute": True,
            "norm": "average",
            "nvals": 100,
            "sampling_interval_iterations": 100,
            "shape": "all",
            "reduction": ["average", "average"],
            "internal_clock_ceiling": False,
            "threshold_sweep": False,
        },
        "rationale": (
            "Temporal convergence tolerance is fixed at 0.1% of physical "
            "characteristic scales so it remains below future percent-level "
            "spatial-uncertainty assessments; grid convergence is not yet established."
        ),
    }


def parse_official_steady_termination(text: str) -> dict[str, Any]:
    iterations = [
        int(value)
        for value in re.findall(r"(?i)Reached\s+steady\s+state\s+(\d+)\s+T", text)
    ]
    return {
        "official_steady_termination": bool(iterations),
        "evidence_phrase": "Reached steady state" if iterations else None,
        "iterations": sorted(set(iterations)),
        "confirmation_iteration": max(iterations) if iterations else None,
    }


def _window_metrics(values: np.ndarray, threshold: float) -> dict[str, Any]:
    absolute = np.abs(values)
    if len(absolute) < 20:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_TRACKING_INVALID",
            "INSUFFICIENT_TRACKING_WINDOW: fewer than 20 tracking deltas",
        )
    windows = {
        str(width): float(np.mean(absolute[-width:]))
        for width in (1, 5, 10, 20)
    }
    if len(absolute) >= 40:
        early_width = 20
        late_width = 20
        window_status = "FULL_WINDOW_DIAGNOSTIC"
    else:
        early_width = 10
        late_width = 10
        window_status = "REDUCED_WINDOW_DIAGNOSTIC"
    early = float(np.mean(absolute[:early_width]))
    late = float(np.mean(absolute[-late_width:]))
    ratio = math.inf if early == 0.0 and late > 0.0 else (0.0 if early == 0.0 else late / early)
    if ratio < 0.5:
        trend = "CLEARLY_DECAYING"
    elif ratio < 0.9:
        trend = "DECAYING"
    elif ratio <= 1.1:
        trend = "PLATEAU"
    else:
        trend = "GROWING"
    return {
        "sample_count": len(absolute) + 1,
        "delta_count": len(absolute),
        "window_status": window_status,
        "last_sample_delta": float(absolute[-1]),
        "mean_absolute_delta_by_last_n_samples": windows,
        "last_20_mean_absolute_delta": windows["20"],
        "last_20_max_absolute_delta": float(np.max(absolute[-20:])),
        "threshold": threshold,
        "ratio_to_threshold": windows["20"] / threshold,
        "early_20_mean_absolute_delta": early,
        "trend_early_window": early_width,
        "trend_late_window": late_width,
        "trend_ratio_late_over_early": ratio,
        "trend_interpretation": trend,
    }


def summarize_convergence(records: list[dict[str, Any]]) -> dict[str, Any]:
    pressure_delta = np.asarray(
        [item["delta_pressure_vs_previous"] for item in records[1:]], dtype=float
    )
    velocity_delta = np.asarray(
        [item["delta_velocity_vs_previous"] for item in records[1:]], dtype=float
    )
    return {
        "status": "PASS",
        "sample_count": len(records),
        "pressure": _window_metrics(pressure_delta, PRESSURE_THRESHOLD_PA),
        "velocity": _window_metrics(velocity_delta, VELOCITY_THRESHOLD_M_S),
        "thresholds_unchanged": True,
        "estimate_policy": "DIAGNOSTIC_ESTIMATE_ONLY_NO_LONG_EXTRAPOLATION",
    }


def create_convergence_figure(records: list[dict[str, Any]], path: Path) -> None:
    iterations = np.asarray([item["iteration"] for item in records[1:]], dtype=int)
    pressure = np.abs(
        np.asarray([item["delta_pressure_vs_previous"] for item in records[1:]], dtype=float)
    )
    velocity = np.abs(
        np.asarray([item["delta_velocity_vs_previous"] for item in records[1:]], dtype=float)
    )
    figure, axes = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True, constrained_layout=True)
    for axis, values, threshold, label, color in (
        (axes[0], pressure, PRESSURE_THRESHOLD_PA, "|Δ average pressure| (Pa)", "tab:blue"),
        (axes[1], velocity, VELOCITY_THRESHOLD_M_S, "|Δ average velocity| (m/s)", "tab:orange"),
    ):
        safe = np.maximum(values, np.finfo(float).tiny)
        axis.semilogy(iterations, safe, color=color, linewidth=1.2)
        axis.axhline(threshold, color="black", linestyle="--", linewidth=1.0, label="fixed threshold")
        axis.set_ylabel(label)
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(loc="best")
    axes[1].set_xlabel("Iteration")
    axes[0].set_title("Restart convergence diagnostic (100-iteration sampling)")
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _classify(
    field_qc: dict[str, Any],
    flow_qc: dict[str, Any],
    convergence: dict[str, Any],
) -> tuple[str, str, dict[str, Any]]:
    finite = bool(field_qc["finite"])
    mach_ok = bool(field_qc["actual_lattice_mach_below_0_05"])
    directions_ok = bool(flow_qc["flow_direction_pass"])
    mass_error = float(flow_qc["mass_conservation_relative_error"])
    pressure = convergence["pressure"]
    velocity = convergence["velocity"]
    both_growing = (
        pressure["trend_interpretation"] == "GROWING"
        and velocity["trend_interpretation"] == "GROWING"
    )
    physically_plausible = finite and mach_ok and directions_ok and mass_error <= 0.20
    mass_close = mass_error <= 0.05
    still_evolving = (
        pressure["ratio_to_threshold"] > 10.0
        or velocity["ratio_to_threshold"] > 10.0
    )
    plateau = (
        pressure["trend_interpretation"] == "PLATEAU"
        or velocity["trend_interpretation"] == "PLATEAU"
    )
    criterion_floor = mass_error <= 0.01 and physically_plausible and plateau and still_evolving

    if not physically_plausible or both_growing:
        classification = "NUMERICAL_OR_BC_PROBLEM"
        next_step = "REVIEW_NUMERICAL_OR_BOUNDARY_CONDITION"
    elif criterion_floor:
        classification = "CONVERGENCE_CRITERION_REVIEW_NEEDED"
        next_step = "REVIEW_CONVERGENCE_CRITERION"
    elif (
        mass_error <= 0.01
        and pressure["ratio_to_threshold"] <= 10.0
        and velocity["ratio_to_threshold"] <= 10.0
        and pressure["trend_ratio_late_over_early"] < 0.9
        and velocity["trend_ratio_late_over_early"] < 0.9
    ):
        classification = "NEAR_STEADY_CONTINUE_FROM_RESTART"
        next_step = "CONTINUE_FROM_RESTART_SHORT"
    else:
        classification = "PHYSICALLY_REASONABLE_BUT_STILL_DEVELOPING"
        next_step = "CONTINUE_FROM_RESTART_LONGER"
    evidence = {
        "physically_plausible": physically_plausible,
        "mass_conservation_close_within_5_percent": mass_close,
        "flow_still_clearly_evolving": still_evolving,
        "convergence_criterion_possibly_below_numerical_floor": (
            "YES" if criterion_floor else ("NO" if still_evolving and not plateau else "INSUFFICIENT_EVIDENCE")
        ),
    }
    return classification, next_step, evidence


def _parse_terminal_iteration(text: str) -> int:
    values = [int(value) for value in re.findall(r"iterations:\s*(\d+)", text)]
    if not values:
        raise FlowError("CFD_FLOW_DIAGNOSTIC_MUSUBI_FAILED", "No iteration evidence in short run")
    return max(values)


def _runtime_numerical_errors(
    text: str,
    returncode: int,
    *,
    wrapper_timeout: bool = False,
) -> dict[str, Any]:
    """Separate solver numerics from diagnostic I/O and process-control failures."""

    main_loop_started = "Starting Musubi MAIN loop" in text
    nan_detected = bool(re.search(r"(?i)(?<![A-Za-z])nan(?![A-Za-z])", text))
    inf_detected = bool(
        re.search(r"(?i)(?<![A-Za-z])[-+]?inf(?:inity)?(?![A-Za-z])", text)
    )
    floating_point_exception = bool(
        re.search(r"(?i)floating[-\s]+point\s+exception", text)
    )
    segfault_detected = bool(
        re.search(r"(?i)(?:segmentation\s+fault|segfault|sigsegv)", text)
    )
    fortran_io_error = bool(
        re.search(r"(?is)hvs_ascii[^\n]*.*?fortran\s+runtime\s+error:\s*end\s+of\s+record", text)
    )
    time_control_termination = bool(
        re.search(r"(?i)reached\s+maximal\s+wall\s+clock(?:\s+running)?\s+time", text)
    )
    negative_density = bool(
        re.search(r"(?i)(?:negative\s+density|density[^\n]*negative)", text)
    )
    numerical_abort = bool(
        re.search(
            r"(?i)(?:numerical\s+(?:failure|abort)|density\s+explosion|"
            r"velocity\s+explosion)",
            text,
        )
    )
    numerical_evidence = (
        nan_detected
        or inf_detected
        or floating_point_exception
        or segfault_detected
        or negative_density
        or numerical_abort
    )
    numerical_runtime_error = main_loop_started and numerical_evidence
    if fortran_io_error and not numerical_runtime_error:
        failure_classification = "DIAGNOSTIC_IO_ERROR"
    elif time_control_termination and not numerical_runtime_error:
        failure_classification = "TIME_CONTROL_TERMINATION"
    elif wrapper_timeout and not numerical_runtime_error:
        failure_classification = "WRAPPER_TIMEOUT"
    elif numerical_runtime_error:
        failure_classification = "NUMERICAL_RUNTIME_ERROR"
    elif returncode != 0:
        failure_classification = "NONZERO_RETURN_WITHOUT_NUMERICAL_EVIDENCE"
    else:
        failure_classification = "NONE"
    return {
        "nan_detected": nan_detected,
        "inf_detected": inf_detected,
        "floating_point_exception_detected": floating_point_exception,
        "segfault_detected": segfault_detected,
        "fortran_io_error_detected": fortran_io_error,
        "time_control_termination_detected": time_control_termination,
        "wrapper_timeout_detected": wrapper_timeout,
        "negative_density_warning_detected": negative_density,
        "solver_numerical_abort_detected": numerical_abort,
        "solver_main_loop_started": main_loop_started,
        "returncode_nonzero": returncode != 0,
        "numerical_runtime_error": numerical_runtime_error,
        "failure_classification": failure_classification,
    }


def _parse_restart_start_iteration(text: str) -> int:
    match = re.search(
        r"Restarting\s+from\s+point\s+in\s+time:\s*.*?iterations:\s*(\d+)",
        text,
        re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        raise FlowError(
            "CFD_FLOW_RESTART_START_INVALID",
            "No explicit restart-start iteration in Musubi runtime log",
        )
    return int(match.group(1))


def classify_convergence_probe(
    convergence: dict[str, Any],
    *,
    numerical_runtime_error: bool,
) -> tuple[str, str]:
    """Apply only the four authorized convergence-probe decisions."""

    pressure = convergence["pressure"]
    velocity = convergence["velocity"]
    pressure_ratio = float(pressure["ratio_to_threshold"])
    velocity_ratio = float(velocity["ratio_to_threshold"])
    pressure_trend = str(pressure["trend_interpretation"])
    velocity_trend = str(velocity["trend_interpretation"])

    if numerical_runtime_error or "GROWING" in (pressure_trend, velocity_trend):
        return "NUMERICAL_PROBLEM", "REVIEW_NUMERICAL_OR_BOUNDARY_CONDITION"
    if pressure_ratio <= 10.0 and velocity_ratio <= 10.0:
        return "NEAR_STEADY", "CONTINUE_FROM_RESTART_SHORT"

    decaying = {"DECAYING", "CLEARLY_DECAYING"}
    if (
        (pressure_ratio > 10.0 or velocity_ratio > 10.0)
        and pressure_trend in decaying
        and velocity_trend in decaying
    ):
        return "STILL_CONVERGING", "CONTINUE_FROM_RESTART_LONGER"
    if (
        pressure_trend == "PLATEAU"
        and velocity_trend == "PLATEAU"
        and (pressure_ratio > 1.0 or velocity_ratio > 1.0)
    ):
        return "CONVERGENCE_CRITERION_REVIEW_NEEDED", "REVIEW_CONVERGENCE_CRITERION"
    raise FlowError(
        "CFD_FLOW_CONVERGENCE_PROBE_UNCLASSIFIED",
        "Tracking does not satisfy any authorized convergence classification",
    )


def _restart_manifest(directory: Path, end_iteration: int) -> dict[str, Any]:
    manifest = _directory_manifest(directory)
    headers = sorted(directory.glob("*lastHeader.lua"))
    if not headers:
        raise FlowError("CFD_FLOW_DIAGNOSTIC_MUSUBI_FAILED", "No diagnostic continuation restart header")
    return {
        "status": "DIAGNOSTIC_CONTINUATION_RESTART",
        "iteration": end_iteration,
        "not_frozen_musubi_solution": True,
        "files": manifest,
    }


def run_diagnostic(project_root: Path) -> dict[str, Any]:
    """Execute exactly one harvest and one 10k-step restart continuation."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = root / "outputs" / "cfd_flow"
    existing = sorted(output_root.glob("musubi_diagnostic_anchor003274_*"))
    if existing:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_RUN_LIMIT_REACHED",
            f"A diagnostic run already exists: {existing[-1]}",
        )
    source_run = output_root / SOURCE_RUN_NAME
    frozen_seeder_run = output_root / FROZEN_SEEDER_NAME
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"musubi_diagnostic_anchor003274_{timestamp}"
    current_snapshot = run_root / "current_snapshot"
    tracking = run_root / "tracking"
    qc_dir = run_root / "qc"
    for directory in (current_snapshot, tracking, qc_dir):
        directory.mkdir(parents=True, exist_ok=False)

    summary: dict[str, Any] = {
        "status": "RUNNING_ZERO_SOLVER_COST_DIAGNOSTIC",
        "run_root": str(run_root),
        "production_code_modified": False,
        "source_run": str(source_run),
        "frozen_seeder_run": str(frozen_seeder_run),
        "seeder_run_count": 0,
        "long_musubi_run_count": 0,
        "short_restart_musubi_run_count": 0,
        "additional_musubi_iterations": 0,
        "diagnostic_harvester_run_count": 0,
        "grid_sweep": False,
        "physics_or_bc_changed": False,
        "performance_optimization": "PERFORMANCE_OPTIMIZATION_DEFERRED",
        "proteus_import_ready_claimed": False,
    }
    summary_path = qc_dir / "diagnostic_summary.json"
    try:
        source_manifest_before = _directory_manifest(source_run)
        stdout_path = source_run / "musubi" / "musubi_stdout.log"
        density_records = parse_total_density_log(
            stdout_path.read_text(encoding="utf-8", errors="replace")
        )
        write_total_density_csv(qc_dir / "diagnostic_total_density.csv", density_records)
        density_summary = summarize_total_density(density_records)
        write_json(qc_dir / "total_density_summary.json", density_summary)
        restart_qc = validate_restart(source_run, frozen_seeder_run)
        write_json(qc_dir / "restart_validation.json", restart_qc)
        input_reference = {
            "status": "PASS",
            "source_run": str(source_run),
            "source_run_summary": str(source_run / "qc" / "run_summary.json"),
            "restart_validation": restart_qc,
            "total_density_summary": density_summary,
            "official_restart_harvest_basis": {
                "pinned_musubi_source": "/home/lzy/apes-pinned/musubi_official",
                "pinned_source_commit": "4e8b277b66226277171ef93bf054d21270812793",
                "runtime_revision": "81f8c4f13772",
                "example": "mus/examples/fluid_incompressible/benchmark/Pipe/PIP_Simple/mus_harvester.lua",
                "contract": "restart.read points to the lastHeader Lua",
            },
        }
        write_json(run_root / "input_reference.json", input_reference)
        summary.update(
            {
                "status": "ZERO_SOLVER_COST_DIAGNOSTIC_PASS",
                "restart_validation": restart_qc,
                "total_density": density_summary,
            }
        )
        write_json(summary_path, summary)

        inputs = load_flow_inputs(config.paths.source_surface_run)
        partition = load_frozen_surface_partition(
            inputs,
            frozen_seeder_run / "geometry" / "geometry_solver_m",
        )
        bc = load_boundary_conditions(inputs.boundary_conditions)
        scaling = compute_lattice_scaling(
            config,
            bc,
            partition.patch("inlet").area_um2 * 1.0e-12,
        )
        environment = inspect_apes_environment(config.apes)
        if environment.status != "PASS":
            raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Pinned APES tools are unavailable")

        raw_snapshot = current_snapshot / "raw_vtk"
        raw_snapshot.mkdir()
        harvest_lua = current_snapshot / "diagnostic_harvester.lua"
        harvest_lua.write_text(
            generate_diagnostic_harvester_lua(
                windows_to_wsl(raw_snapshot, config.apes.wsl_distribution)
            ),
            encoding="utf-8",
        )
        luac_harvest = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=current_snapshot,
            command=[
                str(environment.binaries["lua_compiler"]),
                "-p",
                "diagnostic_harvester.lua",
            ],
            stdout_path=current_snapshot / "luac_harvester_stdout.log",
            stderr_path=current_snapshot / "luac_harvester_stderr.log",
            timeout_s=30,
        )
        if luac_harvest.returncode != 0:
            raise FlowError(HARVEST_FAILED, "Diagnostic harvester Lua syntax failed")
        harvest_config_wsl = windows_to_wsl(harvest_lua, config.apes.wsl_distribution)
        summary["diagnostic_harvester_run_count"] = 1
        write_json(summary_path, summary)
        harvest = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=source_run / "musubi",
            command=[str(environment.binaries["mus_harvesting"]), harvest_config_wsl],
            stdout_path=current_snapshot / "harvester_stdout.log",
            stderr_path=current_snapshot / "harvester_stderr.log",
            timeout_s=300,
        )
        summary["diagnostic_harvester_wall_time_s"] = harvest.wall_time_s
        if harvest.returncode != 0:
            summary.update(
                {
                    "status": HARVEST_FAILED,
                    "harvester_returncode": harvest.returncode,
                    "harvester_failure_stage": (
                        "AFTER_RESTART_READ_BEFORE_COMPLETE_VTK_DATASET_WRITE"
                    ),
                    "short_restart_run_blocked_by_hard_rule": True,
                }
            )
            write_json(summary_path, summary)
            raise FlowError(HARVEST_FAILED, f"mus_harvesting return code {harvest.returncode}")
        raw_vtk = _find_harvested_vtk(raw_snapshot)
        current_vtu = current_snapshot / "current_flow_field.vtu"
        grid, field_qc, flow_qc = current_snapshot_qc(
            raw_vtk,
            current_vtu,
            partition=partition,
            scaling=scaling,
            target_q_m3_s=bc.inlet_flow_m3_s,
        )
        write_json(qc_dir / "current_snapshot_qc.json", field_qc)
        write_json(qc_dir / "current_flow_balance_qc.json", flow_qc)
        create_velocity_snapshot(grid, qc_dir / "current_velocity_snapshot.png")
        summary.update(
            {
                "status": "CURRENT_DIAGNOSTIC_HARVEST_PASS",
                "current_diagnostic_harvest": "PASS",
                "current_snapshot_vtu": str(current_vtu),
                "current_snapshot_qc": field_qc,
                "current_flow_balance_qc": flow_qc,
            }
        )
        write_json(summary_path, summary)

        production_lua_path = source_run / "musubi" / "musubi.lua"
        production_lua = production_lua_path.read_text(encoding="utf-8")
        diagnostic_restart = tracking / "diagnostic_restart"
        diagnostic_restart.mkdir()
        tracking_wsl = windows_to_wsl(tracking, config.apes.wsl_distribution)
        diagnostic_lua = generate_diagnostic_musubi_lua(
            production_lua,
            pressure_folder_wsl=tracking_wsl,
            velocity_folder_wsl=tracking_wsl,
            restart_folder_wsl=windows_to_wsl(
                diagnostic_restart,
                config.apes.wsl_distribution,
            ),
            timing_file_wsl=f"{tracking_wsl.rstrip('/')}/diagnostic_timing.res",
        )
        diagnostic_lua_path = run_root / "diagnostic_musubi.lua"
        diagnostic_lua_path.write_text(diagnostic_lua, encoding="utf-8")
        lua_contract = diagnostic_lua_contract(production_lua, diagnostic_lua)
        write_json(qc_dir / "diagnostic_musubi_contract.json", lua_contract)
        if lua_contract["status"] != "PASS":
            raise FlowError("CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID", "Frozen physics/BC contract failed")
        luac_musubi = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[
                str(environment.binaries["lua_compiler"]),
                "-p",
                "diagnostic_musubi.lua",
            ],
            stdout_path=tracking / "luac_musubi_stdout.log",
            stderr_path=tracking / "luac_musubi_stderr.log",
            timeout_s=30,
        )
        if luac_musubi.returncode != 0:
            raise FlowError("CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID", "Diagnostic Musubi Lua syntax failed")

        diagnostic_config_wsl = windows_to_wsl(
            diagnostic_lua_path,
            config.apes.wsl_distribution,
        )
        mpi_command = [
            "env",
            "OMPI_ALLOW_RUN_AS_ROOT=1",
            "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            str(environment.binaries["mpi_launcher"]),
            "-np",
            str(environment.mpi_ranks),
            str(environment.binaries["musubi"]),
            diagnostic_config_wsl,
        ]
        summary["short_restart_musubi_run_count"] = 1
        summary["restart_continuation_start_iteration"] = START_ITERATION
        write_json(summary_path, summary)
        short_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=source_run / "musubi",
            command=mpi_command,
            stdout_path=tracking / "musubi_stdout.log",
            stderr_path=tracking / "musubi_stderr.log",
            timeout_s=DIAGNOSTIC_WALLCLOCK_S + 60,
        )
        summary["short_run_wall_time_s"] = short_run.wall_time_s
        short_stdout = short_run.stdout_path.read_text(encoding="utf-8", errors="replace")
        short_stderr = short_run.stderr_path.read_text(encoding="utf-8", errors="replace")
        if short_run.returncode != 0:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_MUSUBI_FAILED",
                f"Short Musubi return code {short_run.returncode}",
            )
        end_iteration = _parse_terminal_iteration(short_stdout + "\n" + short_stderr)
        additional = end_iteration - START_ITERATION
        if not 0 < additional <= ADDITIONAL_ITERATIONS:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_RUN_LIMIT_VIOLATED",
                f"Unexpected additional iterations: {additional}",
            )
        if re.search(r"(?i)(?<![A-Za-z])(?:nan|inf)(?![A-Za-z])", short_stdout + "\n" + short_stderr):
            raise FlowError("CFD_FLOW_DIAGNOSTIC_NUMERICAL_FAILED", "NaN/Inf in short-run logs")
        summary["restart_continuation_end_iteration"] = end_iteration
        summary["additional_musubi_iterations"] = additional
        summary["nan_or_inf"] = False

        tracking_records = parse_tracking(tracking, scaling.dt_s)
        write_tracking_csv(qc_dir / "restart_convergence_diagnostic.csv", tracking_records)
        convergence = summarize_convergence(tracking_records)
        write_json(qc_dir / "convergence_distance_summary.json", convergence)
        create_convergence_figure(tracking_records, qc_dir / "convergence_diagnostic.png")
        continuation_restart = _restart_manifest(diagnostic_restart, end_iteration)
        write_json(qc_dir / "diagnostic_restart_manifest.json", continuation_restart)
        classification, next_step, evidence = _classify(field_qc, flow_qc, convergence)

        source_manifest_after = _directory_manifest(source_run)
        historical_read_only = source_manifest_before == source_manifest_after
        if not historical_read_only:
            raise FlowError(
                "CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED",
                "The source formal run changed during diagnostics",
            )
        summary.update(
            {
                "status": DIAGNOSTIC_STATUS,
                "historical_outputs_read_only": historical_read_only,
                "tracking": convergence,
                "convergence_classification": classification,
                "classification_evidence": evidence,
                "diagnostic_restart_saved": True,
                "diagnostic_restart": continuation_restart,
                "figures": [
                    str(qc_dir / "convergence_diagnostic.png"),
                    str(qc_dir / "current_velocity_snapshot.png"),
                ],
                "next_recommendation": next_step,
                "environment": asdict(environment),
            }
        )
        write_json(summary_path, summary)
        return summary
    except Exception as error:
        if isinstance(error, FlowError):
            summary["status"] = error.status
        else:
            summary["status"] = "CFD_FLOW_DIAGNOSTIC_INTERNAL_ERROR"
        summary["failure"] = str(error)
        write_json(summary_path, summary)
        raise


def _read_json_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _finalize_convergence_probe(project_root: Path, run_root: Path) -> dict[str, Any]:
    """Post-process an already completed short run without launching any executable."""

    root = Path(project_root).resolve()
    output_root = root / "outputs" / "cfd_flow"
    source_run = output_root / SOURCE_RUN_NAME
    tracking = run_root / "tracking"
    qc_dir = run_root / "qc"
    summary_path = qc_dir / "convergence_probe_summary.json"
    summary = _read_json_object(summary_path)

    stdout_path = tracking / "musubi_stdout.log"
    stderr_path = tracking / "musubi_stderr.log"
    short_stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    short_stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
    runtime_text = short_stdout + "\n" + short_stderr
    returncode = int(summary["musubi_returncode"])
    numerical_errors = _runtime_numerical_errors(runtime_text, returncode)
    write_json(qc_dir / "runtime_numerical_error_scan.json", numerical_errors)

    logged_start_iteration = _parse_restart_start_iteration(short_stdout)
    if logged_start_iteration != START_ITERATION:
        raise FlowError(
            "CFD_FLOW_RESTART_START_INVALID",
            f"Runtime loaded iteration {logged_start_iteration}, expected {START_ITERATION}",
        )

    original_stdout = source_run / "musubi" / "musubi_stdout.log"
    original_density = parse_total_density_log(
        original_stdout.read_text(encoding="utf-8", errors="replace")
    )
    if returncode != 0:
        full_density_summary = summarize_total_density(original_density)
        density_summary = {
            "status": "PASS_ORIGINAL_RUN_ONLY",
            "interpretation": "DIAGNOSTIC_ONLY_NOT_A_STEADY_HARD_GATE",
            "continuation_sample_count": 0,
            "last_1000_iterations": full_density_summary["windows"]["1000"],
            "last_5000_iterations": full_density_summary["windows"]["5000"],
        }
        write_total_density_csv(qc_dir / "combined_total_density.csv", original_density)
        write_json(qc_dir / "total_density_summary.json", density_summary)
        unavailable = "NOT_AVAILABLE_RUN_ABORTED_BEFORE_FIRST_TRACKING_SAMPLE"
        tracking_unavailable = {
            "status": unavailable,
            "pressure": {
                "last_sample_delta": unavailable,
                "last_20_mean_absolute_delta": unavailable,
                "threshold": PRESSURE_THRESHOLD_PA,
                "ratio_to_threshold": unavailable,
                "trend_ratio_late_over_early": unavailable,
                "trend_interpretation": unavailable,
            },
            "velocity": {
                "last_sample_delta": unavailable,
                "last_20_mean_absolute_delta": unavailable,
                "threshold": VELOCITY_THRESHOLD_M_S,
                "ratio_to_threshold": unavailable,
                "trend_ratio_late_over_early": unavailable,
                "trend_interpretation": unavailable,
            },
        }
        write_json(qc_dir / "convergence_distance_summary.json", tracking_unavailable)
        historical_read_only = _read_json_object(
            qc_dir / "source_manifest_before.json"
        ) == _directory_manifest(source_run)
        if not historical_read_only:
            raise FlowError(
                "CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED",
                "The source formal run changed during the convergence probe",
            )
        summary.update(
            {
                "status": "CFD_FLOW_CONVERGENCE_PROBE_MUSUBI_FAILED",
                "restart_continuation_start_iteration": logged_start_iteration,
                "restart_continuation_end_iteration": logged_start_iteration,
                "additional_musubi_iterations": 0,
                "numerical_runtime_error": numerical_errors,
                "runtime_failure_stage": "ASCII_TRACKING_INITIALIZATION_BEFORE_MAIN_LOOP",
                "runtime_failure": "FORTRAN_END_OF_RECORD",
                "tracking": tracking_unavailable,
                "total_density": density_summary,
                "diagnostic_restart": "NOT_AVAILABLE_RUN_ABORTED_BEFORE_RESTART_WRITE",
                "diagnostic_restart_path": "NOT_AVAILABLE_RUN_ABORTED_BEFORE_RESTART_WRITE",
                "historical_outputs_read_only": True,
                "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
                "convergence_classification": "NUMERICAL_PROBLEM",
                "next_recommendation": "REVIEW_NUMERICAL_OR_BOUNDARY_CONDITION",
                "field_metrics": {
                    "qin": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                    "qout": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                    "mass_error": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                    "mach": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                    "outlet_pressure_spatial_mean": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                },
                "figures": [],
                "failure": "Musubi return code 2 during ASCII tracking initialization",
            }
        )
        write_json(qc_dir / "convergence_probe_failure_qc.json", summary)
        write_json(summary_path, summary)
        return summary

    short_density = parse_total_density_log(short_stdout)
    start_iteration = int(short_density[0]["iteration"])
    end_iteration = int(short_density[-1]["iteration"])
    if start_iteration != START_ITERATION:
        raise FlowError(
            "CFD_FLOW_RESTART_START_INVALID",
            f"Runtime began at iteration {start_iteration}, expected {START_ITERATION}",
        )
    additional = end_iteration - start_iteration
    if not 0 < additional <= ADDITIONAL_ITERATIONS:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_RUN_LIMIT_VIOLATED",
            f"Unexpected additional iterations: {additional}",
        )
    production_lua = (source_run / "musubi" / "musubi.lua").read_text(encoding="utf-8")
    dt_s = _extract_lua_number(production_lua, "dt")
    tracking_records = parse_tracking(tracking, dt_s)
    write_tracking_csv(qc_dir / "restart_convergence_diagnostic.csv", tracking_records)
    convergence = summarize_convergence(tracking_records)
    write_json(qc_dir / "convergence_distance_summary.json", convergence)
    create_convergence_figure(tracking_records, qc_dir / "convergence_diagnostic.png")

    combined_by_iteration = {
        int(item["iteration"]): item for item in (*original_density, *short_density)
    }
    combined_density = [combined_by_iteration[key] for key in sorted(combined_by_iteration)]
    write_total_density_csv(qc_dir / "combined_total_density.csv", combined_density)
    full_density_summary = summarize_total_density(combined_density)
    density_summary = {
        "status": "PASS",
        "interpretation": "DIAGNOSTIC_ONLY_NOT_A_STEADY_HARD_GATE",
        "last_1000_iterations": full_density_summary["windows"]["1000"],
        "last_5000_iterations": full_density_summary["windows"]["5000"],
    }
    write_json(qc_dir / "total_density_summary.json", density_summary)

    diagnostic_restart_dir = tracking / "diagnostic_restart"
    continuation_restart = _restart_manifest(diagnostic_restart_dir, end_iteration)
    last_headers = sorted(diagnostic_restart_dir.glob("*lastHeader.lua"))
    continuation_header = last_headers[-1]
    continuation_iteration = int(
        _extract_lua_number(continuation_header.read_text(encoding="utf-8"), "iter")
    )
    if continuation_iteration != end_iteration:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_MUSUBI_FAILED",
            "Diagnostic restart iteration does not match terminal runtime iteration",
        )
    continuation_restart["last_header"] = str(continuation_header)
    write_json(qc_dir / "diagnostic_restart_manifest.json", continuation_restart)

    classification, next_step = classify_convergence_probe(
        convergence,
        numerical_runtime_error=bool(numerical_errors["numerical_runtime_error"]),
    )
    source_manifest_before = _read_json_object(qc_dir / "source_manifest_before.json")
    historical_read_only = source_manifest_before == _directory_manifest(source_run)
    if not historical_read_only:
        raise FlowError(
            "CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED",
            "The source formal run changed during the convergence probe",
        )

    summary.update(
        {
            "status": PROBE_STATUS,
            "restart_continuation_start_iteration": start_iteration,
            "restart_continuation_end_iteration": end_iteration,
            "additional_musubi_iterations": additional,
            "numerical_runtime_error": numerical_errors,
            "tracking": convergence,
            "total_density": density_summary,
            "diagnostic_restart": continuation_restart,
            "diagnostic_restart_path": str(continuation_header),
            "diagnostic_restart_status": "DIAGNOSTIC_CONTINUATION_RESTART",
            "historical_outputs_read_only": True,
            "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
            "convergence_classification": classification,
            "next_recommendation": next_step,
            "field_metrics": {
                "qin": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "qout": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "mass_error": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "mach": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "outlet_pressure_spatial_mean": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
            },
            "figures": [str(qc_dir / "convergence_diagnostic.png")],
        }
    )
    summary.pop("failure", None)
    write_json(summary_path, summary)
    return summary


def _run_convergence_probe_obsolete(project_root: Path) -> dict[str, Any]:
    """Run at most one 5k-step restart continuation, with ASCII tracking only."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = root / "outputs" / "cfd_flow"
    existing = sorted(output_root.glob("musubi_convergence_probe_anchor003274_*"))
    if existing:
        latest = existing[-1]
        summary_path = latest / "qc" / "convergence_probe_summary.json"
        if summary_path.is_file():
            prior = _read_json_object(summary_path)
            if int(prior.get("short_restart_musubi_run_count", 0)) == 1:
                if prior.get("status") == PROBE_STATUS:
                    return prior
                return _finalize_convergence_probe(root, latest)
        raise FlowError(
            "CFD_FLOW_CONVERGENCE_PROBE_RUN_LIMIT_REACHED",
            f"A convergence-probe directory already exists: {latest}",
        )

    source_run = output_root / SOURCE_RUN_NAME
    frozen_seeder_run = output_root / FROZEN_SEEDER_NAME
    restart_qc = validate_restart(source_run, frozen_seeder_run)
    production_lua_path = source_run / "musubi" / "musubi.lua"
    production_lua = production_lua_path.read_text(encoding="utf-8")
    environment = inspect_apes_environment(config.apes)
    if environment.status != "PASS" or not environment.binaries["musubi"]:
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Pinned Musubi environment is unavailable")
    if environment.mpi_ranks != 8:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
            f"Expected 8 MPI ranks, found {environment.mpi_ranks}",
        )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"musubi_convergence_probe_anchor003274_{timestamp}"
    tracking = run_root / "tracking"
    qc_dir = run_root / "qc"
    diagnostic_restart = tracking / "diagnostic_restart"
    for directory in (diagnostic_restart, qc_dir):
        directory.mkdir(parents=True, exist_ok=False)

    source_manifest_before = _directory_manifest(source_run)
    write_json(qc_dir / "source_manifest_before.json", source_manifest_before)
    summary: dict[str, Any] = {
        "status": "CONVERGENCE_PROBE_PREFLIGHT",
        "run_root": str(run_root),
        "production_code_modified": False,
        "source_run": str(source_run),
        "frozen_seeder_run": str(frozen_seeder_run),
        "seeder_run_count": 0,
        "diagnostic_harvester_run_count": 0,
        "short_restart_musubi_run_count": 0,
        "long_musubi_run_count": 0,
        "additional_musubi_iterations": 0,
        "physics_or_bc_changed": False,
        "vtk_output": False,
        "restart_validation": restart_qc,
        "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
        "environment": asdict(environment),
    }
    summary_path = qc_dir / "convergence_probe_summary.json"
    write_json(summary_path, summary)
    write_json(
        run_root / "input_reference.json",
        {
            "status": "PASS",
            "source_run": str(source_run),
            "restart_header": str(
                source_run / "musubi" / "restart" / "roi003274_steady_lbm_lastHeader.lua"
            ),
            "restart_validation": restart_qc,
            "source_manifest_before": str(qc_dir / "source_manifest_before.json"),
        },
    )

    try:
        tracking_wsl = windows_to_wsl(tracking, config.apes.wsl_distribution)
        diagnostic_lua = generate_diagnostic_musubi_lua(
            production_lua,
            pressure_folder_wsl=tracking_wsl,
            velocity_folder_wsl=tracking_wsl,
            restart_folder_wsl=windows_to_wsl(
                diagnostic_restart,
                config.apes.wsl_distribution,
            ),
            timing_file_wsl=f"{tracking_wsl.rstrip('/')}/diagnostic_timing.res",
        )
        diagnostic_lua_path = run_root / "diagnostic_musubi.lua"
        diagnostic_lua_path.write_text(diagnostic_lua, encoding="utf-8")
        lua_contract = diagnostic_lua_contract(production_lua, diagnostic_lua)
        write_json(qc_dir / "diagnostic_musubi_contract.json", lua_contract)
        if lua_contract["status"] != "PASS":
            raise FlowError("CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID", "Frozen physics/BC contract failed")

        luac = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[str(environment.binaries["lua_compiler"]), "-p", "diagnostic_musubi.lua"],
            stdout_path=tracking / "luac_musubi_stdout.log",
            stderr_path=tracking / "luac_musubi_stderr.log",
            timeout_s=30,
        )
        if luac.returncode != 0:
            raise FlowError("CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID", "Diagnostic Musubi Lua syntax failed")

        diagnostic_config_wsl = windows_to_wsl(
            diagnostic_lua_path,
            config.apes.wsl_distribution,
        )
        mpi_command = [
            "env",
            "OMPI_ALLOW_RUN_AS_ROOT=1",
            "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            str(environment.binaries["mpi_launcher"]),
            "-np",
            "8",
            str(environment.binaries["musubi"]),
            diagnostic_config_wsl,
        ]
        summary.update(
            {
                "status": "CONVERGENCE_PROBE_MUSUBI_RUNNING",
                "short_restart_musubi_run_count": 1,
                "restart_continuation_start_iteration": START_ITERATION,
                "requested_end_iteration": END_ITERATION,
                "requested_additional_iterations": ADDITIONAL_ITERATIONS,
                "maximum_wall_clock_s": DIAGNOSTIC_WALLCLOCK_S,
            }
        )
        write_json(summary_path, summary)
        short_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=source_run / "musubi",
            command=mpi_command,
            stdout_path=tracking / "musubi_stdout.log",
            stderr_path=tracking / "musubi_stderr.log",
            timeout_s=DIAGNOSTIC_WALLCLOCK_S + 60,
        )
        summary.update(
            {
                "status": "CONVERGENCE_PROBE_MUSUBI_FINISHED",
                "short_run_wall_time_s": short_run.wall_time_s,
                "musubi_returncode": short_run.returncode,
            }
        )
        write_json(summary_path, summary)
        return _finalize_convergence_probe(root, run_root)
    except Exception as error:
        summary = _read_json_object(summary_path)
        summary["status"] = (
            error.status if isinstance(error, FlowError) else "CFD_FLOW_DIAGNOSTIC_INTERNAL_ERROR"
        )
        summary["failure"] = str(error)
        write_json(summary_path, summary)
        raise


def _probe_revision(run_root: Path) -> str | None:
    reference = run_root / "input_reference.json"
    if not reference.is_file():
        return None
    return str(_read_json_object(reference).get("diagnostic_revision") or "") or None


def _copy_short_ascii_results(distribution: str, tracking: Path) -> dict[str, Any]:
    """Copy completed tiny result files from the dedicated WSL temporary tree."""

    temporary_root = Path(rf"\\wsl.localhost\{distribution}\tmp\u3d")
    copied: dict[str, list[str]] = {}
    for label in (ASCII_PRESSURE_LABEL, ASCII_VELOCITY_LABEL):
        destination = tracking / label
        destination.mkdir(parents=True, exist_ok=True)
        existing = sorted(destination.glob("*.res"))
        if existing:
            copied[label] = [str(path) for path in existing]
            continue
        source_files = sorted((temporary_root / label).glob("*.res"))
        if not source_files:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_TRACKING_COPY_FAILED",
                f"No /tmp/u3d/{label} ASCII result files",
            )
        copied[label] = []
        for source in source_files:
            target = destination / source.name
            shutil.copy2(source, target)
            copied[label].append(str(target))
    checks = {
        "pressure_ascii_present": bool(copied[ASCII_PRESSURE_LABEL]),
        "velocity_ascii_present": bool(copied[ASCII_VELOCITY_LABEL]),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "copied_files": copied,
    }


def _short_density_records(text: str) -> list[dict[str, float | int]]:
    try:
        return parse_total_density_log(text)
    except FlowError:
        return []


def _continuation_advanced_from_stdout(
    path: Path,
    *,
    start_iteration: int = START_ITERATION,
) -> bool:
    if not path.is_file():
        return False
    text = path.read_text(encoding="utf-8", errors="replace")
    return any(
        int(item["iteration"]) > start_iteration for item in _short_density_records(text)
    )


def _tracking_unavailable(reason: str) -> dict[str, Any]:
    def field(threshold: float) -> dict[str, Any]:
        return {
            "sample_count": 0,
            "last_sample_delta": reason,
            "last_20_mean_absolute_delta": reason,
            "threshold": threshold,
            "ratio_to_threshold": reason,
            "trend_ratio_late_over_early": reason,
            "trend_interpretation": reason,
        }

    return {
        "status": reason,
        "pressure": field(PRESSURE_THRESHOLD_PA),
        "velocity": field(VELOCITY_THRESHOLD_M_S),
    }


def _finalize_short_ascii_probe(project_root: Path, run_root: Path) -> dict[str, Any]:
    """Finalize the one authorized run exclusively from its persisted files."""

    root = Path(project_root).resolve()
    output_root = root / "outputs" / "cfd_flow"
    source_run = output_root / SOURCE_RUN_NAME
    tracking = run_root / "tracking"
    restart_dir = run_root / "restart"
    qc_dir = run_root / "qc"
    summary_path = qc_dir / "convergence_probe_summary.json"
    summary = _read_json_object(summary_path)
    stdout = (tracking / "musubi_stdout.log").read_text(encoding="utf-8", errors="replace")
    stderr = (tracking / "musubi_stderr.log").read_text(encoding="utf-8", errors="replace")
    runtime_text = stdout + "\n" + stderr
    returncode = int(summary["musubi_returncode"])
    runtime = _runtime_numerical_errors(
        runtime_text,
        returncode,
        wrapper_timeout=returncode == 124,
    )
    write_json(qc_dir / "runtime_error_scan.json", runtime)
    start_iteration = _parse_restart_start_iteration(stdout)
    short_density = _short_density_records(stdout)
    continuation_started = any(
        int(item["iteration"]) > START_ITERATION for item in short_density
    )
    if start_iteration != START_ITERATION:
        raise FlowError(
            "CFD_FLOW_RESTART_START_INVALID",
            f"Runtime loaded iteration {start_iteration}, expected {START_ITERATION}",
        )

    if not continuation_started:
        reason = "NOT_AVAILABLE_RESTART_DID_NOT_ADVANCE"
        summary.update(
            {
                "status": "CFD_FLOW_RESTART_CONTINUATION_STILL_BLOCKED",
                "restart_continuation_start_iteration": start_iteration,
                "restart_continuation_started": False,
                "restart_continuation_end_iteration": start_iteration,
                "additional_musubi_iterations": 0,
                "runtime_error_scan": runtime,
                "tracking": _tracking_unavailable(reason),
                "convergence_classification": runtime["failure_classification"],
                "next_recommendation": "STOP_NO_RERUN",
                "diagnostic_restart_path": "NOT_AVAILABLE_RESTART_DID_NOT_ADVANCE",
            }
        )
        write_json(summary_path, summary)
        return summary

    production_lua = (source_run / "musubi" / "musubi.lua").read_text(encoding="utf-8")
    dt_s = _extract_lua_number(production_lua, "dt")
    records = parse_tracking(tracking, dt_s)
    write_tracking_csv(qc_dir / "restart_convergence_diagnostic.csv", records)
    convergence = summarize_convergence(records)
    write_json(qc_dir / "convergence_distance_summary.json", convergence)
    end_iteration = int(records[-1]["iteration"])
    additional = end_iteration - start_iteration
    if not 0 < additional <= ADDITIONAL_ITERATIONS:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_RUN_LIMIT_VIOLATED",
            f"Unexpected additional iterations: {additional}",
        )

    original_density = parse_total_density_log(
        (source_run / "musubi" / "musubi_stdout.log").read_text(
            encoding="utf-8",
            errors="replace",
        )
    )
    combined_by_iteration = {
        int(item["iteration"]): item for item in (*original_density, *short_density)
    }
    combined_density = [combined_by_iteration[key] for key in sorted(combined_by_iteration)]
    write_total_density_csv(qc_dir / "combined_total_density.csv", combined_density)
    full_density = summarize_total_density(combined_density)
    density = {
        "status": "PASS",
        "interpretation": "DIAGNOSTIC_ONLY",
        "last_1000_iterations": full_density["windows"]["1000"],
        "last_5000_iterations": full_density["windows"]["5000"],
    }
    write_json(qc_dir / "total_density_summary.json", density)

    restart = _restart_manifest(restart_dir, end_iteration)
    last_header = sorted(restart_dir.glob("*lastHeader.lua"))[-1]
    restart_iteration = int(
        _extract_lua_number(last_header.read_text(encoding="utf-8"), "iter")
    )
    if restart_iteration != end_iteration:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_MUSUBI_FAILED",
            f"Restart iteration {restart_iteration} does not match tracking {end_iteration}",
        )
    restart["last_header"] = str(last_header)
    write_json(qc_dir / "diagnostic_restart_manifest.json", restart)

    classification, next_step = classify_convergence_probe(
        convergence,
        numerical_runtime_error=bool(runtime["numerical_runtime_error"]),
    )
    source_manifest_before = _read_json_object(qc_dir / "source_manifest_before.json")
    historical_read_only = source_manifest_before == _directory_manifest(source_run)
    if not historical_read_only:
        raise FlowError(
            "CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED",
            "The formal source run changed during the convergence probe",
        )
    restart_sha_before = str(summary["restart_sha256_before"])
    restart_sha_after = sha256_file(
        source_run / "musubi" / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    )
    if restart_sha_before != restart_sha_after:
        raise FlowError(
            "CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED",
            "Formal restart SHA changed during the convergence probe",
        )

    if runtime["wrapper_timeout_detected"]:
        final_status = "CFD_FLOW_CONVERGENCE_PROBE_TIMEOUT"
    elif returncode != 0:
        final_status = "CFD_FLOW_CONVERGENCE_PROBE_RUNTIME_FAILED"
    else:
        final_status = PROBE_STATUS
    summary.update(
        {
            "status": final_status,
            "restart_continuation_start_iteration": start_iteration,
            "restart_continuation_started": True,
            "restart_continuation_end_iteration": end_iteration,
            "additional_musubi_iterations": additional,
            "runtime_error_scan": runtime,
            "tracking": convergence,
            "total_density": density,
            "diagnostic_restart": restart,
            "diagnostic_restart_path": str(last_header),
            "diagnostic_restart_status": "DIAGNOSTIC_CONTINUATION_RESTART",
            "historical_outputs_read_only": True,
            "restart_sha256_after": restart_sha_after,
            "restart_sha_unchanged": True,
            "official_steady_termination": bool(
                re.search(r"(?i)steady[- ]state\s+convergence", runtime_text)
            ),
            "convergence_classification": classification,
            "next_recommendation": next_step,
            "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
            "field_metrics": {
                "qin": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "qout": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "mass_conservation": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "actual_mach": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "spatial_outlet_pressure": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
            },
            "figures": [],
        }
    )
    summary.pop("failure", None)
    write_json(summary_path, summary)
    return summary


def _remove_owned_wsl_temp(
    *,
    distribution: str,
    workdir: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> None:
    if ASCII_TMP_ROOT_WSL != "/tmp/u3d":
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
            "Refusing to remove an unexpected temporary directory",
        )
    result = run_wsl_tool(
        distribution=distribution,
        workdir=workdir,
        command=["rm", "-rf", "--", ASCII_TMP_ROOT_WSL],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        timeout_s=30,
    )
    if result.returncode != 0:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_TEMP_FAILED",
            f"Could not remove {ASCII_TMP_ROOT_WSL}",
        )


def _copy_cleanup_finalize_probe(project_root: Path, run_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    tracking = run_root / "tracking"
    qc_dir = run_root / "qc"
    copy_qc = _copy_short_ascii_results(config.apes.wsl_distribution, tracking)
    write_json(qc_dir / "ascii_copy_qc.json", copy_qc)
    _remove_owned_wsl_temp(
        distribution=config.apes.wsl_distribution,
        workdir=run_root,
        stdout_path=qc_dir / "tmp_cleanup_after_stdout.log",
        stderr_path=qc_dir / "tmp_cleanup_after_stderr.log",
    )
    return _finalize_short_ascii_probe(root, run_root)


def run_convergence_probe(project_root: Path) -> dict[str, Any]:
    """Run the single short-path 5k restart probe, with no automatic retry."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = root / "outputs" / "cfd_flow"
    prior_runs = sorted(output_root.glob("musubi_convergence_probe_anchor003274_*"))
    current_revision_runs = [path for path in prior_runs if _probe_revision(path) == PROBE_REVISION]
    if current_revision_runs:
        latest = current_revision_runs[-1]
        summary_path = latest / "qc" / "convergence_probe_summary.json"
        prior = _read_json_object(summary_path)
        if prior.get("status") == PROBE_STATUS:
            return prior
        if int(prior.get("short_restart_musubi_run_count", 0)) == 1:
            if _continuation_advanced_from_stdout(latest / "tracking" / "musubi_stdout.log"):
                return _copy_cleanup_finalize_probe(root, latest)
            return _finalize_short_ascii_probe(root, latest)
        raise FlowError(
            "CFD_FLOW_CONVERGENCE_PROBE_RUN_LIMIT_REACHED",
            f"Current-revision probe directory already exists: {latest}",
        )

    source_run = output_root / SOURCE_RUN_NAME
    frozen_seeder_run = output_root / FROZEN_SEEDER_NAME
    restart_qc = validate_restart(source_run, frozen_seeder_run)
    restart_header = (
        source_run / "musubi" / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    )
    restart_sha_before = sha256_file(restart_header)
    production_lua = (source_run / "musubi" / "musubi.lua").read_text(encoding="utf-8")
    environment = inspect_apes_environment(config.apes)
    if environment.status != "PASS" or not environment.binaries["musubi"]:
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Pinned Musubi is unavailable")
    if environment.mpi_ranks != 8:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
            f"Expected 8 MPI ranks, found {environment.mpi_ranks}",
        )
    path_preflight = ascii_path_preflight()
    if path_preflight["status"] != "PASS":
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_ASCII_PATH_TOO_LONG",
            f"ASCII filename preflight failed: {path_preflight}",
        )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"musubi_convergence_probe_anchor003274_{timestamp}"
    tracking = run_root / "tracking"
    restart_dir = run_root / "restart"
    qc_dir = run_root / "qc"
    for directory in (
        tracking / ASCII_PRESSURE_LABEL,
        tracking / ASCII_VELOCITY_LABEL,
        restart_dir,
        qc_dir,
    ):
        directory.mkdir(parents=True, exist_ok=False)
    write_json(qc_dir / "ascii_path_preflight.json", path_preflight)
    write_json(qc_dir / "source_manifest_before.json", _directory_manifest(source_run))

    summary: dict[str, Any] = {
        "status": "CONVERGENCE_PROBE_PREFLIGHT",
        "diagnostic_revision": PROBE_REVISION,
        "run_root": str(run_root),
        "production_code_modified": False,
        "source_run": str(source_run),
        "frozen_seeder_run": str(frozen_seeder_run),
        "seeder_run_count": 0,
        "diagnostic_harvester_run_count": 0,
        "short_restart_musubi_run_count": 0,
        "long_musubi_run_count": 0,
        "additional_musubi_iterations": 0,
        "physics_or_bc_changed": False,
        "vtk_output": False,
        "figures": [],
        "restart_validation": restart_qc,
        "restart_sha256_before": restart_sha_before,
        "ascii_path_preflight": path_preflight,
        "solver_intended_budget_s": DIAGNOSTIC_WALLCLOCK_S,
        "wrapper_hard_timeout_s": WRAPPER_HARD_TIMEOUT_S,
        "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
        "environment": asdict(environment),
    }
    summary_path = qc_dir / "convergence_probe_summary.json"
    write_json(summary_path, summary)
    write_json(
        run_root / "input_reference.json",
        {
            "status": "PASS",
            "diagnostic_revision": PROBE_REVISION,
            "source_run": str(source_run),
            "restart_header": str(restart_header),
            "restart_validation": restart_qc,
            "restart_sha256": restart_sha_before,
        },
    )

    try:
        _remove_owned_wsl_temp(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            stdout_path=qc_dir / "tmp_cleanup_before_stdout.log",
            stderr_path=qc_dir / "tmp_cleanup_before_stderr.log",
        )
        make_temp = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[
                "mkdir",
                "-p",
                ASCII_PRESSURE_FOLDER_WSL,
                ASCII_VELOCITY_FOLDER_WSL,
            ],
            stdout_path=qc_dir / "tmp_mkdir_stdout.log",
            stderr_path=qc_dir / "tmp_mkdir_stderr.log",
            timeout_s=30,
        )
        if make_temp.returncode != 0:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_TEMP_FAILED",
                "Could not create the dedicated /tmp/u3d tracking tree",
            )

        diagnostic_lua = generate_diagnostic_musubi_lua(
            production_lua,
            pressure_folder_wsl=ASCII_PRESSURE_FOLDER_WSL,
            velocity_folder_wsl=ASCII_VELOCITY_FOLDER_WSL,
            restart_folder_wsl=windows_to_wsl(
                restart_dir,
                config.apes.wsl_distribution,
            ),
            timing_file_wsl=f"{ASCII_TMP_ROOT_WSL}/timing.res",
        )
        diagnostic_lua_path = run_root / "diagnostic_musubi.lua"
        diagnostic_lua_path.write_text(diagnostic_lua, encoding="utf-8")
        lua_contract = diagnostic_lua_contract(production_lua, diagnostic_lua)
        write_json(qc_dir / "diagnostic_musubi_contract.json", lua_contract)
        if lua_contract["status"] != "PASS":
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
                "Frozen physics/BC or diagnostic time-control contract failed",
            )
        luac = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[str(environment.binaries["lua_compiler"]), "-p", "diagnostic_musubi.lua"],
            stdout_path=tracking / "luac_musubi_stdout.log",
            stderr_path=tracking / "luac_musubi_stderr.log",
            timeout_s=30,
        )
        if luac.returncode != 0:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
                "Diagnostic Musubi Lua syntax failed",
            )

        diagnostic_config_wsl = windows_to_wsl(
            diagnostic_lua_path,
            config.apes.wsl_distribution,
        )
        mpi_command = [
            "env",
            "OMPI_ALLOW_RUN_AS_ROOT=1",
            "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            str(environment.binaries["mpi_launcher"]),
            "-np",
            "8",
            str(environment.binaries["musubi"]),
            diagnostic_config_wsl,
        ]
        summary.update(
            {
                "status": "CONVERGENCE_PROBE_MUSUBI_RUNNING",
                "short_restart_musubi_run_count": 1,
                "restart_continuation_start_iteration": START_ITERATION,
                "requested_end_iteration": END_ITERATION,
                "requested_additional_iterations": ADDITIONAL_ITERATIONS,
            }
        )
        write_json(summary_path, summary)
        short_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=source_run / "musubi",
            command=mpi_command,
            stdout_path=tracking / "musubi_stdout.log",
            stderr_path=tracking / "musubi_stderr.log",
            timeout_s=WRAPPER_HARD_TIMEOUT_S,
        )
        summary.update(
            {
                "status": "CONVERGENCE_PROBE_MUSUBI_FINISHED",
                "short_run_wall_time_s": short_run.wall_time_s,
                "musubi_returncode": short_run.returncode,
            }
        )
        write_json(summary_path, summary)
        if (
            short_run.returncode != 0
            and not _continuation_advanced_from_stdout(short_run.stdout_path)
        ):
            return _finalize_short_ascii_probe(root, run_root)
        return _copy_cleanup_finalize_probe(root, run_root)
    except Exception as error:
        latest_summary = _read_json_object(summary_path)
        latest_summary["status"] = (
            error.status if isinstance(error, FlowError) else "CFD_FLOW_DIAGNOSTIC_INTERNAL_ERROR"
        )
        latest_summary["failure"] = str(error)
        write_json(summary_path, latest_summary)
        raise


def _official_equivalent_run_revision(run_root: Path) -> str | None:
    reference = run_root / "input_reference.json"
    if not reference.is_file():
        return None
    return str(_read_json_object(reference).get("diagnostic_revision") or "") or None


def _finalize_official_equivalent_probe(
    project_root: Path,
    run_root: Path,
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    output_root = root / "outputs" / "cfd_flow"
    previous_probe = output_root / PREVIOUS_PROBE_NAME
    formal_source = output_root / SOURCE_RUN_NAME
    tracking = run_root / "tracking"
    restart_dir = run_root / "restart"
    qc_dir = run_root / "qc"
    summary_path = qc_dir / "official_equivalent_probe_summary.json"
    summary = _read_json_object(summary_path)
    stdout = (tracking / "musubi_stdout.log").read_text(encoding="utf-8", errors="replace")
    stderr = (tracking / "musubi_stderr.log").read_text(encoding="utf-8", errors="replace")
    runtime_text = stdout + "\n" + stderr
    returncode = int(summary["musubi_returncode"])
    runtime = _runtime_numerical_errors(
        runtime_text,
        returncode,
        wrapper_timeout=returncode == 124,
    )
    write_json(qc_dir / "runtime_error_scan.json", runtime)
    start_iteration = _parse_restart_start_iteration(stdout)
    short_density = _short_density_records(stdout)
    continuation_started = any(
        int(item["iteration"]) > OFFICIAL_EQUIVALENT_START_ITERATION
        for item in short_density
    )
    if start_iteration != OFFICIAL_EQUIVALENT_START_ITERATION:
        raise FlowError(
            "CFD_FLOW_RESTART_START_INVALID",
            f"Runtime loaded iteration {start_iteration}, expected "
            f"{OFFICIAL_EQUIVALENT_START_ITERATION}",
        )
    if not continuation_started:
        summary.update(
            {
                "status": "CFD_FLOW_RESTART_CONTINUATION_STILL_BLOCKED",
                "restart_continuation_started": False,
                "restart_continuation_start_iteration": start_iteration,
                "restart_continuation_end_iteration": start_iteration,
                "additional_musubi_iterations": 0,
                "runtime_error_scan": runtime,
                "next_recommendation": "STOP_NO_RERUN",
            }
        )
        write_json(summary_path, summary)
        return summary

    production_lua = (formal_source / "musubi" / "musubi.lua").read_text(encoding="utf-8")
    dt_s = _extract_lua_number(production_lua, "dt")
    previous_records = parse_tracking(previous_probe / "tracking", dt_s)
    current_records = parse_tracking(tracking, dt_s)
    write_tracking_csv(qc_dir / "new_convergence_tracking.csv", current_records)
    combined, continuity = combine_continuous_tracking(previous_records, current_records)
    write_json(qc_dir / "tracking_continuity_qc.json", continuity)
    write_combined_tracking_csv(qc_dir / "combined_convergence_tracking.csv", combined)
    residuals = official_equivalent_residuals(combined)
    write_official_equivalent_residuals_csv(
        qc_dir / "official_equivalent_residuals.csv",
        residuals,
    )
    equivalent = summarize_official_equivalent_residuals(combined, residuals)
    write_json(qc_dir / "official_equivalent_residual_summary.json", equivalent)
    end_iteration = int(current_records[-1]["iteration"])
    additional = end_iteration - start_iteration
    if not 0 < additional <= OFFICIAL_EQUIVALENT_ADDITIONAL_ITERATIONS:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_RUN_LIMIT_VIOLATED",
            f"Unexpected additional iterations: {additional}",
        )

    pressure_adjacent = np.abs(
        np.asarray(
            [item["pressure_adjacent_delta"] for item in combined[1:]],
            dtype=float,
        )
    )
    velocity_adjacent = np.abs(
        np.asarray(
            [item["velocity_adjacent_delta"] for item in combined[1:]],
            dtype=float,
        )
    )
    adjacent = {
        "status": "SHORT_TIMESCALE_DIAGNOSTIC_ONLY",
        "cannot_set_final_steady_classification": True,
        "pressure_final_adjacent_delta": float(pressure_adjacent[-1]),
        "pressure_last20_mean_absolute_delta": float(np.mean(pressure_adjacent[-20:])),
        "pressure_adjacent_delta_ratio": (
            float(np.mean(pressure_adjacent[-20:])) / PRESSURE_THRESHOLD_PA
        ),
        "velocity_final_adjacent_delta": float(velocity_adjacent[-1]),
        "velocity_last20_mean_absolute_delta": float(np.mean(velocity_adjacent[-20:])),
        "velocity_adjacent_delta_ratio": (
            float(np.mean(velocity_adjacent[-20:])) / VELOCITY_THRESHOLD_M_S
        ),
    }
    write_json(qc_dir / "adjacent_delta_diagnostic.json", adjacent)

    original_density = parse_total_density_log(
        (formal_source / "musubi" / "musubi_stdout.log").read_text(
            encoding="utf-8",
            errors="replace",
        )
    )
    previous_density = parse_total_density_log(
        (previous_probe / "tracking" / "musubi_stdout.log").read_text(
            encoding="utf-8",
            errors="replace",
        )
    )
    density_by_iteration = {
        int(item["iteration"]): item
        for item in (*original_density, *previous_density, *short_density)
    }
    combined_density = [density_by_iteration[key] for key in sorted(density_by_iteration)]
    write_total_density_csv(qc_dir / "combined_total_density.csv", combined_density)
    full_density = summarize_total_density(combined_density)
    density = {
        "status": "DIAGNOSTIC_ONLY",
        "last_1000_iterations": full_density["windows"]["1000"],
        "last_5000_iterations": full_density["windows"]["5000"],
        "last_10000_iterations": full_density["windows"]["10000"],
    }
    write_json(qc_dir / "total_density_summary.json", density)

    restart = _restart_manifest(restart_dir, end_iteration)
    last_header = sorted(restart_dir.glob("*lastHeader.lua"))[-1]
    saved_iteration = int(
        _extract_lua_number(last_header.read_text(encoding="utf-8"), "iter")
    )
    if saved_iteration != end_iteration:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_MUSUBI_FAILED",
            f"Restart iteration {saved_iteration} does not match tracking {end_iteration}",
        )
    restart["last_header"] = str(last_header)
    write_json(qc_dir / "diagnostic_restart_manifest.json", restart)

    classification, next_step = classify_official_equivalent(
        equivalent,
        numerical_runtime_error=bool(runtime["numerical_runtime_error"]),
    )
    source_manifest_before = _read_json_object(qc_dir / "source_manifest_before.json")
    source_read_only = source_manifest_before == _directory_manifest(previous_probe)
    source_restart = previous_probe / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    restart_sha_after = sha256_file(source_restart)
    restart_sha_unchanged = restart_sha_after == str(summary["source_restart_sha256_before"])
    if not source_read_only or not restart_sha_unchanged:
        raise FlowError(
            "CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED",
            "Previous probe or its restart changed during continuation",
        )
    if runtime["wrapper_timeout_detected"]:
        final_status = "CFD_FLOW_OFFICIAL_EQUIVALENT_CONVERGENCE_PROBE_TIMEOUT"
    elif returncode != 0:
        final_status = "CFD_FLOW_OFFICIAL_EQUIVALENT_CONVERGENCE_PROBE_RUNTIME_FAILED"
    else:
        final_status = OFFICIAL_EQUIVALENT_STATUS
    wall_time = float(summary["short_run_wall_time_s"])
    summary.update(
        {
            "status": final_status,
            "restart_continuation_started": True,
            "restart_continuation_start_iteration": start_iteration,
            "restart_continuation_end_iteration": end_iteration,
            "additional_musubi_iterations": additional,
            "iterations_per_second": additional / wall_time,
            "runtime_error_scan": runtime,
            "tracking_continuity": continuity,
            "official_equivalent": equivalent,
            "adjacent_delta_diagnostic": adjacent,
            "total_density": density,
            "diagnostic_restart": restart,
            "diagnostic_restart_path": str(last_header),
            "diagnostic_restart_status": "DIAGNOSTIC_CONTINUATION_RESTART",
            "source_outputs_read_only": True,
            "source_restart_sha256_after": restart_sha_after,
            "source_restart_sha_unchanged": True,
            "old_heuristic_classification": {
                "classification": "NEAR_STEADY",
                "status": "HEURISTIC_ADJACENT_DELTA_CLASSIFICATION",
                "invalid_as_formal_steady_conclusion": True,
            },
            "convergence_classification": classification,
            "next_recommendation": next_step,
            "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
            "field_metrics": {
                "qin": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "qout": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "mass_conservation": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "actual_mach": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "spatial_pressure": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
            },
            "figures": [],
        }
    )
    summary.pop("failure", None)
    write_json(summary_path, summary)
    return summary


def _copy_cleanup_finalize_official_equivalent(
    project_root: Path,
    run_root: Path,
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    tracking = run_root / "tracking"
    qc_dir = run_root / "qc"
    copy_qc = _copy_short_ascii_results(config.apes.wsl_distribution, tracking)
    write_json(qc_dir / "ascii_copy_qc.json", copy_qc)
    _remove_owned_wsl_temp(
        distribution=config.apes.wsl_distribution,
        workdir=run_root,
        stdout_path=qc_dir / "tmp_cleanup_after_stdout.log",
        stderr_path=qc_dir / "tmp_cleanup_after_stderr.log",
    )
    return _finalize_official_equivalent_probe(root, run_root)


def run_official_equivalent_probe(project_root: Path) -> dict[str, Any]:
    """Run one 6k continuation and evaluate pinned nvals=100 semantics offline."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = root / "outputs" / "cfd_flow"
    existing = sorted(output_root.glob("musubi_official_equivalent_probe_anchor003274_*"))
    revision_runs = [
        path
        for path in existing
        if _official_equivalent_run_revision(path) == OFFICIAL_EQUIVALENT_REVISION
    ]
    if revision_runs:
        latest = revision_runs[-1]
        summary_path = latest / "qc" / "official_equivalent_probe_summary.json"
        prior = _read_json_object(summary_path)
        if prior.get("status") == OFFICIAL_EQUIVALENT_STATUS:
            return prior
        if int(prior.get("short_restart_musubi_run_count", 0)) == 1:
            if _continuation_advanced_from_stdout(
                latest / "tracking" / "musubi_stdout.log",
                start_iteration=OFFICIAL_EQUIVALENT_START_ITERATION,
            ):
                return _copy_cleanup_finalize_official_equivalent(root, latest)
            return _finalize_official_equivalent_probe(root, latest)
        raise FlowError(
            "CFD_FLOW_OFFICIAL_EQUIVALENT_RUN_LIMIT_REACHED",
            f"Current-revision probe directory already exists: {latest}",
        )

    previous_probe = output_root / PREVIOUS_PROBE_NAME
    formal_source = output_root / SOURCE_RUN_NAME
    frozen_seeder_run = output_root / FROZEN_SEEDER_NAME
    source_restart = previous_probe / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    restart_qc = validate_continuation_restart(
        source_restart,
        solver_workdir=formal_source / "musubi",
        frozen_seeder_run=frozen_seeder_run,
        expected_iteration=OFFICIAL_EQUIVALENT_START_ITERATION,
    )
    source_restart_sha = sha256_file(source_restart)
    production_lua = (formal_source / "musubi" / "musubi.lua").read_text(encoding="utf-8")
    environment = inspect_apes_environment(config.apes)
    if environment.status != "PASS" or not environment.binaries["musubi"]:
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Pinned Musubi is unavailable")
    if environment.mpi_ranks != 8:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
            f"Expected 8 MPI ranks, found {environment.mpi_ranks}",
        )
    path_preflight = ascii_path_preflight()
    if path_preflight["status"] != "PASS":
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_ASCII_PATH_TOO_LONG",
            f"ASCII filename preflight failed: {path_preflight}",
        )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"musubi_official_equivalent_probe_anchor003274_{timestamp}"
    tracking = run_root / "tracking"
    restart_dir = run_root / "restart"
    qc_dir = run_root / "qc"
    for directory in (
        tracking / ASCII_PRESSURE_LABEL,
        tracking / ASCII_VELOCITY_LABEL,
        restart_dir,
        qc_dir,
    ):
        directory.mkdir(parents=True, exist_ok=False)
    write_json(qc_dir / "ascii_path_preflight.json", path_preflight)
    write_json(qc_dir / "source_manifest_before.json", _directory_manifest(previous_probe))
    semantics = {
        "status": "PASS",
        "equivalence_basis": "PINNED_TREELM_NORM_AVERAGE_NVALS_100_ABSOLUTE",
        "norm": "average",
        "nvals": 100,
        "absolute": True,
        "spatial_reduction": "DOMAIN_AVERAGE",
        "history_window_excludes_current": True,
        "musubi_internal_convergence_history_persists_across_restart": False,
        "source": "pinned tem/source/tem_convergence_module.f90",
    }
    write_json(qc_dir / "pinned_convergence_semantics.json", semantics)
    summary: dict[str, Any] = {
        "status": "OFFICIAL_EQUIVALENT_PROBE_PREFLIGHT",
        "diagnostic_revision": OFFICIAL_EQUIVALENT_REVISION,
        "run_root": str(run_root),
        "production_code_modified": False,
        "source_run": str(previous_probe),
        "solver_workdir": str(formal_source / "musubi"),
        "frozen_seeder_run": str(frozen_seeder_run),
        "seeder_run_count": 0,
        "diagnostic_harvester_run_count": 0,
        "short_restart_musubi_run_count": 0,
        "long_musubi_run_count": 0,
        "additional_musubi_iterations": 0,
        "physics_or_bc_changed": False,
        "vtk_output": False,
        "figures": [],
        "restart_validation": restart_qc,
        "source_restart_sha256_before": source_restart_sha,
        "ascii_path_preflight": path_preflight,
        "pinned_convergence_semantics": semantics,
        "solver_expected_budget_s": OFFICIAL_EQUIVALENT_SOLVER_BUDGET_S,
        "wrapper_hard_timeout_s": OFFICIAL_EQUIVALENT_WRAPPER_TIMEOUT_S,
        "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
        "environment": asdict(environment),
    }
    summary_path = qc_dir / "official_equivalent_probe_summary.json"
    write_json(summary_path, summary)
    write_json(
        run_root / "input_reference.json",
        {
            "status": "PASS",
            "diagnostic_revision": OFFICIAL_EQUIVALENT_REVISION,
            "previous_probe": str(previous_probe),
            "source_restart": str(source_restart),
            "source_restart_sha256": source_restart_sha,
            "restart_validation": restart_qc,
        },
    )

    try:
        _remove_owned_wsl_temp(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            stdout_path=qc_dir / "tmp_cleanup_before_stdout.log",
            stderr_path=qc_dir / "tmp_cleanup_before_stderr.log",
        )
        make_temp = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[
                "mkdir",
                "-p",
                ASCII_PRESSURE_FOLDER_WSL,
                ASCII_VELOCITY_FOLDER_WSL,
            ],
            stdout_path=qc_dir / "tmp_mkdir_stdout.log",
            stderr_path=qc_dir / "tmp_mkdir_stderr.log",
            timeout_s=30,
        )
        if make_temp.returncode != 0:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_TEMP_FAILED",
                "Could not create /tmp/u3d tracking folders",
            )
        restart_read_wsl = windows_to_wsl(source_restart, config.apes.wsl_distribution)
        diagnostic_lua = generate_diagnostic_musubi_lua(
            production_lua,
            pressure_folder_wsl=ASCII_PRESSURE_FOLDER_WSL,
            velocity_folder_wsl=ASCII_VELOCITY_FOLDER_WSL,
            restart_folder_wsl=windows_to_wsl(
                restart_dir,
                config.apes.wsl_distribution,
            ),
            timing_file_wsl=f"{ASCII_TMP_ROOT_WSL}/timing.res",
            start_iteration=OFFICIAL_EQUIVALENT_START_ITERATION,
            end_iteration=OFFICIAL_EQUIVALENT_END_ITERATION,
            restart_read=restart_read_wsl,
        )
        diagnostic_lua_path = run_root / "diagnostic_musubi.lua"
        diagnostic_lua_path.write_text(diagnostic_lua, encoding="utf-8")
        lua_contract = diagnostic_lua_contract(
            production_lua,
            diagnostic_lua,
            start_iteration=OFFICIAL_EQUIVALENT_START_ITERATION,
            end_iteration=OFFICIAL_EQUIVALENT_END_ITERATION,
            restart_read=restart_read_wsl,
        )
        write_json(qc_dir / "diagnostic_musubi_contract.json", lua_contract)
        if lua_contract["status"] != "PASS":
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
                "Frozen physics/BC or pinned convergence contract failed",
            )
        luac = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[str(environment.binaries["lua_compiler"]), "-p", "diagnostic_musubi.lua"],
            stdout_path=tracking / "luac_musubi_stdout.log",
            stderr_path=tracking / "luac_musubi_stderr.log",
            timeout_s=30,
        )
        if luac.returncode != 0:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
                "Diagnostic Musubi Lua syntax failed",
            )
        diagnostic_config_wsl = windows_to_wsl(
            diagnostic_lua_path,
            config.apes.wsl_distribution,
        )
        mpi_command = [
            "env",
            "OMPI_ALLOW_RUN_AS_ROOT=1",
            "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            str(environment.binaries["mpi_launcher"]),
            "-np",
            "8",
            str(environment.binaries["musubi"]),
            diagnostic_config_wsl,
        ]
        summary.update(
            {
                "status": "OFFICIAL_EQUIVALENT_PROBE_MUSUBI_RUNNING",
                "short_restart_musubi_run_count": 1,
                "restart_continuation_start_iteration": OFFICIAL_EQUIVALENT_START_ITERATION,
                "requested_end_iteration": OFFICIAL_EQUIVALENT_END_ITERATION,
                "requested_additional_iterations": OFFICIAL_EQUIVALENT_ADDITIONAL_ITERATIONS,
            }
        )
        write_json(summary_path, summary)
        short_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=formal_source / "musubi",
            command=mpi_command,
            stdout_path=tracking / "musubi_stdout.log",
            stderr_path=tracking / "musubi_stderr.log",
            timeout_s=OFFICIAL_EQUIVALENT_WRAPPER_TIMEOUT_S,
        )
        summary.update(
            {
                "status": "OFFICIAL_EQUIVALENT_PROBE_MUSUBI_FINISHED",
                "short_run_wall_time_s": short_run.wall_time_s,
                "musubi_returncode": short_run.returncode,
            }
        )
        write_json(summary_path, summary)
        if (
            short_run.returncode != 0
            and not _continuation_advanced_from_stdout(
                short_run.stdout_path,
                start_iteration=OFFICIAL_EQUIVALENT_START_ITERATION,
            )
        ):
            return _finalize_official_equivalent_probe(root, run_root)
        return _copy_cleanup_finalize_official_equivalent(root, run_root)
    except Exception as error:
        latest_summary = _read_json_object(summary_path)
        latest_summary["status"] = (
            error.status if isinstance(error, FlowError) else "CFD_FLOW_DIAGNOSTIC_INTERNAL_ERROR"
        )
        latest_summary["failure"] = str(error)
        write_json(summary_path, latest_summary)
        raise


def _project_steady_run_revision(run_root: Path) -> str | None:
    reference = run_root / "input_reference.json"
    if not reference.is_file():
        return None
    return str(_read_json_object(reference).get("diagnostic_revision") or "") or None


def _finalize_project_steady_confirmation(
    project_root: Path,
    run_root: Path,
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    output_root = root / "outputs" / "cfd_flow"
    source_run = output_root / PROJECT_STEADY_SOURCE_RUN_NAME
    formal_source = output_root / SOURCE_RUN_NAME
    tracking = run_root / "tracking"
    restart_dir = run_root / "restart"
    qc_dir = run_root / "qc"
    summary_path = qc_dir / "project_steady_confirmation.json"
    summary = _read_json_object(summary_path)
    criterion = _read_json_object(qc_dir / "convergence_criterion_review.json")
    stdout = (tracking / "musubi_stdout.log").read_text(encoding="utf-8", errors="replace")
    stderr = (tracking / "musubi_stderr.log").read_text(encoding="utf-8", errors="replace")
    runtime_text = stdout + "\n" + stderr
    returncode = int(summary["musubi_returncode"])
    runtime = _runtime_numerical_errors(
        runtime_text,
        returncode,
        wrapper_timeout=returncode == 124,
    )
    steady_evidence = parse_official_steady_termination(runtime_text)
    write_json(qc_dir / "runtime_error_scan.json", runtime)
    write_json(qc_dir / "official_steady_termination_evidence.json", steady_evidence)
    start_iteration = _parse_restart_start_iteration(stdout)
    short_density = _short_density_records(stdout)
    continuation_started = any(
        int(item["iteration"]) > PROJECT_STEADY_START_ITERATION
        for item in short_density
    )
    if start_iteration != PROJECT_STEADY_START_ITERATION:
        raise FlowError(
            "CFD_FLOW_RESTART_START_INVALID",
            f"Runtime loaded iteration {start_iteration}, expected {PROJECT_STEADY_START_ITERATION}",
        )
    if not continuation_started:
        summary.update(
            {
                "status": "CFD_FLOW_RESTART_CONTINUATION_STILL_BLOCKED",
                "restart_continuation_started": False,
                "restart_continuation_start_iteration": start_iteration,
                "restart_continuation_end_iteration": start_iteration,
                "additional_musubi_iterations": 0,
                "runtime_error_scan": runtime,
                "official_steady_termination": steady_evidence,
                "next_recommendation": "STOP_NO_RERUN",
            }
        )
        write_json(summary_path, summary)
        return summary

    production_lua = (formal_source / "musubi" / "musubi.lua").read_text(encoding="utf-8")
    dt_s = _extract_lua_number(production_lua, "dt")
    records = parse_tracking(tracking, dt_s)
    write_tracking_csv(qc_dir / "project_steady_tracking.csv", records)
    residuals = official_equivalent_residuals(records)
    write_official_equivalent_residuals_csv(
        qc_dir / "project_steady_residuals.csv",
        residuals,
    )
    residual_summary = summarize_official_equivalent_residuals(records, residuals)
    project_pressure_threshold = float(criterion["derived_thresholds"]["pressure_pa"])
    project_velocity_threshold = float(criterion["derived_thresholds"]["velocity_m_s"])
    final_pressure_residual = float(residuals[-1]["pressure_residual_pa"])
    final_velocity_residual = float(residuals[-1]["velocity_residual_m_s"])
    final_project_metrics = {
        "pressure_residual_pa": final_pressure_residual,
        "pressure_threshold_pa": project_pressure_threshold,
        "pressure_ratio_to_project_threshold": (
            final_pressure_residual / project_pressure_threshold
        ),
        "velocity_residual_m_s": final_velocity_residual,
        "velocity_threshold_m_s": project_velocity_threshold,
        "velocity_ratio_to_project_threshold": (
            final_velocity_residual / project_velocity_threshold
        ),
        "offline_same_process_equivalent": True,
    }
    write_json(qc_dir / "final_project_criterion_metrics.json", final_project_metrics)
    end_iteration = int(records[-1]["iteration"])
    additional = end_iteration - start_iteration
    if not 0 < additional <= PROJECT_STEADY_ADDITIONAL_ITERATIONS:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_RUN_LIMIT_VIOLATED",
            f"Unexpected additional iterations: {additional}",
        )
    confirmation_iteration = steady_evidence["confirmation_iteration"]
    if confirmation_iteration is not None and int(confirmation_iteration) != end_iteration:
        raise FlowError(
            "CFD_FLOW_PROJECT_STEADY_EVIDENCE_INVALID",
            "Official steady iteration does not match final tracking iteration",
        )

    restart = _restart_manifest(restart_dir, end_iteration)
    last_header = sorted(restart_dir.glob("*lastHeader.lua"))[-1]
    saved_iteration = int(
        _extract_lua_number(last_header.read_text(encoding="utf-8"), "iter")
    )
    if saved_iteration != end_iteration:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_MUSUBI_FAILED",
            f"Restart iteration {saved_iteration} does not match tracking {end_iteration}",
        )
    official_pass = bool(steady_evidence["official_steady_termination"])
    if official_pass:
        restart.update(
            {
                "status": "FROZEN_PROJECT_STEADY_RESTART",
                "not_final_cfd_solution": True,
                "field_export_pending": True,
                "mass_conservation_field_qc_pending": True,
                "grid_convergence_pending": True,
            }
        )
        restart_manifest_path = qc_dir / "steady_restart_manifest.json"
    else:
        restart_manifest_path = qc_dir / "diagnostic_restart_manifest.json"
    restart["last_header"] = str(last_header)
    write_json(restart_manifest_path, restart)

    source_manifest_before = _read_json_object(qc_dir / "source_manifest_before.json")
    source_read_only = source_manifest_before == _directory_manifest(source_run)
    source_restart = source_run / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    source_restart_sha_after = sha256_file(source_restart)
    source_restart_unchanged = (
        source_restart_sha_after == str(summary["source_restart_sha256_before"])
    )
    if not source_read_only or not source_restart_unchanged:
        raise FlowError(
            "CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED",
            "Source official-equivalent run changed during steady confirmation",
        )
    if runtime["numerical_runtime_error"]:
        final_status = "CFD_FLOW_PROJECT_STEADY_NUMERICAL_PROBLEM"
        next_step = "REVIEW_NUMERICAL_OR_BOUNDARY_CONDITION"
    elif runtime["wrapper_timeout_detected"]:
        final_status = "CFD_FLOW_PROJECT_STEADY_CONFIRMATION_TIMEOUT"
        next_step = "REVIEW_RESIDUAL_TREND_BEFORE_ANY_FURTHER_SOLVER_RUN"
    elif official_pass:
        final_status = PROJECT_STEADY_CONFIRMED_STATUS
        next_step = "FIX_OR_BYPASS_MUS_HARVESTING_AND_EXPORT_FROZEN_STEADY_FLOW_FIELD"
    else:
        final_status = PROJECT_STEADY_NOT_REACHED_STATUS
        next_step = "REVIEW_RESIDUAL_TREND_BEFORE_ANY_FURTHER_SOLVER_RUN"
    wall_time = float(summary["short_run_wall_time_s"])
    summary.update(
        {
            "status": final_status,
            "restart_continuation_started": True,
            "restart_continuation_start_iteration": start_iteration,
            "restart_continuation_end_iteration": end_iteration,
            "additional_musubi_iterations": additional,
            "iterations_per_second": additional / wall_time,
            "runtime_error_scan": runtime,
            "official_steady_termination": steady_evidence,
            "criterion_review": criterion,
            "same_process_official_equivalent_residual_summary": residual_summary,
            "final_project_criterion_metrics": final_project_metrics,
            "frozen_project_steady_restart": official_pass,
            "restart": restart,
            "restart_path": str(last_header),
            "source_outputs_read_only": True,
            "source_restart_sha256_after": source_restart_sha_after,
            "source_restart_sha_unchanged": True,
            "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
            "grid_convergence": "NOT_RUN",
            "field_metrics": {
                "qin": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "qout": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "mass_conservation": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "actual_mach": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
                "spatial_pressure": "DEFERRED_UNTIL_FIELD_EXPORT_AVAILABLE",
            },
            "next_recommendation": next_step,
        }
    )
    summary.pop("failure", None)
    write_json(summary_path, summary)
    return summary


def _copy_cleanup_finalize_project_steady(
    project_root: Path,
    run_root: Path,
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    tracking = run_root / "tracking"
    qc_dir = run_root / "qc"
    copy_qc = _copy_short_ascii_results(config.apes.wsl_distribution, tracking)
    write_json(qc_dir / "ascii_copy_qc.json", copy_qc)
    _remove_owned_wsl_temp(
        distribution=config.apes.wsl_distribution,
        workdir=run_root,
        stdout_path=qc_dir / "tmp_cleanup_after_stdout.log",
        stderr_path=qc_dir / "tmp_cleanup_after_stderr.log",
    )
    return _finalize_project_steady_confirmation(root, run_root)


def run_project_steady_confirmation(project_root: Path) -> dict[str, Any]:
    """Run the one pre-registered project-level 0.1% steady confirmation."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = root / "outputs" / "cfd_flow"
    existing = sorted(output_root.glob("musubi_project_steady_confirmation_anchor003274_*"))
    revision_runs = [
        path for path in existing if _project_steady_run_revision(path) == PROJECT_STEADY_REVISION
    ]
    if revision_runs:
        latest = revision_runs[-1]
        summary_path = latest / "qc" / "project_steady_confirmation.json"
        prior = _read_json_object(summary_path)
        if prior.get("status") in (
            PROJECT_STEADY_CONFIRMED_STATUS,
            PROJECT_STEADY_NOT_REACHED_STATUS,
        ):
            return prior
        if int(prior.get("short_restart_musubi_run_count", 0)) == 1:
            if _continuation_advanced_from_stdout(
                latest / "tracking" / "musubi_stdout.log",
                start_iteration=PROJECT_STEADY_START_ITERATION,
            ):
                return _copy_cleanup_finalize_project_steady(root, latest)
            return _finalize_project_steady_confirmation(root, latest)
        raise FlowError(
            "CFD_FLOW_PROJECT_STEADY_RUN_LIMIT_REACHED",
            f"Current-revision confirmation directory already exists: {latest}",
        )

    source_run = output_root / PROJECT_STEADY_SOURCE_RUN_NAME
    formal_source = output_root / SOURCE_RUN_NAME
    frozen_seeder_run = output_root / FROZEN_SEEDER_NAME
    source_restart = source_run / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    restart_qc = validate_continuation_restart(
        source_restart,
        solver_workdir=formal_source / "musubi",
        frozen_seeder_run=frozen_seeder_run,
        expected_iteration=PROJECT_STEADY_START_ITERATION,
    )
    source_restart_sha = sha256_file(source_restart)
    scaling_path = frozen_seeder_run / "qc" / "lbm_scaling_qc.json"
    scaling_qc = _read_json_object(scaling_path)
    prior_residual = _read_json_object(
        source_run / "qc" / "official_equivalent_residual_summary.json"
    )
    criterion = build_project_criterion_review(
        outlet_gauge_pressures_pa=scaling_qc["outlet_gauge_pressures_pa"],
        inlet_target_mean_velocity_m_s=float(scaling_qc["velocity_mean_m_s"]),
        current_pressure_residual_pa=float(prior_residual["pressure"]["final_residual"]),
        current_velocity_residual_m_s=float(prior_residual["velocity"]["final_residual"]),
    )
    production_lua = (formal_source / "musubi" / "musubi.lua").read_text(encoding="utf-8")
    environment = inspect_apes_environment(config.apes)
    if environment.status != "PASS" or not environment.binaries["musubi"]:
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Pinned Musubi is unavailable")
    if environment.mpi_ranks != 8:
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
            f"Expected 8 MPI ranks, found {environment.mpi_ranks}",
        )
    path_preflight = ascii_path_preflight()
    if path_preflight["status"] != "PASS":
        raise FlowError(
            "CFD_FLOW_DIAGNOSTIC_ASCII_PATH_TOO_LONG",
            f"ASCII filename preflight failed: {path_preflight}",
        )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"musubi_project_steady_confirmation_anchor003274_{timestamp}"
    tracking = run_root / "tracking"
    restart_dir = run_root / "restart"
    qc_dir = run_root / "qc"
    for directory in (
        tracking / ASCII_PRESSURE_LABEL,
        tracking / ASCII_VELOCITY_LABEL,
        restart_dir,
        qc_dir,
    ):
        directory.mkdir(parents=True, exist_ok=False)
    write_json(qc_dir / "convergence_criterion_review.json", criterion)
    write_json(qc_dir / "ascii_path_preflight.json", path_preflight)
    write_json(qc_dir / "source_manifest_before.json", _directory_manifest(source_run))
    summary: dict[str, Any] = {
        "status": "PROJECT_STEADY_CONFIRMATION_PREFLIGHT",
        "diagnostic_revision": PROJECT_STEADY_REVISION,
        "run_root": str(run_root),
        "production_code_modified": False,
        "source_run": str(source_run),
        "solver_workdir": str(formal_source / "musubi"),
        "frozen_seeder_run": str(frozen_seeder_run),
        "seeder_run_count": 0,
        "diagnostic_harvester_run_count": 0,
        "short_restart_musubi_run_count": 0,
        "long_musubi_run_count": 0,
        "additional_musubi_iterations": 0,
        "physics_or_bc_changed": False,
        "production_config_modified": False,
        "tolerance_sweep": False,
        "vtk_output": False,
        "figures": [],
        "restart_validation": restart_qc,
        "source_restart_sha256_before": source_restart_sha,
        "criterion_review": criterion,
        "ascii_path_preflight": path_preflight,
        "solver_expected_budget_s": PROJECT_STEADY_EXPECTED_BUDGET_S,
        "wrapper_hard_timeout_s": PROJECT_STEADY_WRAPPER_TIMEOUT_S,
        "harvester_issue": "MUS_HARVESTING_SIGSEGV_DEFERRED",
        "grid_convergence": "NOT_RUN",
        "environment": asdict(environment),
    }
    summary_path = qc_dir / "project_steady_confirmation.json"
    write_json(summary_path, summary)
    write_json(
        run_root / "input_reference.json",
        {
            "status": "PASS",
            "diagnostic_revision": PROJECT_STEADY_REVISION,
            "source_run": str(source_run),
            "source_restart": str(source_restart),
            "source_restart_sha256": source_restart_sha,
            "restart_validation": restart_qc,
            "lattice_scaling_qc": str(scaling_path),
        },
    )

    try:
        _remove_owned_wsl_temp(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            stdout_path=qc_dir / "tmp_cleanup_before_stdout.log",
            stderr_path=qc_dir / "tmp_cleanup_before_stderr.log",
        )
        make_temp = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[
                "mkdir",
                "-p",
                ASCII_PRESSURE_FOLDER_WSL,
                ASCII_VELOCITY_FOLDER_WSL,
            ],
            stdout_path=qc_dir / "tmp_mkdir_stdout.log",
            stderr_path=qc_dir / "tmp_mkdir_stderr.log",
            timeout_s=30,
        )
        if make_temp.returncode != 0:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_TEMP_FAILED",
                "Could not create /tmp/u3d tracking folders",
            )
        restart_read_wsl = windows_to_wsl(source_restart, config.apes.wsl_distribution)
        pressure_threshold = float(criterion["derived_thresholds"]["pressure_pa"])
        velocity_threshold = float(criterion["derived_thresholds"]["velocity_m_s"])
        diagnostic_lua = generate_diagnostic_musubi_lua(
            production_lua,
            pressure_folder_wsl=ASCII_PRESSURE_FOLDER_WSL,
            velocity_folder_wsl=ASCII_VELOCITY_FOLDER_WSL,
            restart_folder_wsl=windows_to_wsl(
                restart_dir,
                config.apes.wsl_distribution,
            ),
            timing_file_wsl=f"{ASCII_TMP_ROOT_WSL}/timing.res",
            start_iteration=PROJECT_STEADY_START_ITERATION,
            end_iteration=PROJECT_STEADY_END_ITERATION,
            restart_read=restart_read_wsl,
            pressure_convergence_threshold_pa=pressure_threshold,
            velocity_convergence_threshold_m_s=velocity_threshold,
        )
        diagnostic_lua_path = run_root / "diagnostic_musubi.lua"
        diagnostic_lua_path.write_text(diagnostic_lua, encoding="utf-8")
        lua_contract = diagnostic_lua_contract(
            production_lua,
            diagnostic_lua,
            start_iteration=PROJECT_STEADY_START_ITERATION,
            end_iteration=PROJECT_STEADY_END_ITERATION,
            restart_read=restart_read_wsl,
            pressure_convergence_threshold_pa=pressure_threshold,
            velocity_convergence_threshold_m_s=velocity_threshold,
        )
        write_json(qc_dir / "diagnostic_musubi_contract.json", lua_contract)
        if lua_contract["status"] != "PASS":
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
                "Frozen physics/BC or project convergence contract failed",
            )
        luac = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=run_root,
            command=[str(environment.binaries["lua_compiler"]), "-p", "diagnostic_musubi.lua"],
            stdout_path=tracking / "luac_musubi_stdout.log",
            stderr_path=tracking / "luac_musubi_stderr.log",
            timeout_s=30,
        )
        if luac.returncode != 0:
            raise FlowError(
                "CFD_FLOW_DIAGNOSTIC_CONFIG_INVALID",
                "Diagnostic Musubi Lua syntax failed",
            )
        diagnostic_config_wsl = windows_to_wsl(
            diagnostic_lua_path,
            config.apes.wsl_distribution,
        )
        mpi_command = [
            "env",
            "OMPI_ALLOW_RUN_AS_ROOT=1",
            "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            str(environment.binaries["mpi_launcher"]),
            "-np",
            "8",
            str(environment.binaries["musubi"]),
            diagnostic_config_wsl,
        ]
        summary.update(
            {
                "status": "PROJECT_STEADY_CONFIRMATION_MUSUBI_RUNNING",
                "short_restart_musubi_run_count": 1,
                "restart_continuation_start_iteration": PROJECT_STEADY_START_ITERATION,
                "requested_end_iteration": PROJECT_STEADY_END_ITERATION,
                "requested_additional_iterations": PROJECT_STEADY_ADDITIONAL_ITERATIONS,
            }
        )
        write_json(summary_path, summary)
        short_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=formal_source / "musubi",
            command=mpi_command,
            stdout_path=tracking / "musubi_stdout.log",
            stderr_path=tracking / "musubi_stderr.log",
            timeout_s=PROJECT_STEADY_WRAPPER_TIMEOUT_S,
        )
        summary.update(
            {
                "status": "PROJECT_STEADY_CONFIRMATION_MUSUBI_FINISHED",
                "short_run_wall_time_s": short_run.wall_time_s,
                "musubi_returncode": short_run.returncode,
            }
        )
        write_json(summary_path, summary)
        if (
            short_run.returncode != 0
            and not _continuation_advanced_from_stdout(
                short_run.stdout_path,
                start_iteration=PROJECT_STEADY_START_ITERATION,
            )
        ):
            return _finalize_project_steady_confirmation(root, run_root)
        return _copy_cleanup_finalize_project_steady(root, run_root)
    except Exception as error:
        latest_summary = _read_json_object(summary_path)
        latest_summary["status"] = (
            error.status if isinstance(error, FlowError) else "CFD_FLOW_DIAGNOSTIC_INTERNAL_ERROR"
        )
        latest_summary["failure"] = str(error)
        write_json(summary_path, latest_summary)
        raise
