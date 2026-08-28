"""Low-cost read-only diagnostics for the frozen Musubi restart."""

from __future__ import annotations

import csv
import math
import re
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
RESTART_INVALID = "CFD_FLOW_RESTART_INVALID"
HARVEST_FAILED = "CFD_FLOW_DIAGNOSTIC_HARVEST_FAILED"
SOURCE_RUN_NAME = "musubi_solver_recovery_anchor003274_20260828_164912"
FROZEN_SEEDER_NAME = "musubi_recovery_anchor003274_20260828_162530"
START_ITERATION = 162_464
ADDITIONAL_ITERATIONS = 10_000
END_ITERATION = START_ITERATION + ADDITIONAL_ITERATIONS
DIAGNOSTIC_WALLCLOCK_S = 600
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
    tracking_folder_wsl: str,
    restart_folder_wsl: str,
    timing_file_wsl: str,
) -> str:
    """Change only restart, stop ceiling, tracking, and output destinations."""

    text = production_lua
    text = _replace_once(
        text,
        "-- Generated production configuration; official Musubi syntax.",
        "-- Diagnostic restart continuation; frozen production physics and BC.",
    )
    text = _replace_once(text, "timing_file = 'mus_timing.res'", f"timing_file = '{timing_file_wsl}'")
    text = _replace_once(text, "maximum_iterations = 1000000", f"maximum_iterations = {END_ITERATION}")
    text = _replace_once(
        text,
        "max = { iter = maximum_iterations, clock = 3600 }",
        f"max = {{ iter = maximum_iterations, clock = {DIAGNOSTIC_WALLCLOCK_S} }}",
    )
    tracking_folder = tracking_folder_wsl.rstrip("/") + "/"
    restart_folder = restart_folder_wsl.rstrip("/") + "/"
    diagnostic_output = f"""tracking = {{
  {{
    label = 'domain_average_pressure',
    folder = '{tracking_folder}',
    variable = {{ 'pressure_phy' }},
    shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {START_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = 100 }} }},
    output = {{ format = 'ascii' }}
  }},
  {{
    label = 'domain_average_velocity_magnitude',
    folder = '{tracking_folder}',
    variable = {{ 'vel_mag_phy' }},
    shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {START_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = 100 }} }},
    output = {{ format = 'ascii' }}
  }}
}}

restart = {{
  read = 'restart/roi003274_steady_lbm_lastHeader.lua',
  write = '{restart_folder}'
}}
"""
    text = _replace_once(text, "restart = { write = 'restart/' }\n", diagnostic_output)
    return text


def diagnostic_lua_contract(production_lua: str, diagnostic_lua: str) -> dict[str, Any]:
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
        "restart_read_present": "read = 'restart/roi003274_steady_lbm_lastHeader.lua'" in diagnostic_lua,
        "iteration_ceiling_exact": f"maximum_iterations = {END_ITERATION}" in diagnostic_lua,
        "wallclock_ceiling_exact": f"clock = {DIAGNOSTIC_WALLCLOCK_S}" in diagnostic_lua,
        "pressure_tracking_ascii": (
            "label = 'domain_average_pressure'" in diagnostic_lua
            and "variable = { 'pressure_phy' }" in diagnostic_lua
        ),
        "velocity_tracking_ascii": (
            "label = 'domain_average_velocity_magnitude'" in diagnostic_lua
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


def _load_tracking_series(directory: Path, label: str, dt_s: float) -> dict[int, float]:
    files = sorted(path for path in directory.rglob("*.res") if label in path.name)
    if not files:
        raise FlowError("CFD_FLOW_DIAGNOSTIC_TRACKING_INVALID", f"No tracking file for {label}")
    values: dict[int, float] = {}
    for path in files:
        for time_value, result in _numeric_tracking_rows(path):
            iteration = int(round(time_value / dt_s)) if abs(time_value) < 1_000.0 else int(round(time_value))
            values[iteration] = result
    if len(values) < 2:
        raise FlowError("CFD_FLOW_DIAGNOSTIC_TRACKING_INVALID", f"Too few samples for {label}")
    return values


def parse_tracking(directory: Path, dt_s: float) -> list[dict[str, float | int | None]]:
    pressure = _load_tracking_series(directory, "domain_average_pressure", dt_s)
    velocity = _load_tracking_series(directory, "domain_average_velocity_magnitude", dt_s)
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


def _window_metrics(values: np.ndarray, threshold: float) -> dict[str, Any]:
    absolute = np.abs(values)
    if len(absolute) < 50:
        raise FlowError("CFD_FLOW_DIAGNOSTIC_TRACKING_INVALID", "Fewer than 50 tracking deltas")
    windows = {
        str(width): float(np.mean(absolute[-width:]))
        for width in (1, 5, 10, 20, 50)
    }
    early = float(np.mean(absolute[:20]))
    late = windows["20"]
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
        "last_sample_delta": float(values[-1]),
        "mean_absolute_delta_by_last_n_samples": windows,
        "last_20_mean_absolute_delta": late,
        "last_20_max_absolute_delta": float(np.max(absolute[-20:])),
        "threshold": threshold,
        "ratio_to_threshold": late / threshold,
        "early_20_mean_absolute_delta": early,
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
            tracking_folder_wsl=tracking_wsl,
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
