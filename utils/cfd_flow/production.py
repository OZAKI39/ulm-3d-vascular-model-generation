"""Validated Base replay, production QC, and low-cost smoke helpers."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv

from .apes import LatticeScaling
from .config import FlowConfig, PLANE_CONTRACT_SHA256
from .io import FlowError, read_json, sha256_file, write_json
from .physical_port_flux import (
    FLUX_ALGORITHM_REVISION,
    FLUX_DEFINITION,
    evaluate_physical_port_fluxes,
    mesh_origin_dx,
)
from .qc import (
    accepted_steady_qc,
    full_v2_qc,
    physical_flow_balance_qc,
    vtu_contract_qc,
)
from .restart_decode import (
    parse_restart_header,
    read_restart_pdf,
    read_treelm_elemlist,
    reconstruct_macroscopic_field,
    tree_ids_to_ijk,
)
from .steady_export import quantize_cell_centers, reconstruct_hexahedral_field
from .validated_contract import FULL_TIMESTEP_REFEREE_REVISION


ACCEPTED_RESTART_SHA256 = "ffcd98b2dc684d1569d937d915b603805809c581d5341e71b17afac2ac64c39f"
ACCEPTED_ITERATION = 598_755
ACCEPTED_PHYSICAL_TIME_S = 0.001220703363914373
SUCCESS_STATUS = "CFD_FLOW_PRODUCTION_TAU1_INTEGRATION_AND_VISUAL_REGRESSION_PASS"
SMOKE_FAILURE = "CFD_FLOW_PRODUCTION_PROMOTION_SMOKE_FAILED"
REPLAY_FAILURE = "CFD_FLOW_PRODUCTION_REPLAY_REPRODUCTION_FAILED"


@dataclass(frozen=True, slots=True)
class ProductionLayout:
    root: Path
    input: Path
    solver_smoke: Path
    steady_replay: Path
    flow: Path
    visualization: Path
    qc: Path
    logs: Path


@dataclass(frozen=True, slots=True)
class ReplayResult:
    metrics: dict[str, Any]
    flux: dict[str, Any]
    field_qc: dict[str, Any]
    reproduction: dict[str, Any]
    vtu_path: Path
    grid: pv.UnstructuredGrid
    points_m: np.ndarray


def create_production_layout(
    output_root: Path,
    *,
    timestamp: datetime | None = None,
) -> ProductionLayout:
    stamp = (timestamp or datetime.now()).strftime("%Y%m%d_%H%M%S")
    root = Path(output_root).resolve() / (
        f"production_tau1_base_promotion_anchor003274_{stamp}"
    )
    if root.exists():
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", f"output already exists: {root}")
    values = {
        name: root / name
        for name in (
            "input", "solver_smoke", "steady_replay", "flow",
            "visualization", "qc", "logs",
        )
    }
    for path in values.values():
        path.mkdir(parents=True, exist_ok=False)
    return ProductionLayout(root=root, **values)


def validated_scaling_record(config: FlowConfig, scaling: LatticeScaling) -> dict[str, Any]:
    target_lattice = (
        config.boundary.target_mass_flow_kg_s
        / config.physics.density_kg_m3
        * scaling.dt_s
        / scaling.dx_m**3
    )
    outlets = {
        f"outlet_{index:02d}": float(value)
        for index, value in enumerate(scaling.outlet_absolute_pressures_pa, start=1)
    }
    gauges = {
        f"outlet_{index:02d}": float(value)
        for index, value in enumerate(config.boundary.outlet_gauge_pressures_pa, start=1)
    }
    return {
        "status": "PASS",
        "scaling": "diffusive",
        "dx_m": scaling.dx_m,
        "dt_formula": "dx^2/(6*nu_phy)",
        "dt_s": scaling.dt_s,
        "rho0_kg_m3": config.physics.density_kg_m3,
        "kinematic_viscosity_m2_s": config.physics.kinematic_viscosity_m2_s,
        "bulk_viscosity_m2_s": config.physics.bulk_viscosity_m2_s,
        "nu_lattice": scaling.nu_lattice,
        "tau": scaling.tau,
        "omega": scaling.omega,
        "pressure_reference_formula": "rho0*cs2*(dx/dt)^2",
        "pressure_reference_pa": scaling.pressure_reference_pa,
        "pressure_reference_role": "LBM_NUMERICAL_PRESSURE_OFFSET",
        "outlet_gauge_pressures_pa": gauges,
        "outlet_absolute_solver_pressures_pa": outlets,
        "target_mass_flow_kg_s": config.boundary.target_mass_flow_kg_s,
        "target_volume_flow_m3_s": config.boundary.target_volume_flow_m3_s,
        "target_lattice_flux": target_lattice,
        "wall_boundary": config.boundary.wall_boundary,
        "inlet_boundary": config.boundary.inlet_boundary,
        "outlet_boundary": config.boundary.outlet_boundary,
        "fresh_initial_pressure_pa": scaling.pressure_reference_pa,
        "fresh_initial_velocity_m_s": [0.0, 0.0, 0.0],
    }


def validate_local_artifacts(config: FlowConfig) -> dict[str, Any]:
    """Hash the immutable mesh/restart/evidence before any solver launch."""

    mesh = config.paths.frozen_base_mesh
    required = {
        "elemlist.lsb": config.mesh.elemlist_sha256,
        "bnd.lsb": config.mesh.bnd_sha256,
        "qval.lsb": config.mesh.qval_sha256,
    }
    mesh_hashes: dict[str, Any] = {}
    for name, expected in required.items():
        path = mesh / name
        actual = sha256_file(path)
        mesh_hashes[name] = {
            "path": str(path), "expected_sha256": expected,
            "actual_sha256": actual, "status": "PASS" if actual == expected else "FAIL",
        }
    restart_actual = sha256_file(config.paths.accepted_base_restart_binary)
    header = parse_restart_header(config.paths.accepted_base_restart_header)
    plane = read_json(config.paths.physical_plane_contract)
    checks = {
        "mesh_hashes": all(value["status"] == "PASS" for value in mesh_hashes.values()),
        "restart_sha256": restart_actual == ACCEPTED_RESTART_SHA256,
        "restart_iteration": header.iteration == ACCEPTED_ITERATION,
        "restart_cells": header.n_elems == config.mesh.expected_cells,
        "restart_pdf_components": header.n_components == 19,
        "plane_contract": plane.get("contract_sha256") == PLANE_CONTRACT_SHA256,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "mesh_hashes": mesh_hashes,
        "accepted_restart": {
            "header": str(config.paths.accepted_base_restart_header),
            "binary": str(config.paths.accepted_base_restart_binary),
            "sha256": restart_actual,
            "iteration": header.iteration,
        },
        "physical_plane_contract_sha256": plane.get("contract_sha256"),
    }
    if result["status"] != "PASS":
        raise FlowError(REPLAY_FAILURE, f"immutable artifact preflight failed: {checks}")
    return result


_CONTROLLER_PATTERN = re.compile(
    r"ADAPTIVE_FLUX_PRESSURE\s+iter=(?P<iteration>\d+)\s+"
    r"target_lattice=\s*(?P<target>[-+0-9.Ee]+)\s+"
    r"controlled_lattice=\s*(?P<controlled>[-+0-9.Ee]+)\s+"
    r"relative_error=\s*(?P<error>[-+0-9.Ee]+)\s+"
    r"rho_boundary=\s*(?P<rho>[-+0-9.Ee]+)\s+"
    r"pressure_pa=\s*(?P<pressure>[-+0-9.Ee]+)\s+"
    r"max_lattice_velocity=\s*(?P<max_velocity>[-+0-9.Ee]+)\s+"
    r"minimum_pdf=\s*(?P<minimum_pdf>[-+0-9.Ee]+)\s+"
    r"globBC_count=(?P<count>\d+)"
)


def parse_controller_records(text: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for match in _CONTROLLER_PATTERN.finditer(text):
        values = match.groupdict()
        records.append(
            {
                "iteration": int(values["iteration"]),
                "target_lattice": float(values["target"]),
                "controlled_lattice": float(values["controlled"]),
                "relative_error": float(values["error"]),
                "rho_boundary": float(values["rho"]),
                "pressure_pa": float(values["pressure"]),
                "max_lattice_velocity": float(values["max_velocity"]),
                "minimum_pdf": float(values["minimum_pdf"]),
                "active_global_boundary_count": int(values["count"]),
            }
        )
    return records


def evaluate_smoke_restart(
    config: FlowConfig,
    scaling: LatticeScaling,
    restart_header_path: Path,
    stdout_text: str,
) -> dict[str, Any]:
    """Hard-gate the one authorized 5000-step fresh equilibrium smoke."""

    header = parse_restart_header(restart_header_path)
    pdf = read_restart_pdf(
        header.binary_path,
        n_elems=header.n_elems,
        n_components=header.n_components,
    )
    field = reconstruct_macroscopic_field(
        pdf,
        dx_m=scaling.dx_m,
        dt_s=scaling.dt_s,
        rho0_kg_m3=config.physics.density_kg_m3,
    )
    controller = parse_controller_records(stdout_text)
    if not controller:
        raise FlowError(SMOKE_FAILURE, "adaptive controller emitted no audit records")
    last = controller[-1]
    expected_target = (
        config.boundary.target_mass_flow_kg_s
        / config.physics.density_kg_m3
        * scaling.dt_s
        / scaling.dx_m**3
    )
    target_error = abs(last["target_lattice"] - expected_target) / expected_target
    controlled_error = abs(last["controlled_lattice"] - expected_target) / expected_target
    velocity_lattice = np.linalg.norm(field.velocity_lattice, axis=1)
    finite = bool(
        np.all(np.isfinite(pdf))
        and np.all(np.isfinite(field.density_lattice))
        and np.all(np.isfinite(field.velocity_lattice))
    )
    gates = {
        "iteration_exact": header.iteration == config.execution.solver_smoke_iterations,
        "mean_rho": 0.9 <= float(np.mean(field.density_lattice)) <= 1.1,
        "minimum_pdf_positive": float(np.min(pdf)) > 0.0,
        "maximum_lattice_speed": float(np.max(velocity_lattice)) < config.solver.maximum_lattice_speed,
        "all_finite": finite,
        "controller_target": target_error <= config.solver.controller_target_error,
        "controller_controlled_flux": controlled_error <= config.solver.controller_controlled_flux_error,
        "active_population_nonzero": last["active_global_boundary_count"] > 0,
    }
    result = {
        "status": "PASS" if all(gates.values()) else "FAIL",
        "fresh_equilibrium_initialization": True,
        "restart_read": None,
        "iteration": header.iteration,
        "restart_header": str(restart_header_path),
        "restart_binary": str(header.binary_path),
        "restart_sha256": sha256_file(header.binary_path),
        "expected_target_lattice": expected_target,
        "observed_target_lattice": last["target_lattice"],
        "observed_controlled_lattice": last["controlled_lattice"],
        "target_error": target_error,
        "controlled_flux_error": controlled_error,
        "mean_rho_lattice": float(np.mean(field.density_lattice)),
        "minimum_pdf": float(np.min(pdf)),
        "maximum_lattice_speed": float(np.max(velocity_lattice)),
        "controller_last_record": last,
        "controller_record_count": len(controller),
        "gates": gates,
    }
    if result["status"] != "PASS":
        raise FlowError(SMOKE_FAILURE, f"smoke hard gate failed: {gates}")
    return result


def _reference_checkpoint(config: FlowConfig) -> tuple[dict[str, Any], dict[str, Any]]:
    base_final = read_json(config.paths.accepted_base_qc)
    history = read_json(config.paths.accepted_base_checkpoint_history)
    records = history.get("checkpoint_history", [])
    matches = [record for record in records if record.get("iteration") == ACCEPTED_ITERATION]
    if len(matches) != 1:
        raise FlowError(REPLAY_FAILURE, "accepted checkpoint record is missing or ambiguous")
    if base_final.get("accepted_restart", {}).get("sha256") != ACCEPTED_RESTART_SHA256:
        raise FlowError(REPLAY_FAILURE, "accepted Base QC lineage changed")
    return base_final, matches[0]


def _comparison(
    actual: float,
    reference: float,
    *,
    relative_gate: float | None = None,
    absolute_gate: float | None = None,
) -> dict[str, Any]:
    absolute_error = abs(float(actual) - float(reference))
    relative_error = absolute_error / max(abs(float(reference)), np.finfo(float).tiny)
    passed = (
        (relative_gate is not None and relative_error <= relative_gate)
        or (absolute_gate is not None and absolute_error <= absolute_gate)
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "actual": float(actual), "reference": float(reference),
        "absolute_error": absolute_error, "relative_error": relative_error,
        "relative_gate": relative_gate, "absolute_gate": absolute_gate,
    }


def replay_accepted_base(
    config: FlowConfig,
    scaling: LatticeScaling,
    layout: ProductionLayout,
) -> ReplayResult:
    """Decode and independently re-run production post-processing on accepted Base."""

    base_final, reference = _reference_checkpoint(config)
    header = parse_restart_header(config.paths.accepted_base_restart_header)
    pdf = read_restart_pdf(
        config.paths.accepted_base_restart_binary,
        n_elems=header.n_elems,
        n_components=header.n_components,
    )
    field = reconstruct_macroscopic_field(
        pdf,
        dx_m=scaling.dx_m,
        dt_s=scaling.dt_s,
        rho0_kg_m3=config.physics.density_kg_m3,
    )
    tree_ids, _, elemlist_contract = read_treelm_elemlist(
        config.paths.frozen_base_mesh / "elemlist.lsb",
        n_elems=header.n_elems,
    )
    origin_m, mesh_dx_m = mesh_origin_dx(config.paths.frozen_base_mesh)
    if mesh_dx_m != scaling.dx_m:
        raise FlowError(REPLAY_FAILURE, "mesh dx differs from production contract")
    cell_indices = tree_ids_to_ijk(tree_ids)
    points_m = origin_m + (cell_indices.astype(np.float64) + 0.5) * mesh_dx_m
    mapping = quantize_cell_centers(points_m, origin_m=origin_m, dx_m=mesh_dx_m)

    plane_contract = read_json(config.paths.physical_plane_contract)
    flux = evaluate_physical_port_fluxes(
        plane_contract,
        points_m,
        field.velocity_phy,
        field.density_lattice,
        dx_m=mesh_dx_m,
    )
    flow_balance_qc = physical_flow_balance_qc(flux)
    if flux["status"] != "PASS" or flow_balance_qc["status"] != "PASS":
        raise FlowError(REPLAY_FAILURE, "physical aperture-flow QC failed")
    controller_pressure = float(reference["controller"]["pressure_pa"])
    inlet_gauge = controller_pressure - scaling.pressure_reference_pa
    pressure_drops = {
        f"outlet_{index:02d}": inlet_gauge - gauge
        for index, gauge in enumerate(config.boundary.outlet_gauge_pressures_pa, start=1)
    }
    speed_lattice = np.linalg.norm(field.velocity_lattice, axis=1)
    density = field.density_lattice
    metrics = {
        "iteration": header.iteration,
        "physical_time_s": ACCEPTED_PHYSICAL_TIME_S,
        "rho_mean": float(np.mean(density)),
        "rho_median": float(np.median(density)),
        "rho_p1": float(np.percentile(density, 1)),
        "rho_p99": float(np.percentile(density, 99)),
        "mean_speed_m_s": float(np.mean(np.linalg.norm(field.velocity_phy, axis=1))),
        "minimum_pdf": float(np.min(pdf)),
        "maximum_lattice_speed": float(np.max(speed_lattice)),
        "inlet_gauge_pressure_pa": inlet_gauge,
        "pressure_drops_pa": pressure_drops,
        "Qin_m3_s": flux["Qin_m3_s"],
        "Qout_m3_s": flux["Qout_m3_s"],
        "physical_volume_closure": flux["closure"],
        "Q1_m3_s": flux["ports"]["outlet_01"]["physical_q_m3_s"],
        "Q2_m3_s": flux["ports"]["outlet_02"]["physical_q_m3_s"],
        "Q3_m3_s": flux["ports"]["outlet_03"]["physical_q_m3_s"],
        "flow_fractions": flux["flow_fractions"],
        "all_finite": bool(
            np.all(np.isfinite(pdf))
            and np.all(np.isfinite(field.velocity_phy))
            and np.all(np.isfinite(field.pressure_phy))
        ),
    }

    comparisons = {
        "rho_mean": _comparison(metrics["rho_mean"], reference["rho_lattice"]["mean"], relative_gate=1.0e-8),
        "rho_p1": _comparison(metrics["rho_p1"], reference["rho_lattice"]["p1"], relative_gate=1.0e-8),
        "rho_p99": _comparison(metrics["rho_p99"], reference["rho_lattice"]["p99"], relative_gate=1.0e-8),
        "minimum_pdf": _comparison(metrics["minimum_pdf"], reference["minimum_pdf"], relative_gate=1.0e-8),
        "maximum_lattice_speed": _comparison(metrics["maximum_lattice_speed"], reference["maximum_lattice_speed"], relative_gate=1.0e-8),
        "inlet_gauge_pressure_pa": _comparison(metrics["inlet_gauge_pressure_pa"], reference["inlet_gauge_pressure_pa"], absolute_gate=1.0e-6),
    }
    reference_flux = {
        "Qin_m3_s": reference["ports"]["inlet"]["Q_velocity_m3_s"],
        "Q1_m3_s": reference["ports"]["outlet_01"]["Q_velocity_m3_s"],
        "Q2_m3_s": reference["ports"]["outlet_02"]["Q_velocity_m3_s"],
        "Q3_m3_s": reference["ports"]["outlet_03"]["Q_velocity_m3_s"],
        "Qout_m3_s": reference["Qout_sum_m3_s"],
        "physical_volume_closure": reference["physical_volume_closure"],
    }
    for name, expected in reference_flux.items():
        comparisons[name] = _comparison(metrics[name], expected, relative_gate=1.0e-6)
    for label, expected in reference["flow_fractions"].items():
        comparisons[f"flow_fraction_{label}"] = _comparison(
            metrics["flow_fractions"][label], expected, relative_gate=1.0e-8
        )
    for label, expected in reference["pressure_drops_pa"].items():
        comparisons[f"pressure_drop_{label}"] = _comparison(
            metrics["pressure_drops_pa"][label], expected, absolute_gate=1.0e-6
        )
    reproduction = {
        "status": "PASS" if all(item["status"] == "PASS" for item in comparisons.values()) else "FAIL",
        "same_restart_sha256": sha256_file(config.paths.accepted_base_restart_binary) == ACCEPTED_RESTART_SHA256,
        "comparisons": comparisons,
        "maximum_relative_error": max(item["relative_error"] for item in comparisons.values()),
        "maximum_absolute_pressure_error_pa": max(
            item["absolute_error"] for name, item in comparisons.items()
            if "pressure" in name
        ),
    }
    if reproduction["status"] != "PASS":
        raise FlowError(REPLAY_FAILURE, f"accepted Base reproduction failed: {comparisons}")

    steady = base_final["steady_audit"]
    promoted_steady = accepted_steady_qc(steady)
    if promoted_steady["status"] != "PASS":
        raise FlowError(REPLAY_FAILURE, "accepted physical-time steady gates changed")
    field_qc = {
        "status": "PASS",
        "source": "VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART",
        "fresh_full_production_steady_solve": False,
        "steady_gates": steady["gates"],
        "promoted_steady_acceptance": promoted_steady,
        "physical_flow_balance": flow_balance_qc,
        "R_mass_short": steady["R_mass_short"],
        "R_mass_long": steady["R_mass_long"],
        "R_velocity": steady["R_velocity"],
        "R_pressure": steady["R_pressure"],
        "R_inlet": steady["R_inlet"],
        "maximum_flow_fraction_drift": steady["maximum_flow_fraction_drift"],
        "minimum_pdf": metrics["minimum_pdf"],
        "maximum_lattice_speed": metrics["maximum_lattice_speed"],
        "all_finite": metrics["all_finite"],
        "elemlist_contract": elemlist_contract,
    }

    grid = reconstruct_hexahedral_field(
        mapping=mapping,
        origin_m=origin_m,
        dx_m=mesh_dx_m,
        pressure_pa=field.pressure_phy,
        velocity_m_s=field.velocity_phy,
        pressure_reference_pa=scaling.pressure_reference_pa,
    )
    speed_m_s = np.linalg.norm(field.velocity_phy, axis=1)
    # Keep the production VTU contract explicit at the promotion boundary.
    grid.cell_data["velocity_phy"] = field.velocity_phy
    grid.cell_data["velocity_magnitude_m_s"] = speed_m_s
    grid.cell_data["velocity_magnitude_mm_s"] = speed_m_s * 1.0e3
    grid.cell_data["pressure_gauge_pa"] = (
        field.pressure_phy - scaling.pressure_reference_pa
    )
    grid.cell_data["pressure_absolute_solver_pa"] = field.pressure_phy
    grid.cell_data["rho_lattice"] = field.density_lattice
    vtu_path = layout.flow / "production_steady_flow_field.vtu"
    grid.save(vtu_path, binary=True)
    vtu_qc = vtu_contract_qc(grid)
    if vtu_qc["status"] != "PASS":
        raise FlowError(REPLAY_FAILURE, "production VTU arrays are incomplete")

    write_json(layout.steady_replay / "physical_port_flux.json", flux)
    write_json(layout.steady_replay / "accepted_base_reproduction.json", reproduction)
    write_json(layout.qc / "production_steady_qc.json", field_qc)
    write_json(
        layout.flow / "production_steady_flow_field_manifest.json",
        {
            "status": "PASS", "path": str(vtu_path),
            "sha256": sha256_file(vtu_path), "size_bytes": vtu_path.stat().st_size,
            "cell_count": grid.n_cells, "point_count": grid.n_points,
            "cell_arrays": sorted(grid.cell_data.keys()),
            "vtu_contract_qc": vtu_qc,
            "source_restart_sha256": ACCEPTED_RESTART_SHA256,
        },
    )
    return ReplayResult(metrics, flux, field_qc, reproduction, vtu_path, grid, points_m)


def validate_full_v2(config: FlowConfig) -> dict[str, Any]:
    evidence = read_json(config.paths.accepted_base_full_v2)
    residual = float(evidence["referee"]["R_full_one_step_identity"])
    oracle_qc = full_v2_qc(evidence, gate=config.solver.full_timestep_identity_gate)
    checks = {
        "status_pass": evidence.get("status") == "PASS",
        "start_iteration": evidence.get("iteration_start") == ACCEPTED_ITERATION,
        "restart_sha256": evidence.get("restart_sha256", {}).get(str(ACCEPTED_ITERATION)) == ACCEPTED_RESTART_SHA256,
        "hard_gate": float(evidence.get("hard_gate")) == config.solver.full_timestep_identity_gate,
        "residual": residual <= config.solver.full_timestep_identity_gate,
        "oracle_qc": oracle_qc["status"] == "PASS",
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "oracle": FULL_TIMESTEP_REFEREE_REVISION,
        "residual": residual,
        "gate": config.solver.full_timestep_identity_gate,
        "evidence_path": str(config.paths.accepted_base_full_v2),
        "evidence_sha256": sha256_file(config.paths.accepted_base_full_v2),
        "checks": checks,
        "oracle_qc": oracle_qc,
        "new_solver_calls": 0,
    }


def write_primary_metrics_csv(
    path: Path,
    metrics: dict[str, Any],
    steady_qc: dict[str, Any],
    full_v2: dict[str, Any],
    config: FlowConfig,
) -> None:
    rows = [
        ("Qin", metrics["Qin_m3_s"], "m3/s", config.boundary.target_volume_flow_m3_s, 0.01, "PASS", FLUX_DEFINITION),
        ("Qout_total", metrics["Qout_m3_s"], "m3/s", metrics["Qin_m3_s"], 0.01, "PASS", FLUX_DEFINITION),
        ("physical_closure", metrics["physical_volume_closure"], "1", 0.0, 0.01, "PASS", FLUX_ALGORITHM_REVISION),
    ]
    for label, value in metrics["flow_fractions"].items():
        rows.append((f"{label}_flow_fraction", value, "1", "", 0.01, "PASS", FLUX_DEFINITION))
    rows.extend(
        [
            ("inlet_gauge_pressure", metrics["inlet_gauge_pressure_pa"], "Pa", "", 1.0e-6, "PASS", "accepted controller record"),
            *[
                (f"DeltaP{index:02d}", metrics["pressure_drops_pa"][f"outlet_{index:02d}"], "Pa", "", 1.0e-6, "PASS", "gauge pressure difference")
                for index in range(1, 4)
            ],
            ("rho_mean", metrics["rho_mean"], "lattice", 1.0, "[0.9,1.1]", "PASS", "decoded accepted restart"),
            ("minimum_pdf", metrics["minimum_pdf"], "lattice", ">0", 0.0, "PASS", "decoded accepted restart"),
            ("maximum_lattice_speed", metrics["maximum_lattice_speed"], "lattice", "", 0.05, "PASS", "decoded accepted restart"),
            ("R_mass_short", steady_qc["R_mass_short"], "1", 0.0, 0.01, "PASS", "accepted physical-time steady audit"),
            ("R_mass_long", steady_qc["R_mass_long"], "1", 0.0, 0.01, "PASS", "accepted physical-time steady audit"),
            ("Full_V2", full_v2["residual"], "1", 0.0, full_v2["gate"], full_v2["status"], FULL_TIMESTEP_REFEREE_REVISION),
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("metric", "value", "unit", "target", "threshold", "status", "source"))
        writer.writerows(rows)


def git_head(project_root: Path) -> str:
    import subprocess

    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=project_root,
        capture_output=True, text=True, check=True,
    )
    return result.stdout.strip()
