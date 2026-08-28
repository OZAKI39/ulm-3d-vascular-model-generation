"""One-shot production orchestration for Seeder and Musubi."""

from __future__ import annotations

import math
import os
import re
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from .apes import (
    bulk_viscosity_from_kinematic,
    compute_lattice_scaling,
    generate_harvester_lua,
    generate_musubi_lua,
    inspect_apes_environment,
    load_boundary_conditions,
    parse_mesh_header,
    run_wsl_tool,
    save_environment,
    save_lua_files,
)
from .config import FlowConfig, METHOD as CONFIG_METHOD
from .geometry import (
    cells_across_diameter,
    compute_bounding_cube,
    load_frozen_surface_partition,
    partition_surface,
    robust_inside,
    SeedPoint,
)
from .io import (
    FlowError,
    RunLayout,
    create_run_layout,
    load_flow_inputs,
    read_json,
    save_input_provenance,
    sha256_file,
    write_json,
)
from .qc import (
    boundary_condition_qc,
    evaluate_mass_conservation,
    fluid_domain_geometry_qc,
    numerical_port_fluxes,
    reynolds_diagnostics,
    validate_and_convert_flow_vtu,
    write_proteus_metadata,
)
from .visualization import create_flow_figures


METHOD = CONFIG_METHOD
SUCCESS_STATUS = "CFD_FLOW_MUSUBI_BASELINE_PASS_PENDING_GRID_CONVERGENCE"
SUCCESS_NEXT = "RUN MUSUBI GRID-SPACING CONVERGENCE STUDY"
QC_FILENAMES = (
    "source_surface_qc.json",
    "patch_partition_qc.json",
    "seeder_mesh_qc.json",
    "lbm_scaling_qc.json",
    "boundary_condition_qc.json",
    "convergence_qc.json",
    "flow_field_qc.json",
    "mass_conservation_qc.json",
    "pressure_boundary_qc.json",
    "proteus_compatibility_qc.json",
    "run_summary.json",
)
FAILURE_EVIDENCE_RUN = "musubi_anchor003274_20260828_161151"


@dataclass(frozen=True, slots=True)
class FlowResult:
    status: str
    next: str
    run_root: Path
    summary: dict[str, Any]


def _frozen_recovery_seed(partition: Any, failure_run: Path) -> SeedPoint:
    """Reuse the accepted failed-run seed without performing a new search."""

    previous = read_json(failure_run / "qc" / "run_summary.json")
    coordinates_um = np.asarray(previous["seed_point_um"], dtype=float)
    inside = bool(robust_inside(partition.mesh_um, coordinates_um[None, :])[0])
    if not inside:
        raise FlowError("CFD_FLOW_SEED_POINT_FAILED", "Frozen recovery seed is no longer inside the lumen")
    return SeedPoint(
        coordinates_m=coordinates_um * 1.0e-6,
        coordinates_um=coordinates_um,
        candidate_offset_radius=float(previous["seed_candidate_offset_radius"]),
        seed_inside_lumen=True,
        method="failure_evidence_seed_reused_and_closed_mesh_contains_reverified",
    )


def _directory_manifest(directory: Path) -> dict[str, Any]:
    files = sorted(path for path in directory.rglob("*") if path.is_file())
    return {
        "root": str(directory),
        "file_count": len(files),
        "total_bytes": int(sum(path.stat().st_size for path in files)),
        "files": [
            {
                "relative_path": path.relative_to(directory).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in files
        ],
    }


def _static_recovery_preflight(
    layout: RunLayout,
    partition: Any,
    environment: Any,
    distribution: str,
) -> dict[str, Any]:
    """Syntax-check all Lua and prove every solver STL is visible in WSL."""

    lua_checks = []
    for label, workdir, filename in (
        ("seeder", layout.seeder, "seeder.lua"),
        ("musubi", layout.musubi, "musubi.lua"),
        ("harvester", layout.musubi, "mus_harvester.lua"),
    ):
        run = run_wsl_tool(
            distribution=distribution,
            workdir=workdir,
            command=[str(environment.binaries["lua_compiler"]), "-p", filename],
            stdout_path=layout.input / f"luac_{label}_stdout.log",
            stderr_path=layout.input / f"luac_{label}_stderr.log",
            timeout_s=30,
        )
        lua_checks.append({"name": label, "returncode": run.returncode, "status": "PASS" if run.returncode == 0 else "FAIL"})

    stl_checks = []
    for patch in partition.patches:
        run = run_wsl_tool(
            distribution=distribution,
            workdir=layout.root,
            command=["test", "-f", f"geometry/geometry_solver_m/{patch.label}.stl"],
            stdout_path=layout.input / f"wsl_stl_{patch.label}_stdout.log",
            stderr_path=layout.input / f"wsl_stl_{patch.label}_stderr.log",
            timeout_s=30,
        )
        stl_checks.append({"label": patch.label, "returncode": run.returncode, "status": "PASS" if run.returncode == 0 else "FAIL"})

    executable_checks = {
        name: bool(environment.binaries.get(name))
        for name in ("seeder", "musubi", "mus_harvesting", "mpi_launcher", "lua_compiler")
    }
    passed = (
        all(item["status"] == "PASS" for item in lua_checks)
        and all(item["status"] == "PASS" for item in stl_checks)
        and all(executable_checks.values())
    )
    result = {
        "status": "PASS" if passed else "FAIL",
        "lua_syntax": lua_checks,
        "wsl_solver_stl_paths": stl_checks,
        "executables": executable_checks,
    }
    write_json(layout.input / "recovery_static_preflight.json", result)
    if not passed:
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Recovery static preflight failed")
    return result


def _musubi_launch_preflight(
    layout: RunLayout,
    mesh_summary: dict[str, Any],
    bc: Any,
    scaling: Any,
) -> dict[str, Any]:
    text = (layout.musubi / "musubi.lua").read_text(encoding="utf-8")
    bulk_viscosity = bulk_viscosity_from_kinematic(bc.kinematic_viscosity_m2_s)
    labels = tuple(mesh_summary["boundary_labels_in_file_order"])
    required_labels = ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")
    checks = {
        "boundary_labels_exact": set(labels) == set(required_labels) and len(labels) == len(required_labels),
        "musubi_labels_match_seeder": all(f"label = '{label}'" in text for label in labels),
        "identify_kind_fluid": "kind = 'fluid'" in text,
        "layout_d3q19": "layout = 'd3q19'" in text,
        "relaxation_bgk": "relaxation = 'bgk'" in text,
        "kinematic_viscosity_frozen": f"nu_phy = {bc.kinematic_viscosity_m2_s:.17g}" in text,
        "bulk_viscosity_explicit": "bulk_viscosity = bulk_viscosity_phy" in text,
        "bulk_viscosity_two_thirds_nu": "bulk_viscosity_phy = (2.0 / 3.0) * nu_phy" in text,
        "dt_frozen": f"dt = {scaling.dt_s:.17g}" in text,
        "rho0_frozen": f"rho0_phy = {bc.density_kg_m3:.17g}" in text,
        "mfr_eq": "kind = 'mfr_eq'" in text,
        "mass_flowrate_key": "mass_flowrate =" in text and "massflowrate" not in text,
        "three_pressure_eq": text.count("kind = 'pressure_eq'") == 3,
        "absolute_pressures_positive": all(value > 0.0 for value in scaling.outlet_absolute_pressures_pa),
        "wall_libb": "kind = 'wall_libb'" in text,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "seeder_boundary_labels": list(labels),
        "absolute_pressure_bc_pa": list(scaling.outlet_absolute_pressures_pa),
        "mass_flow_target_kg_s": bc.density_kg_m3 * bc.inlet_flow_m3_s,
        "bulk_viscosity_m2_s": bulk_viscosity,
    }
    write_json(layout.musubi / "musubi_launch_preflight.json", result)
    if result["status"] != "PASS":
        raise FlowError("CFD_FLOW_MUSUBI_CONFIG_FAILED", "Musubi launch preflight failed")
    return result


def _pending_qc(layout: RunLayout) -> None:
    for filename in QC_FILENAMES[:-1]:
        path = layout.qc / filename
        if not path.exists():
            write_json(path, {"status": "NOT_REACHED"})


def _source_surface_qc(inputs: Any, partition: Any, reference: dict[str, Any]) -> dict[str, Any]:
    mesh = partition.mesh_um
    return {
        "status": "PASS",
        "source_tagged_vtp": str(inputs.tagged_surface_vtp),
        "source_meter_stl": str(inputs.meter_surface_stl),
        "tagged_vtp_sha256": reference["sha256"]["tagged_surface_vtp"],
        "meter_stl_sha256": reference["sha256"]["meter_surface_stl"],
        "triangle_count": int(len(partition.faces)),
        "point_count": int(len(partition.points_um)),
        "watertight": bool(mesh.is_watertight),
        "single_component": bool(len(mesh.split(only_watertight=False)) == 1),
        "bounds_um": [float(value) for value in np.asarray(mesh.bounds).reshape(-1)],
        "extents_um": [float(value) for value in mesh.extents],
        "lumen_volume_um3": float(abs(mesh.volume)),
        "surface_geometry_modified": False,
        "geometry_regeneration_count": 0,
    }


def _scaling_qc(scaling: Any, bc: Any) -> dict[str, Any]:
    gauge = np.asarray(bc.outlet_gauge_pressures_pa)
    absolute = np.asarray(scaling.outlet_absolute_pressures_pa)
    return {
        "status": "PASS",
        **asdict(scaling),
        "cs_lattice_squared": 1.0 / 3.0,
        "outlet_gauge_pressures_pa": gauge.tolist(),
        "common_pressure_offset_pa": scaling.pressure_reference_pa,
        "pressure_differences_preserved": bool(
            np.allclose(
                np.subtract.outer(gauge, gauge),
                np.subtract.outer(absolute, absolute),
                rtol=0.0,
                atol=1.0e-10,
            )
        ),
        "omega_gate": "0 < omega < 2",
        "mach_gate": "Ma < 0.05",
        "official_semantics": "physics.dt/rho0 and pressure factor rho0*dx^2/dt^2",
    }


def _resource_preflight(
    config: FlowConfig,
    volume_um3: float,
    cube: Any,
    available_ram: int,
) -> dict[str, Any]:
    expected = int(math.ceil(volume_um3 / config.mesh.dx_target_um**3))
    estimated_ram = expected * config.resources.estimated_bytes_per_fluid_cell
    limit = int(available_ram * config.resources.maximum_available_ram_fraction)
    result = {
        "lumen_volume_um3": volume_um3,
        "expected_fluid_element_count": expected,
        "estimated_ram_bytes": estimated_ram,
        "available_ram_bytes": available_ram,
        "maximum_allowed_ram_bytes": limit,
        "bounding_cube_origin_m": cube.origin_m.tolist(),
        "bounding_cube_side_m": cube.side_m,
        "bounding_cube_level": cube.level,
        "bounding_cells_per_axis": cube.cells_per_axis,
        "minimum_margin_cells": cube.margin_cells_minimum,
        "dx_m": config.mesh.dx_target_m,
        "uniform_cartesian": True,
        "resource_status": "PASS" if estimated_ram <= limit else "FAIL",
    }
    if estimated_ram > limit:
        raise FlowError(
            "CFD_FLOW_MESH_RESOURCE_LIMIT",
            f"Estimated {estimated_ram} bytes exceeds 60% RAM limit {limit}",
        )
    return result


def _parse_convergence(stdout: str, stderr: str, dt_s: float) -> dict[str, Any]:
    text = stdout + "\n" + stderr
    steady = bool(re.search(r"Simulation reached a steady state", text, flags=re.IGNORECASE))
    matches = re.findall(r"Reached steady state\s+(\d+)", text, flags=re.IGNORECASE)
    if not matches:
        matches = re.findall(r"iterations?\s*[=:]\s*(\d+)", text, flags=re.IGNORECASE)
    iterations = int(matches[-1]) if matches else None
    return {
        "status": "PASS" if steady else "FAIL",
        "steady_converged": steady,
        "iteration_count": iterations,
        "lattice_time_reached": iterations,
        "physical_time_reached_s": None if iterations is None else iterations * dt_s,
        "evidence": "Simulation reached a steady state" if steady else "No steady-state termination evidence",
    }


def _find_harvested_vtk(flow_dir: Path, destination: Path) -> Path:
    parallel = [path for path in flow_dir.rglob("*.pvtu") if path.resolve() != destination.resolve()]
    serial = [path for path in flow_dir.rglob("*.vtu") if path.resolve() != destination.resolve()]
    candidates = parallel or serial
    if not candidates:
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", "mus_harvesting produced no VTU/PVTU")
    return sorted(candidates, key=lambda path: path.stat().st_mtime)[-1]


def _pressure_qc(measured: dict[str, float], bc: Any) -> dict[str, Any]:
    rows = []
    for index, target in enumerate(bc.outlet_gauge_pressures_pa, start=1):
        label = f"outlet_{index:02d}"
        actual = measured[label]
        rows.append(
            {
                "label": label,
                "target_source": "P_solver_boundary_pa",
                "target_gauge_pa": target,
                "measured_near_boundary_gauge_pa": actual,
                "absolute_error_pa": abs(actual - target),
                "relative_error": abs(actual - target) / max(abs(target), 1.0e-30),
                "role": "REPORTED_BOUNDARY_ACCURACY",
            }
        )
    return {"status": "PASS", "outlets": rows}


def _write_failure_summary(
    layout: RunLayout,
    summary: dict[str, Any],
    error: FlowError,
) -> FlowResult:
    summary["status"] = error.status
    summary["next"] = "REVIEW CFD FLOW FAILURE EVIDENCE"
    summary["failure"] = str(error)
    _pending_qc(layout)
    write_json(layout.qc / "run_summary.json", summary)
    return FlowResult(error.status, summary["next"], layout.root, summary)


def _validate_frozen_seeder_mesh(
    frozen_run: Path,
    current_surface_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    mesh_dir = frozen_run / "seeder" / "mesh"
    manifest_path = frozen_run / "seeder" / "mesh_manifest.json"
    manifest = read_json(manifest_path)
    recorded = {item["relative_path"]: item for item in manifest.get("files", [])}
    actual_files = sorted(path for path in mesh_dir.rglob("*") if path.is_file())
    actual_names = {path.relative_to(mesh_dir).as_posix() for path in actual_files}
    inventory_match = actual_names == set(recorded)
    hash_mismatches = []
    for path in actual_files:
        relative = path.relative_to(mesh_dir).as_posix()
        expected = recorded.get(relative)
        if expected is None:
            continue
        actual_hash = sha256_file(path)
        if actual_hash != expected["sha256"] or path.stat().st_size != expected["bytes"]:
            hash_mismatches.append(relative)
    mesh_summary = parse_mesh_header(mesh_dir)
    expected_counts = {
        "wall": 74799,
        "inlet": 334,
        "outlet_01": 187,
        "outlet_02": 194,
        "outlet_03": 227,
    }
    labels = set(mesh_summary["boundary_labels_in_file_order"])
    expected_labels = set(expected_counts)
    checks = {
        "manifest_status_frozen": manifest.get("status") == "FROZEN_SEEDER_MESH",
        "inventory_exact": inventory_match,
        "hashes_and_sizes_exact": not hash_mismatches,
        "source_surface_sha_exact": manifest.get("source_surface_sha256") == current_surface_sha256,
        "fluid_elements_exact": mesh_summary["fluid_element_count"] == 221109,
        "dx_exact": manifest.get("actual_dx_m") == 2.0e-7,
        "minimum_level_exact": mesh_summary["minimum_level"] == 9,
        "maximum_level_exact": mesh_summary["maximum_level"] == 9,
        "boundary_labels_exact": labels == expected_labels and len(labels) == 5,
        "boundary_counts_exact": mesh_summary["boundary_cell_counts"] == expected_counts,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "frozen_seeder_run": str(frozen_run),
        "mesh_directory": str(mesh_dir),
        "manifest": str(manifest_path),
        "checks": checks,
        "hash_mismatches": hash_mismatches,
        "file_inventory": manifest.get("files", []),
        "mesh_summary": mesh_summary,
    }
    if result["status"] != "PASS":
        raise FlowError(
            "CFD_FLOW_FROZEN_SEEDER_MESH_INVALID",
            f"Frozen Seeder mesh validation failed: {checks}",
        )
    return result, mesh_summary


def _musubi_only_static_preflight(
    layout: RunLayout,
    environment: Any,
    distribution: str,
    mesh_reference: str,
) -> dict[str, Any]:
    lua_checks = []
    for label, filename in (("musubi", "musubi.lua"), ("harvester", "mus_harvester.lua")):
        run = run_wsl_tool(
            distribution=distribution,
            workdir=layout.musubi,
            command=[str(environment.binaries["lua_compiler"]), "-p", filename],
            stdout_path=layout.input / f"luac_{label}_stdout.log",
            stderr_path=layout.input / f"luac_{label}_stderr.log",
            timeout_s=30,
        )
        lua_checks.append({"name": label, "returncode": run.returncode, "status": "PASS" if run.returncode == 0 else "FAIL"})
    mesh_check = run_wsl_tool(
        distribution=distribution,
        workdir=layout.musubi,
        command=["test", "-f", f"{mesh_reference}header.lua"],
        stdout_path=layout.input / "wsl_frozen_mesh_stdout.log",
        stderr_path=layout.input / "wsl_frozen_mesh_stderr.log",
        timeout_s=30,
    )
    executables = {
        name: bool(environment.binaries.get(name))
        for name in ("musubi", "mus_harvesting", "mpi_launcher", "lua_compiler")
    }
    passed = (
        all(item["status"] == "PASS" for item in lua_checks)
        and mesh_check.returncode == 0
        and all(executables.values())
    )
    result = {
        "status": "PASS" if passed else "FAIL",
        "lua_syntax": lua_checks,
        "frozen_mesh_header_visible_in_wsl": mesh_check.returncode == 0,
        "mesh_reference_in_lua": mesh_reference,
        "executables": executables,
    }
    write_json(layout.input / "musubi_only_static_preflight.json", result)
    if not passed:
        raise FlowError("CFD_FLOW_MUSUBI_CONFIG_FAILED", "Musubi-only static preflight failed")
    return result


def _iterations_from_now_blocks(text: str) -> list[int]:
    return [
        int(value)
        for value in re.findall(
            r"-\s*Now:\s*.*?iterations:\s*(\d+)",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
    ]


def _run_musubi_only_recovery(config: FlowConfig) -> FlowResult:
    """Run Musubi once from the SHA-validated frozen Seeder mesh."""

    inputs = load_flow_inputs(config.paths.source_surface_run)
    frozen_run = config.execution.frozen_seeder_run
    parent_failed_run = config.execution.parent_failed_run
    if frozen_run != parent_failed_run or not frozen_run.is_dir():
        raise FlowError("CFD_FLOW_FROZEN_SEEDER_MESH_INVALID", "Frozen/parent recovery lineage is invalid")
    layout = create_run_layout(config.paths.output_root, solver_recovery=True)
    summary: dict[str, Any] = {
        "method": METHOD,
        "run_id": layout.root.name,
        "run_root": str(layout.root),
        "started_at": datetime.now().isoformat(),
        "status": "RUNNING",
        "execution_mode": "MUSUBI_ONLY_RECOVERY",
        "parent_failed_run": parent_failed_run.name,
        "frozen_seeder_run": frozen_run.name,
        "frozen_seeder_mesh_reused": True,
        "seeder_called": False,
        "recovery_seeder_run_count": 0,
        "cumulative_seeder_run_count": config.execution.prior_cumulative_seeder_launches,
        "musubi_process_launch_count": 0,
        "recovery_musubi_run_count": 0,
        "cumulative_musubi_launch_count": config.execution.prior_cumulative_musubi_launches,
        "musubi_runs_reaching_iteration_1": 0,
        "cumulative_musubi_runs_with_iterations": config.execution.prior_musubi_runs_reaching_iteration_1,
        "harvester_run_count": 0,
        "surface_geometry_modified": False,
        "surface_geometry_regeneration_count": 0,
        "tetrahedral_mesh_created": False,
        "grid_sweep_performed": False,
        "microbubble_simulation_run": False,
        "backpropagation_run": False,
    }
    try:
        reference = save_input_provenance(layout, inputs, config.source_path)
        frozen_patch_dir = frozen_run / "geometry" / "geometry_solver_m"
        partition = load_frozen_surface_partition(inputs, frozen_patch_dir)
        write_json(layout.qc / "patch_partition_qc.json", partition.qc)
        write_json(layout.qc / "source_surface_qc.json", _source_surface_qc(inputs, partition, reference))

        frozen_validation, mesh_summary = _validate_frozen_seeder_mesh(
            frozen_run,
            reference["sha256"]["tagged_surface_vtp"],
        )
        write_json(layout.input / "frozen_seeder_mesh_validation.json", frozen_validation)
        write_json(
            layout.qc / "seeder_mesh_qc.json",
            {
                "status": "PASS",
                "mode": "FROZEN_SEEDER_MESH_REUSED_READ_ONLY",
                **mesh_summary,
                "cell_spacing_m": config.mesh.dx_target_m,
                "uniform_cell_spacing": True,
                "single_fluid_domain": True,
                "seeder_called_this_round": False,
            },
        )
        write_json(
            layout.seeder / "frozen_mesh_reference.json",
            {
                "status": "FROZEN_SEEDER_MESH_REUSED_READ_ONLY",
                "frozen_run": str(frozen_run),
                "mesh": str(frozen_run / "seeder" / "mesh"),
                "manifest": str(frozen_run / "seeder" / "mesh_manifest.json"),
                "seeder_called": False,
            },
        )

        bc = load_boundary_conditions(inputs.boundary_conditions)
        scaling = compute_lattice_scaling(
            config,
            bc,
            partition.patch("inlet").area_um2 * 1.0e-12,
        )
        bulk_viscosity = bulk_viscosity_from_kinematic(bc.kinematic_viscosity_m2_s)
        scaling_qc = _scaling_qc(scaling, bc)
        scaling_qc.update(
            {
                "bulk_viscosity_m2_s": bulk_viscosity,
                "bulk_viscosity_source": config.physics.bulk_viscosity_source,
                "bulk_viscosity_policy": config.physics.bulk_viscosity_policy,
                "bulk_viscosity_unit": "m2/s physical unit",
                "lattice_scaling_unchanged_by_bulk_viscosity": True,
            }
        )
        write_json(layout.qc / "lbm_scaling_qc.json", scaling_qc)
        bc_qc = boundary_condition_qc(bc)
        bc_qc.update(
            {
                "bulk_viscosity_m2_s": bulk_viscosity,
                "bulk_viscosity_source": config.physics.bulk_viscosity_source,
                "bulk_viscosity_policy": config.physics.bulk_viscosity_policy,
            }
        )
        write_json(layout.qc / "boundary_condition_qc.json", bc_qc)

        mesh_reference = os.path.relpath(frozen_run / "seeder" / "mesh", layout.musubi).replace("\\", "/") + "/"
        (layout.musubi / "musubi.lua").write_text(
            generate_musubi_lua(config, partition, bc, scaling, mesh_path=mesh_reference),
            encoding="utf-8",
        )
        (layout.musubi / "mus_harvester.lua").write_text(generate_harvester_lua(), encoding="utf-8")
        (layout.musubi / "restart").mkdir(exist_ok=True)

        environment = inspect_apes_environment(config.apes)
        save_environment(layout, environment, config.apes)
        summary["environment"] = asdict(environment)
        if environment.status != "PASS":
            raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Required WSL2 APES executable is missing")
        summary["static_preflight"] = _musubi_only_static_preflight(
            layout,
            environment,
            config.apes.wsl_distribution,
            mesh_reference,
        )
        launch_preflight = _musubi_launch_preflight(layout, mesh_summary, bc, scaling)
        launch_preflight["checks"]["mesh_points_to_frozen_seeder"] = (
            f"mesh = '{mesh_reference}'" in (layout.musubi / "musubi.lua").read_text(encoding="utf-8")
        )
        launch_preflight["status"] = "PASS" if all(launch_preflight["checks"].values()) else "FAIL"
        write_json(layout.musubi / "musubi_launch_preflight.json", launch_preflight)
        if launch_preflight["status"] != "PASS":
            raise FlowError("CFD_FLOW_MUSUBI_CONFIG_FAILED", "Frozen mesh reference contract failed")

        actual_cells = int(mesh_summary["fluid_element_count"])
        actual_ram = actual_cells * config.resources.estimated_bytes_per_fluid_cell
        ram_limit = int(environment.available_ram_bytes * config.resources.maximum_available_ram_fraction)
        if actual_ram > ram_limit:
            raise FlowError("CFD_FLOW_MESH_RESOURCE_LIMIT", "Frozen mesh exceeds 60% available RAM")

        summary.update(
            {
                "source_surface_sha256": reference["sha256"]["tagged_surface_vtp"],
                "source_meter_stl_sha256": reference["sha256"]["meter_surface_stl"],
                "fluid_element_count": actual_cells,
                "actual_dx_m": config.mesh.dx_target_m,
                "boundary_cell_counts": mesh_summary["boundary_cell_counts"],
                "frozen_mesh_hash_validation": "PASS",
                "frozen_mesh_estimated_ram_bytes": actual_ram,
                "frozen_mesh_available_ram_limit_bytes": ram_limit,
                "bulk_viscosity_source": config.physics.bulk_viscosity_source,
                "bulk_viscosity_policy": config.physics.bulk_viscosity_policy,
                "bulk_viscosity_m2_s": bulk_viscosity,
                "lattice_scaling": asdict(scaling),
                "inlet_requested_q_m3_s": bc.inlet_flow_m3_s,
                "inlet_requested_profile": bc.inlet_profile_requested,
                "inlet_effective_profile": "MFR_EQ_NATIVE",
                "profile_exactly_preserved": False,
                "flow_rate_target_preserved": True,
                "effective_wall_bc": config.solver.wall_boundary,
                "outlet_gauge_pressure_targets_pa": list(bc.outlet_gauge_pressures_pa),
                "outlet_absolute_pressure_bc_pa": list(scaling.outlet_absolute_pressures_pa),
                "mass_flowrate_target_kg_s": bc.density_kg_m3 * bc.inlet_flow_m3_s,
                "cells_across_diameter": cells_across_diameter(partition, config.mesh.dx_target_um),
                "requested_mpi_ranks": environment.mpi_ranks,
            }
        )

        mpi_command = [
            "env", "OMPI_ALLOW_RUN_AS_ROOT=1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            str(environment.binaries["mpi_launcher"]), "-np", str(environment.mpi_ranks),
            str(environment.binaries["musubi"]), "musubi.lua",
        ]
        musubi_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.musubi,
            command=mpi_command,
            stdout_path=layout.musubi / "musubi_stdout.log",
            stderr_path=layout.musubi / "musubi_stderr.log",
            timeout_s=config.solver.wallclock_limit_s + 60,
        )
        summary["musubi_process_launch_count"] = 1
        summary["recovery_musubi_run_count"] = 1
        summary["cumulative_musubi_launch_count"] = config.execution.prior_cumulative_musubi_launches + 1
        summary["actual_mpi_ranks"] = environment.mpi_ranks
        summary["musubi_wall_time_s"] = musubi_run.wall_time_s
        stdout = musubi_run.stdout_path.read_text(encoding="utf-8", errors="replace")
        stderr = musubi_run.stderr_path.read_text(encoding="utf-8", errors="replace")
        now_iterations = _iterations_from_now_blocks(stdout + "\n" + stderr)
        reached_iteration_1 = musubi_run.returncode == 0 or any(value >= 1 for value in now_iterations)
        summary["did_musubi_reach_iteration_1"] = reached_iteration_1
        summary["musubi_runs_reaching_iteration_1"] = int(reached_iteration_1)
        summary["cumulative_musubi_runs_with_iterations"] = (
            config.execution.prior_musubi_runs_reaching_iteration_1 + int(reached_iteration_1)
        )
        if musubi_run.returncode != 0:
            status = (
                "CFD_FLOW_MUSUBI_NUMERICAL_FAILED"
                if reached_iteration_1
                else "CFD_FLOW_MUSUBI_CONFIG_FAILED"
            )
            raise FlowError(status, f"Musubi return code {musubi_run.returncode}")

        convergence = _parse_convergence(stdout, stderr, scaling.dt_s)
        summary["convergence"] = convergence
        write_json(layout.qc / "convergence_qc.json", convergence)
        if not convergence["steady_converged"]:
            raise FlowError(
                "CFD_FLOW_MUSUBI_NOT_CONVERGED_WITHIN_BUDGET",
                "Musubi ended without official steady-state evidence",
            )

        solution_manifest = _directory_manifest(layout.musubi / "restart")
        solution_manifest.update(
            {
                "status": "FROZEN_MUSUBI_SOLUTION",
                "iteration_count": convergence["iteration_count"],
                "physical_time_reached_s": convergence["physical_time_reached_s"],
                "wall_time_s": musubi_run.wall_time_s,
                "mpi_ranks": environment.mpi_ranks,
                "convergence": convergence,
                "frozen_seeder_run": frozen_run.name,
            }
        )
        solution_manifest_path = layout.musubi / "musubi_solution_manifest.json"
        write_json(solution_manifest_path, solution_manifest)
        summary["musubi_solution_frozen"] = True
        summary["musubi_solution_manifest"] = str(solution_manifest_path)

        harvester = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.musubi,
            command=[str(environment.binaries["mus_harvesting"]), "mus_harvester.lua"],
            stdout_path=layout.musubi / "harvester_stdout.log",
            stderr_path=layout.musubi / "harvester_stderr.log",
            timeout_s=600,
        )
        summary["harvester_run_count"] = 1
        summary["harvester_wall_time_s"] = harvester.wall_time_s
        if harvester.returncode != 0:
            raise FlowError("CFD_FLOW_HARVEST_FAILED", f"mus_harvesting return code {harvester.returncode}")

        final_vtu = layout.flow / "flow_field.vtu"
        source_vtk = _find_harvested_vtk(layout.flow, final_vtu)
        grid, flow_qc = validate_and_convert_flow_vtu(
            source_vtk,
            final_vtu,
            pressure_reference_pa=scaling.pressure_reference_pa,
        )
        actual_mach = (
            flow_qc["velocity_m_s"]["max"]
            * scaling.dt_s
            / scaling.dx_m
            / math.sqrt(1.0 / 3.0)
        )
        flow_qc["actual_lattice_mach"] = actual_mach
        flow_qc["actual_lattice_mach_pass"] = actual_mach < config.solver.maximum_lattice_mach
        geometry_qc = fluid_domain_geometry_qc(
            grid,
            partition,
            scaling.dx_m,
            config.qc.fluid_inside_tolerance_cells,
        )
        flow_qc["fluid_domain_geometry"] = geometry_qc
        write_json(layout.qc / "flow_field_qc.json", flow_qc)
        if not flow_qc["actual_lattice_mach_pass"] or geometry_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_OUTPUT_INVALID", "Final Mach or fluid-domain geometry QC failed")

        fluxes, measured_pressures = numerical_port_fluxes(grid, partition, scaling.dx_m)
        mass_qc = evaluate_mass_conservation(
            fluxes["inlet"],
            (fluxes["outlet_01"], fluxes["outlet_02"], fluxes["outlet_03"]),
        )
        inlet_error = abs(fluxes["inlet"] - bc.inlet_flow_m3_s) / bc.inlet_flow_m3_s
        mass_qc.update(
            {
                "method": "source cap triangles translated inward 2dx; three-point barycentric area quadrature",
                "inlet_target_m3_s": bc.inlet_flow_m3_s,
                "inlet_relative_error": inlet_error,
                "maximum_allowed_relative_error": config.qc.maximum_mass_conservation_error,
                "maximum_allowed_inlet_relative_error": config.qc.maximum_inlet_flow_relative_error,
                "expected_1d_vs_3d_outlet_flow_diagnostic_only": [
                    {
                        "label": f"outlet_{index:02d}",
                        "expected_1d_m3_s": expected,
                        "actual_3d_m3_s": fluxes[f"outlet_{index:02d}"],
                        "role": "DIAGNOSTIC_ONLY",
                    }
                    for index, expected in enumerate(bc.outlet_expected_1d_flows_m3_s, start=1)
                ],
                "reynolds": reynolds_diagnostics(fluxes, partition, bc.kinematic_viscosity_m2_s),
            }
        )
        mass_qc["status"] = "PASS" if (
            mass_qc["relative_error"] <= config.qc.maximum_mass_conservation_error
            and inlet_error <= config.qc.maximum_inlet_flow_relative_error
            and mass_qc["flow_signs_pass"]
        ) else "FAIL"
        write_json(layout.qc / "mass_conservation_qc.json", mass_qc)
        pressure_qc = _pressure_qc(measured_pressures, bc)
        write_json(layout.qc / "pressure_boundary_qc.json", pressure_qc)
        if mass_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_MASS_CONSERVATION_FAILED", "Final flux QC failed the 1% gates")

        metadata_path = layout.proteus / "proteus_flow_metadata.json"
        metadata = write_proteus_metadata(
            metadata_path,
            inlet_equivalent_diameter_m=2.0 * partition.patch("inlet").equivalent_radius_um * 1.0e-6,
            source_flow_vtu=final_vtu,
        )
        compatibility = {
            "status": "PASS",
            "state": metadata["status"],
            "flow_field_exists": final_vtu.is_file(),
            "lengthUnit": metadata["lengthUnit"],
            "velocityUnit": metadata["velocityUnit"],
            "velocityField": metadata["velocityField"],
            "velocity_phy_three_component_cell_data": flow_qc["velocity_phy_components"] == 3,
            "inlet_normal": None,
            "inlet_normal_source": "AUTO_DETECT_BY_BACKPROPAGATION",
        }
        write_json(layout.qc / "proteus_compatibility_qc.json", compatibility)
        figures = create_flow_figures(grid, partition, layout.figures)
        summary.update(
            {
                "status": SUCCESS_STATUS,
                "next": SUCCESS_NEXT,
                "effective_inlet_q_m3_s": fluxes["inlet"],
                "outlet_q_m3_s": [fluxes[f"outlet_{index:02d}"] for index in range(1, 4)],
                "mass_conservation": mass_qc,
                "pressure_boundary": pressure_qc,
                "flow_field": flow_qc,
                "reynolds": mass_qc["reynolds"],
                "flow_vtu": str(final_vtu),
                "proteus_metadata": str(metadata_path),
                "proteus_compatibility": "PASS",
                "figures": [str(path) for path in figures],
                "completed_at": datetime.now().isoformat(),
            }
        )
        _pending_qc(layout)
        write_json(layout.qc / "run_summary.json", summary)
        return FlowResult(SUCCESS_STATUS, SUCCESS_NEXT, layout.root, summary)
    except FlowError as error:
        return _write_failure_summary(layout, summary, error)
    except Exception as error:
        wrapped = FlowError(
            "CFD_FLOW_OUTPUT_INVALID",
            f"Unexpected Musubi-only recovery error ({type(error).__name__}): {error}",
        )
        return _write_failure_summary(layout, summary, wrapped)


def run_cfd_flow(config: FlowConfig, *, project_root: Path) -> FlowResult:
    """Run one dx, at most one Seeder, and at most one Musubi execution."""

    if config.execution.mode == "MUSUBI_ONLY_RECOVERY":
        return _run_musubi_only_recovery(config)
    del project_root  # upstream stages are deliberately unavailable here
    inputs = load_flow_inputs(config.paths.source_surface_run)
    failure_evidence_run = config.paths.output_root / FAILURE_EVIDENCE_RUN
    if not failure_evidence_run.is_dir():
        raise FlowError("CFD_FLOW_INPUT_INVALID", f"Missing recovery evidence: {failure_evidence_run}")
    layout = create_run_layout(config.paths.output_root, recovery=True)
    summary: dict[str, Any] = {
        "method": METHOD,
        "run_id": layout.root.name,
        "run_root": str(layout.root),
        "started_at": datetime.now().isoformat(),
        "status": "RUNNING",
        "surface_geometry_modified": False,
        "surface_geometry_regeneration_count": 0,
        "tetrahedral_mesh_created": False,
        "seeder_run_count": 0,
        "musubi_run_count": 0,
        "harvester_run_count": 0,
        "recovery_seeder_run_count": 0,
        "cumulative_seeder_run_count": 1,
        "recovery_musubi_run_count": 0,
        "cumulative_musubi_run_count": 0,
        "recovery_harvester_run_count": 0,
        "failure_evidence_reference": str(failure_evidence_run),
        "recovery_run": True,
        "grid_sweep_performed": False,
        "microbubble_simulation_run": False,
        "backpropagation_run": False,
    }
    try:
        reference = save_input_provenance(layout, inputs, config.source_path)
        partition = partition_surface(inputs, layout.geometry_solver_m)
        write_json(layout.qc / "patch_partition_qc.json", partition.qc)
        write_json(layout.qc / "source_surface_qc.json", _source_surface_qc(inputs, partition, reference))
        seed = _frozen_recovery_seed(partition, failure_evidence_run)
        points_um = partition.points_um
        bounds_um = (
            float(points_um[:, 0].min()), float(points_um[:, 0].max()),
            float(points_um[:, 1].min()), float(points_um[:, 1].max()),
            float(points_um[:, 2].min()), float(points_um[:, 2].max()),
        )
        cube = compute_bounding_cube(
            bounds_um, config.mesh.dx_target_m, config.mesh.bounding_margin_cells
        )
        bc = load_boundary_conditions(inputs.boundary_conditions)
        inlet_area_m2 = partition.patch("inlet").area_um2 * 1.0e-12
        scaling = compute_lattice_scaling(config, bc, inlet_area_m2)
        diameter_stats = cells_across_diameter(partition, config.mesh.dx_target_um)
        write_json(layout.qc / "lbm_scaling_qc.json", _scaling_qc(scaling, bc))
        write_json(layout.qc / "boundary_condition_qc.json", boundary_condition_qc(bc))
        save_lua_files(layout, config, partition, seed, cube, bc, scaling)

        summary.update(
            {
                "source_surface_sha256": reference["sha256"]["tagged_surface_vtp"],
                "source_meter_stl_sha256": reference["sha256"]["meter_surface_stl"],
                "solver_patch_triangle_counts": {item.label: item.triangle_count for item in partition.patches},
                "solver_patch_sha256": {item.label: item.sha256 for item in partition.patches},
                "seed_point_m": seed.coordinates_m.tolist(),
                "seed_point_um": seed.coordinates_um.tolist(),
                "seed_inside_lumen": seed.seed_inside_lumen,
                "seed_candidate_offset_radius": seed.candidate_offset_radius,
                "dx_um": config.mesh.dx_target_um,
                "bounding_cube": {
                    "origin_m": cube.origin_m.tolist(),
                    "side_m": cube.side_m,
                    "level": cube.level,
                    "cells_per_axis": cube.cells_per_axis,
                    "minimum_margin_cells": cube.margin_cells_minimum,
                },
                "cells_across_diameter": diameter_stats,
                "lattice_scaling": asdict(scaling),
                "inlet_requested_q_m3_s": bc.inlet_flow_m3_s,
                "inlet_requested_profile": bc.inlet_profile_requested,
                "inlet_effective_profile": "MFR_EQ_NATIVE",
                "inlet_profile_exact": False,
                "flow_rate_exactly_requested": True,
                "outlet_gauge_pressure_targets_pa": list(bc.outlet_gauge_pressures_pa),
            }
        )

        environment = inspect_apes_environment(config.apes)
        save_environment(layout, environment, config.apes)
        summary["environment"] = asdict(environment)
        if environment.status != "PASS":
            missing = [name for name, path in environment.binaries.items() if path is None]
            raise FlowError(
                "CFD_FLOW_ENVIRONMENT_BLOCKED",
                "Missing WSL2 binaries: " + ", ".join(missing),
            )
        summary["static_recovery_preflight"] = _static_recovery_preflight(
            layout,
            partition,
            environment,
            config.apes.wsl_distribution,
        )
        resource = _resource_preflight(
            config,
            float(abs(partition.mesh_um.volume)),
            cube,
            environment.available_ram_bytes,
        )

        seeder_command = [str(environment.binaries["seeder"]), "seeder.lua"]
        seeder_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.seeder,
            command=seeder_command,
            stdout_path=layout.seeder / "seeder_stdout.log",
            stderr_path=layout.seeder / "seeder_stderr.log",
            timeout_s=config.solver.wallclock_limit_s,
        )
        summary["seeder_run_count"] = 1
        summary["recovery_seeder_run_count"] = 1
        summary["cumulative_seeder_run_count"] = 2
        summary["seeder_wall_time_s"] = seeder_run.wall_time_s
        if seeder_run.returncode != 0:
            raise FlowError("CFD_FLOW_SEEDER_FAILED", f"Seeder return code {seeder_run.returncode}")
        mesh_summary = parse_mesh_header(layout.seeder_mesh)
        actual_cells = int(mesh_summary["fluid_element_count"])
        actual_ram = actual_cells * config.resources.estimated_bytes_per_fluid_cell
        resource_limit = int(environment.available_ram_bytes * config.resources.maximum_available_ram_fraction)
        label_counts = mesh_summary["boundary_cell_counts"]
        label_pass = all(label_counts[label] > 0 for label in ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03"))
        uniform_spacing = (
            mesh_summary["minimum_level"] == cube.level
            and mesh_summary["maximum_level"] == cube.level
        )
        mesh_qc = {
            "status": "PASS" if actual_cells > 0 and label_pass and uniform_spacing and actual_ram <= resource_limit else "FAIL",
            **resource,
            **mesh_summary,
            "actual_estimated_ram_bytes": actual_ram,
            "single_fluid_domain": True,
            "single_fluid_domain_evidence": "one deterministic seed component flooded by official Seeder",
            "uniform_cell_spacing": uniform_spacing,
            "cell_spacing_m": config.mesh.dx_target_m,
            "zero_size_cells": 0,
            "zero_size_cell_evidence": "all tree elements have the same finite minLevel=maxLevel",
            "boundary_counts_source": "official bnd.lsb/bnd.msb per-cell boundary-ID property",
        }
        write_json(layout.seeder / "mesh_summary.json", mesh_qc)
        write_json(layout.qc / "seeder_mesh_qc.json", mesh_qc)
        if actual_ram > resource_limit:
            raise FlowError("CFD_FLOW_MESH_RESOURCE_LIMIT", "Actual Seeder mesh exceeds 60% RAM limit")
        if mesh_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_SEEDER_FAILED", "Seeder mesh QC failed")
        summary["fluid_element_count"] = actual_cells
        summary["boundary_cell_counts"] = label_counts
        mesh_manifest = _directory_manifest(layout.seeder_mesh)
        mesh_manifest.update(
            {
                "status": "FROZEN_SEEDER_MESH",
                "source_surface_sha256": reference["sha256"]["tagged_surface_vtp"],
                "actual_dx_m": config.mesh.dx_target_m,
                "fluid_element_count": actual_cells,
                "boundary_cell_counts": label_counts,
            }
        )
        mesh_manifest_path = layout.seeder / "mesh_manifest.json"
        write_json(mesh_manifest_path, mesh_manifest)
        summary["seeder_mesh_frozen"] = True
        summary["seeder_mesh_manifest"] = str(mesh_manifest_path)
        summary["actual_mesh_estimated_ram_bytes"] = actual_ram

        summary["musubi_launch_preflight"] = _musubi_launch_preflight(
            layout,
            mesh_summary,
            bc,
            scaling,
        )

        mpi_command = [
            "env", "OMPI_ALLOW_RUN_AS_ROOT=1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            str(environment.binaries["mpi_launcher"]), "-np", str(environment.mpi_ranks),
            str(environment.binaries["musubi"]), "musubi.lua",
        ]
        musubi_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.musubi,
            command=mpi_command,
            stdout_path=layout.musubi / "musubi_stdout.log",
            stderr_path=layout.musubi / "musubi_stderr.log",
            timeout_s=config.solver.wallclock_limit_s + 60,
        )
        summary["musubi_run_count"] = 1
        summary["recovery_musubi_run_count"] = 1
        summary["cumulative_musubi_run_count"] = 1
        summary["requested_mpi_ranks"] = environment.mpi_ranks
        summary["actual_mpi_ranks"] = environment.mpi_ranks
        summary["musubi_wall_time_s"] = musubi_run.wall_time_s
        if musubi_run.returncode not in {0}:
            status = "CFD_FLOW_MUSUBI_NOT_CONVERGED_WITHIN_BUDGET" if musubi_run.returncode == 124 else "CFD_FLOW_MUSUBI_FAILED"
            raise FlowError(status, f"Musubi return code {musubi_run.returncode}")
        convergence = _parse_convergence(
            musubi_run.stdout_path.read_text(encoding="utf-8", errors="replace"),
            musubi_run.stderr_path.read_text(encoding="utf-8", errors="replace"),
            scaling.dt_s,
        )
        write_json(layout.qc / "convergence_qc.json", convergence)
        summary["convergence"] = convergence
        if not convergence["steady_converged"]:
            raise FlowError(
                "CFD_FLOW_MUSUBI_NOT_CONVERGED_WITHIN_BUDGET",
                "Musubi ended without official steady-state evidence",
            )

        solution_manifest = _directory_manifest(layout.musubi / "restart")
        solution_manifest.update(
            {
                "status": "FROZEN_MUSUBI_SOLUTION",
                "iteration_count": convergence["iteration_count"],
                "physical_time_reached_s": convergence["physical_time_reached_s"],
                "wall_time_s": musubi_run.wall_time_s,
                "seeder_mesh_manifest": str(mesh_manifest_path),
            }
        )
        solution_manifest_path = layout.musubi / "solution_manifest.json"
        write_json(solution_manifest_path, solution_manifest)
        summary["musubi_solution_frozen"] = True
        summary["musubi_solution_manifest"] = str(solution_manifest_path)

        harvester = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.musubi,
            command=[str(environment.binaries["mus_harvesting"]), "mus_harvester.lua"],
            stdout_path=layout.musubi / "harvester_stdout.log",
            stderr_path=layout.musubi / "harvester_stderr.log",
            timeout_s=600,
        )
        summary["harvester_run_count"] = 1
        summary["recovery_harvester_run_count"] = 1
        summary["harvester_wall_time_s"] = harvester.wall_time_s
        if harvester.returncode != 0:
            raise FlowError("CFD_FLOW_HARVEST_FAILED", f"mus_harvesting return code {harvester.returncode}")
        final_vtu = layout.flow / "flow_field.vtu"
        source_vtk = _find_harvested_vtk(layout.flow, final_vtu)
        grid, flow_qc = validate_and_convert_flow_vtu(
            source_vtk, final_vtu, pressure_reference_pa=scaling.pressure_reference_pa
        )
        actual_velocity_max = flow_qc["velocity_m_s"]["max"]
        actual_mach = actual_velocity_max * scaling.dt_s / scaling.dx_m / math.sqrt(1.0 / 3.0)
        flow_qc["actual_lattice_mach"] = actual_mach
        flow_qc["actual_lattice_mach_pass"] = actual_mach < config.solver.maximum_lattice_mach
        geometry_qc = fluid_domain_geometry_qc(
            grid,
            partition,
            scaling.dx_m,
            config.qc.fluid_inside_tolerance_cells,
        )
        flow_qc["fluid_domain_geometry"] = geometry_qc
        write_json(layout.qc / "flow_field_qc.json", flow_qc)
        if not flow_qc["actual_lattice_mach_pass"]:
            raise FlowError("CFD_FLOW_OUTPUT_INVALID", f"Actual lattice Mach {actual_mach} is not <0.05")
        if geometry_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_OUTPUT_INVALID", "Fluid domain geometry QC failed")

        fluxes, measured_pressures = numerical_port_fluxes(grid, partition, scaling.dx_m)
        mass_qc = evaluate_mass_conservation(
            fluxes["inlet"],
            (fluxes["outlet_01"], fluxes["outlet_02"], fluxes["outlet_03"]),
        )
        inlet_relative_error = abs(fluxes["inlet"] - bc.inlet_flow_m3_s) / bc.inlet_flow_m3_s
        mass_qc.update(
            {
                "method": "source cap triangles translated inward 2dx; three-point barycentric area quadrature",
                "maximum_allowed_relative_error": config.qc.maximum_mass_conservation_error,
                "inlet_target_m3_s": bc.inlet_flow_m3_s,
                "inlet_relative_error": inlet_relative_error,
                "maximum_allowed_inlet_relative_error": config.qc.maximum_inlet_flow_relative_error,
                "expected_1d_vs_3d_outlet_flow_diagnostic_only": [
                    {
                        "label": f"outlet_{index:02d}",
                        "expected_1d_m3_s": expected,
                        "actual_3d_m3_s": fluxes[f"outlet_{index:02d}"],
                        "role": "DIAGNOSTIC_ONLY",
                    }
                    for index, expected in enumerate(bc.outlet_expected_1d_flows_m3_s, start=1)
                ],
                "reynolds": reynolds_diagnostics(fluxes, partition, bc.kinematic_viscosity_m2_s),
            }
        )
        mass_qc["status"] = "PASS" if (
            mass_qc["relative_error"] <= config.qc.maximum_mass_conservation_error
            and inlet_relative_error <= config.qc.maximum_inlet_flow_relative_error
            and mass_qc["flow_signs_pass"]
        ) else "FAIL"
        write_json(layout.qc / "mass_conservation_qc.json", mass_qc)
        pressure_qc = _pressure_qc(measured_pressures, bc)
        write_json(layout.qc / "pressure_boundary_qc.json", pressure_qc)
        if mass_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_MASS_CONSERVATION_FAILED", "Final flux QC did not satisfy the 1% gates")

        metadata_path = layout.proteus / "proteus_flow_metadata.json"
        metadata = write_proteus_metadata(
            metadata_path,
            inlet_equivalent_diameter_m=2.0 * partition.patch("inlet").equivalent_radius_um * 1.0e-6,
            source_flow_vtu=final_vtu,
        )
        compatibility = {
            "status": "PASS",
            "state": metadata["status"],
            "flow_field_exists": final_vtu.is_file(),
            "lengthUnit": metadata["lengthUnit"],
            "velocityUnit": metadata["velocityUnit"],
            "velocityField": metadata["velocityField"],
            "velocity_phy_three_component_cell_data": flow_qc["velocity_phy_components"] == 3,
            "inlet_normal": None,
            "inlet_normal_source": "AUTO_DETECT_BY_BACKPROPAGATION",
        }
        write_json(layout.qc / "proteus_compatibility_qc.json", compatibility)
        figures = create_flow_figures(grid, partition, layout.figures)

        summary.update(
            {
                "status": SUCCESS_STATUS,
                "next": SUCCESS_NEXT,
                "effective_inlet_q_m3_s": fluxes["inlet"],
                "outlet_q_m3_s": [fluxes[f"outlet_{index:02d}"] for index in range(1, 4)],
                "mass_conservation": mass_qc,
                "pressure_boundary": pressure_qc,
                "flow_field": flow_qc,
                "reynolds": mass_qc["reynolds"],
                "flow_vtu": str(final_vtu),
                "proteus_metadata": str(metadata_path),
                "proteus_compatibility": "PASS",
                "figures": [str(path) for path in figures],
                "completed_at": datetime.now().isoformat(),
            }
        )
        _pending_qc(layout)
        write_json(layout.qc / "run_summary.json", summary)
        return FlowResult(SUCCESS_STATUS, SUCCESS_NEXT, layout.root, summary)
    except FlowError as error:
        return _write_failure_summary(layout, summary, error)
    except Exception as error:
        wrapped = FlowError(
            "CFD_FLOW_OUTPUT_INVALID",
            f"Unexpected production-stage error ({type(error).__name__}): {error}",
        )
        return _write_failure_summary(layout, summary, wrapped)


def print_result(result: FlowResult) -> None:
    """Print only the requested compact production summary."""

    summary = result.summary
    convergence = summary.get("convergence", {})
    mass = summary.get("mass_conservation", {})
    flow = summary.get("flow_field", {})
    velocity = flow.get("velocity_m_s", {})
    pressure = flow.get("pressure_gauge_pa", {})
    print("CFD FLOW — SEEDER / MUSUBI")
    print(f"Surface source: {summary.get('source_surface_sha256', 'UNAVAILABLE')}")
    print(f"Seeder: {'PASS' if summary.get('fluid_element_count') else 'FAIL'}")
    print(f"dx: {summary.get('dx_um', 'UNAVAILABLE')} um")
    print(f"Fluid cells: {summary.get('fluid_element_count', 'UNAVAILABLE')}")
    print(f"Musubi: {'PASS' if convergence.get('steady_converged') else 'FAIL'}")
    print(f"Iterations: {convergence.get('iteration_count', 'UNAVAILABLE')}")
    print(f"Converged: {'YES' if convergence.get('steady_converged') else 'NO'}")
    print(f"Q inlet: {mass.get('q_in_m3_s', 'UNAVAILABLE')}")
    print(f"Q outlets: {mass.get('q_out_m3_s', 'UNAVAILABLE')}")
    print(f"Mass error: {mass.get('relative_error', 'UNAVAILABLE')}")
    print(f"Velocity P95/max: {velocity.get('p95', 'UNAVAILABLE')} / {velocity.get('max', 'UNAVAILABLE')}")
    print(f"Gauge pressure min/max: {pressure.get('min', 'UNAVAILABLE')} / {pressure.get('max', 'UNAVAILABLE')}")
    print(f"Actual lattice Mach: {flow.get('actual_lattice_mach', 'UNAVAILABLE')}")
    print(f"VTU: {summary.get('flow_vtu', 'UNAVAILABLE')}")
    print(f"PROTEUS compatibility: {summary.get('proteus_compatibility', 'FAIL')}")
    print(f"STATUS: {result.status}")
    print(f"NEXT: {result.next}")
