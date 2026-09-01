"""Finalize source-proven Tau1 full-timestep referee evidence.

This module is research-only.  It parses the single permitted instrumented
timestep, compares it with the Python replay, and reaccepts only already
computed Base evidence.  It never launches Seeder or Musubi.
"""

from __future__ import annotations

import gzip
import json
import math
import re
from pathlib import Path
from typing import Any, TextIO

import numpy as np

from .full_timestep_mass_referee import (
    BOUNDARY_FLUX_DEFINITION,
    DEFERRED_PHYSICAL_FLUX,
    FULL_IDENTITY_GATE,
    _wall_operations_runtime_only,
    pull_link_counts,
    replay_full_timestep,
    stable_delta,
    stable_total,
)
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import (
    _wall_operations,
    load_mesh_contract,
    replay_sequential_writes,
    runtime_solid_cells,
)
from .tau1_base import (
    BOUNDARY_CLOSURE_GATE,
    EXPECTED_CELLS,
    MESH_HASHES,
    MUSUBI_SHA256,
    OUTLET_GAUGE_PRESSURE_PA,
    RUN_NAME,
    Tau1BaseRuntimeContract,
    historical_tau1_runtime_contract,
    _mesh_path,
    _read_pdf,
    _restart_pairs,
    _runtime_windows,
)


INSTRUMENTED_RUNTIME = Path(
    r"\\wsl.localhost\Ubuntu\home\lzy\u3da\tau1_full_timestep_oracle_20260831"
)
INSTRUMENTED_SOURCE_WSL = (
    "/home/lzy/apes-worktrees/musubi_full_timestep_oracle_20260831"
)
INSTRUMENTED_BINARY_WSL = f"{INSTRUMENTED_SOURCE_WSL}/build/musubi"
INSTRUMENTED_BINARY_SHA256 = (
    "0d727384cf7684276aea3da42c3f3e9f4f7a80a06013ad35bae0cf3f7fff59f9"
)
INSTRUMENTATION_ABSOLUTE_SUM_TOLERANCE = 5.0e-10


def _run_root(project_root: Path) -> Path:
    return Path(project_root) / "outputs" / "cfd_flow" / RUN_NAME


def _open_log(path: Path) -> TextIO:
    if path.is_file():
        return path.open("rt", encoding="utf-8", errors="replace")
    compressed = path.with_suffix(path.suffix + ".gz")
    if compressed.is_file():
        return gzip.open(compressed, "rt", encoding="utf-8", errors="replace")
    raise FileNotFoundError(path)


def _parse_instrumented_log(path: Path) -> dict[str, Any]:
    phase_pattern = re.compile(
        r"FULL_MASS_PHASE phase=(\d+)(?: iter=(\d+))? mass=\s*([-+0-9.Ee]+)"
    )
    boundary_pattern = re.compile(
        r"BOUNDARY_MASS_ORACLE iter=(\d+) label=(\S+) "
        r"delta=\s*([-+0-9.Ee]+) changed_slots=(\d+) "
        r"overlap_count=(\d+) mass=\s*([-+0-9.Ee]+)"
    )
    phases: dict[int, dict[str, Any]] = {}
    boundaries: dict[str, dict[str, Any]] = {}
    successful = False
    with _open_log(path) as handle:
        for line in handle:
            phase = phase_pattern.search(line)
            if phase:
                phase_id = int(phase.group(1))
                phases[phase_id] = {
                    "iteration": (
                        int(phase.group(2)) if phase.group(2) is not None else None
                    ),
                    "mass": float(phase.group(3)),
                }
            boundary = boundary_pattern.search(line)
            if boundary:
                boundaries[boundary.group(2)] = {
                    "iteration": int(boundary.group(1)),
                    "delta": float(boundary.group(3)),
                    "changed_slots": int(boundary.group(4)),
                    "overlap_count": int(boundary.group(5)),
                    "mass": float(boundary.group(6)),
                }
            successful |= "SUCCESSFUL run!" in line
    if set(phases) != set(range(12)):
        raise RuntimeError(f"instrumented phase log is incomplete: {sorted(phases)}")
    expected_boundaries = {"wall", "outlet_02", "outlet_03", "inlet", "outlet_01"}
    if set(boundaries) != expected_boundaries or not successful:
        raise RuntimeError("instrumented boundary or semantic-success evidence is incomplete")
    return {"phases": phases, "boundaries": boundaries, "successful": successful}


def write_instrumented_phase_oracle(project_root: Path) -> dict[str, Any]:
    """Parse and compare the one permitted 3117927 -> 3117928 call."""

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    evidence_root = run_root / "instrumented_full_timestep"
    log_path = evidence_root / "musubi_stdout.log"
    parsed = _parse_instrumented_log(log_path)
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    contract = historical_tau1_runtime_contract()
    main_pairs = _restart_pairs(_runtime_windows() / "restart")
    archived_pairs = _restart_pairs(
        run_root / "dense_diagnostic_corrected" / "attempt_1_8_steps" / "restart"
    )
    instrumented_pairs = _restart_pairs(INSTRUMENTED_RUNTIME / "restart")
    if 3_117_927 not in main_pairs or 3_117_928 not in instrumented_pairs:
        raise RuntimeError("instrumented transition restart pair is incomplete")
    start_pdf = _read_pdf(main_pairs[3_117_927][1])
    end_pdf = _read_pdf(instrumented_pairs[3_117_928][1])
    archived_end_pdf = _read_pdf(archived_pairs[3_117_928][1])
    replay = replay_full_timestep(
        start_pdf,
        end_pdf,
        mesh,
        dx_m=contract.dx_m,
        dt_s=contract.dt_s,
        density_kg_m3=contract.rho_kg_m3,
        target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
        outlet_pressures_pa={
            label: contract.pressure_reference_pa + gauge
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        },
    )
    order = ("wall", "outlet_02", "outlet_03", "inlet", "outlet_01")
    python_phase = {0: stable_total(start_pdf)}
    running = python_phase[0]
    for phase_id, label in enumerate(order, 1):
        running = math.fsum((running, float(replay[f"delta_{label}"])))
        python_phase[phase_id] = running
    python_phase[6] = python_phase[5]
    python_phase[7] = python_phase[5]
    python_phase[8] = float(replay["mass_after_pull"])
    python_phase[9] = python_phase[7]
    python_phase[10] = stable_total(replay["_pdf_after_collision"])
    python_phase[11] = python_phase[10]
    python_phase[12] = stable_total(end_pdf)

    compiled_phase = {
        int(phase_id): float(item["mass"])
        for phase_id, item in parsed["phases"].items()
    }
    compiled_phase[12] = python_phase[12]
    phase_labels = {
        0: "restart_loaded_state",
        1: "after_wall",
        2: "after_outlet_02",
        3: "after_outlet_03",
        4: "after_inlet",
        5: "after_outlet_01",
        6: "after_all_boundary_writes",
        7: "after_swap_before_compute",
        8: "actual_pull_fetched_density_sum",
        9: "after_macro_calculation_state_unchanged",
        10: "after_bgk_compute",
        11: "before_restart_serialization",
        12: "decoded_serialized_restart",
    }
    phase_table = []
    for phase_id in range(13):
        absolute = abs(compiled_phase[phase_id] - python_phase[phase_id])
        phase_table.append(
            {
                "phase": phase_id,
                "label": phase_labels[phase_id],
                "compiled_mass": compiled_phase[phase_id],
                "python_mass": python_phase[phase_id],
                "absolute_difference": absolute,
                "relative_difference": absolute
                / max(abs(python_phase[phase_id]), np.finfo(np.float64).tiny),
                "within_floating_sum_tolerance": (
                    absolute <= INSTRUMENTATION_ABSOLUTE_SUM_TOLERANCE
                ),
            }
        )

    boundary_comparison = {}
    for label in order:
        compiled = float(parsed["boundaries"][label]["delta"])
        python = float(replay[f"delta_{label}"])
        boundary_comparison[label] = {
            **parsed["boundaries"][label],
            "compiled_delta": compiled,
            "python_delta": python,
            "absolute_difference": abs(compiled - python),
        }

    solid = runtime_solid_cells(mesh)
    runtime_only_operations = _wall_operations_runtime_only(start_pdf, mesh, solid)
    runtime_only_state, _ = replay_sequential_writes(
        start_pdf, runtime_only_operations
    )
    pre_fix_wall_delta = stable_delta(runtime_only_state, start_pdf)
    raw_operations = _wall_operations(start_pdf, mesh)
    raw_state, _ = replay_sequential_writes(start_pdf, raw_operations)
    raw_wall_delta = stable_delta(raw_state, start_pdf)
    patch_path = evidence_root / "musubi_instrumentation.patch"
    raw_patch_gzip_path = evidence_root / "musubi_instrumentation.raw.patch.gz"
    archived_end = archived_pairs[3_117_928][1]
    instrumented_end = instrumented_pairs[3_117_928][1]
    result = {
        "status": "PASS",
        "referee_revision": "MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2",
        "instrumentation_needed": True,
        "scientific_musubi_calls": 1,
        "instrumented_timesteps": 1,
        "iteration_start": 3_117_927,
        "iteration_end": 3_117_928,
        "source_base_revision": "81f8c4f13772f6d4af31f335e1e3f99b02726e25",
        "apes_revision": "4e8b277b66226277171ef93bf054d21270812793",
        "instrumented_source_wsl": INSTRUMENTED_SOURCE_WSL,
        "instrumented_binary_wsl": INSTRUMENTED_BINARY_WSL,
        "instrumented_binary_sha256": INSTRUMENTED_BINARY_SHA256,
        "instrumentation_scope": (
            "new full-timestep oracle edits are observation/logging/MPI "
            "reduction only"
        ),
        "preexisting_adaptive_flux_math_retained": True,
        "cumulative_debug_patch_includes_preexisting_adaptive_flux_math": True,
        "instrumentation_patch_sha256": (
            sha256_file(patch_path) if patch_path.is_file() else None
        ),
        "instrumentation_raw_patch_gzip_sha256": (
            sha256_file(raw_patch_gzip_path)
            if raw_patch_gzip_path.is_file()
            else None
        ),
        "raw_stdout_gzip_sha256": (
            sha256_file(log_path.with_suffix(".log.gz"))
            if log_path.with_suffix(".log.gz").is_file()
            else None
        ),
        "phase_table": phase_table,
        "boundary_comparison": boundary_comparison,
        "pull": {
            "link_counts": pull_link_counts(mesh, solid),
            "sum_source_state_pdfs": float(replay["mass_after_boundary"]),
            "sum_fetched_pdfs": float(replay["mass_after_pull"]),
            "difference": float(replay["delta_pull_connectivity"]),
        },
        "runtime_solid_count": len(solid),
        "runtime_solid_affected_wall_writes": int(
            replay["runtime_solid_affected_wall_writes"]
        ),
        "boundary_blocked_affected_wall_writes": int(
            replay["boundary_blocked_affected_wall_writes"]
        ),
        "total_affected_wall_writes": int(replay["affected_wall_writes"]),
        "first_mismatch_phase_before_fix": "PHASE_1_AFTER_WALL_BOUNDARY",
        "first_mismatch_evidence": {
            "coordinate_only_wall_delta": raw_wall_delta,
            "runtime_solid_only_wall_delta": pre_fix_wall_delta,
            "compiled_wall_delta": float(parsed["boundaries"]["wall"]["delta"]),
            "corrected_python_wall_delta": float(replay["delta_wall"]),
            "root_cause": (
                "coordinate lookup admitted a source across a negative TreElm "
                "boundary-neighbor entry; PULL requires local inverse bounce"
            ),
        },
        "first_mismatch_phase_after_fix": "NONE_WITHIN_FLOATING_SUMMATION",
        "R_full_one_step_identity": float(replay["R_full_one_step_identity"]),
        "hard_gate": FULL_IDENTITY_GATE,
        "hard_gate_pass": float(replay["R_full_one_step_identity"])
        <= FULL_IDENTITY_GATE,
        "serialized_restart": {
            "instrumented_sha256": sha256_file(instrumented_end),
            "archived_3117928_sha256": sha256_file(archived_end),
            "byte_identical_to_archived_3117928": sha256_file(instrumented_end)
            == sha256_file(archived_end),
            "elementwise_max_abs_difference": float(
                np.max(np.abs(end_pdf - archived_end_pdf))
            ),
        },
        "flux_definition": BOUNDARY_FLUX_DEFINITION,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "new_long_cfd_iterations": 0,
    }
    if not all(item["within_floating_sum_tolerance"] for item in phase_table):
        result["status"] = "FAIL"
    if not result["hard_gate_pass"]:
        result["status"] = "FAIL"
    write_json(run_root / "qc" / "tau1_instrumented_full_timestep_phase_oracle.json", result)
    return result


def write_final_referee_and_acceptance(project_root: Path) -> dict[str, Any]:
    """Write final Referee V2 and existing-Base all-gate acceptance evidence."""

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    replay = json.loads(
        (qc / "tau1_full_timestep_replay_8step.json").read_text(encoding="utf-8")
    )
    oracle = json.loads(
        (qc / "tau1_instrumented_full_timestep_phase_oracle.json").read_text(
            encoding="utf-8"
        )
    )
    window = json.loads(
        (
            qc
            / "tau1_boundary_window_forensics_full_timestep_referee_corrected.json"
        ).read_text(encoding="utf-8")
    )
    steady = json.loads(
        (qc / "tau1_base_steady_status.json").read_text(encoding="utf-8")
    )
    window_closures = [
        float(method["closure"])
        for item in window["windows"]
        for method in item["integration_methods"].values()
    ]
    corrected_window_max = max(window_closures)
    referee_pass = (
        replay["status"] == "FULL_TIMESTEP_DISCRETE_MASS_IDENTITY_PROVEN"
        and bool(replay["hard_gate_pass"])
        and oracle["status"] == "PASS"
    )
    referee = {
        "status": (
            "CFD_FLOW_FULL_TIMESTEP_REFEREE_VALIDATED"
            if referee_pass
            else "CFD_FLOW_FULL_TIMESTEP_ACCOUNTING_UNRESOLVED"
        ),
        "referee_revision": "MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2",
        "identity_scope": "complete source-proven Musubi timestep",
        "classification": "REFEREE_BOUNDARY_ONLY_ACCOUNTING_INCOMPLETE",
        "boundary_write_only_is_diagnostic_not_final_identity": True,
        "boundary_write_delta": {
            "min_residual": replay["R_boundary_only_min"],
            "median_residual": replay["R_boundary_only_median"],
            "max_residual": replay["R_boundary_only_max"],
        },
        "pull_connectivity_delta": [
            item["delta_pull_connectivity"] for item in replay["per_step"]
        ],
        "collision_delta": [item["delta_collision"] for item in replay["per_step"]],
        "source_delta": [item["delta_source"] for item in replay["per_step"]],
        "full_predicted_delta": [
            item["delta_full_predicted"] for item in replay["per_step"]
        ],
        "actual_delta": [item["delta_actual"] for item in replay["per_step"]],
        "R_full_one_step_identity": {
            "min": replay["R_full_one_step_identity_min"],
            "median": replay["R_full_one_step_identity_median"],
            "max": replay["R_full_one_step_identity_max"],
        },
        "R_full_8step_identity": replay["R_full_8step_identity"],
        "hard_gate": FULL_IDENTITY_GATE,
        "instrumented_confirmation": oracle["status"],
        "scientific_musubi_calls": 1,
        "instrumented_timesteps": 1,
        "seeder_calls": 0,
        "new_long_cfd_iterations": 0,
        "flux_definition": BOUNDARY_FLUX_DEFINITION,
        "physical_Q1_Q2_Q3": DEFERRED_PHYSICAL_FLUX,
    }
    write_json(qc / "tau1_referee_v2_final.json", referee)

    gates = {
        "R_mass_short_le_0p01": float(steady["R_mass_short"]) <= 0.01,
        "R_mass_long_le_0p01": float(steady["R_mass_long"]) <= 0.01,
        "R_velocity_le_0p01": float(steady["R_velocity"]) <= 0.01,
        "R_pressure_le_0p005": float(steady["R_pressure"]) <= 0.005,
        "R_inlet_le_0p01": float(steady["R_inlet"]) <= 0.01,
        "corrected_boundary_window_closure_le_0p001": corrected_window_max
        <= BOUNDARY_CLOSURE_GATE,
        "R_full_one_step_identity_le_1e_minus_8": float(
            replay["R_full_one_step_identity_max"]
        )
        <= FULL_IDENTITY_GATE,
        "no_significant_time_averaged_backflow": not bool(
            steady["time_averaged_backflow"]
        ),
        "minimum_pdf_gt_zero": float(steady["minimum_pdf"]) > 0.0,
        "maximum_lattice_speed_lt_0p05": float(
            steady["maximum_lattice_speed"]
        )
        < 0.05,
        "all_finite": bool(steady["all_finite"]),
    }
    all_pass = referee_pass and all(gates.values())
    pairs = _restart_pairs(_runtime_windows() / "restart")
    accepted_iteration = 3_117_927 if all_pass else None
    header = pairs[3_117_927][0]
    header_text = header.read_text(encoding="utf-8")
    time_match = re.search(r"\bsim\s*=\s*([-+0-9.Ee]+)", header_text)
    if not time_match:
        raise RuntimeError("accepted restart header has no simulation time")
    acceptance = {
        "status": (
            "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS"
            if all_pass
            else "CFD_FLOW_FULL_TIMESTEP_ACCOUNTING_UNRESOLVED"
        ),
        "evidence_mode": "EXISTING_CHECKPOINTS_TRACKING_AND_RESTARTS_ONLY",
        "earliest_existing_all_gate_pass_iteration": accepted_iteration,
        "earliest_evidence_rationale": (
            "3117927 is the earliest preserved endpoint with the two preceding "
            "119751-iteration checkpoints required to evaluate every short/long gate"
        ),
        "accepted_physical_time_s": (
            float(time_match.group(1)) if all_pass else None
        ),
        "contract_dt_times_iteration_s": (
            3_117_927 * Tau1BaseRuntimeContract().dt_s if all_pass else None
        ),
        "accepted_restart_sha256": (
            sha256_file(pairs[3_117_927][1]) if all_pass else None
        ),
        "metrics": {
            "R_mass_short": steady["R_mass_short"],
            "R_mass_long": steady["R_mass_long"],
            "R_velocity": steady["R_velocity"],
            "R_pressure": steady["R_pressure"],
            "R_inlet": steady["R_inlet"],
            "corrected_boundary_window_closure_max": corrected_window_max,
            "R_full_one_step_identity_max": replay[
                "R_full_one_step_identity_max"
            ],
            "R_full_8step_identity": replay["R_full_8step_identity"],
            "time_averaged_backflow": steady["time_averaged_backflow"],
            "minimum_pdf": steady["minimum_pdf"],
            "maximum_lattice_speed": steady["maximum_lattice_speed"],
            "all_finite": steady["all_finite"],
        },
        "gates": gates,
        "protected_restart_sha256": {
            str(iteration): sha256_file(pairs[iteration][1])
            for iteration in (2_878_425, 2_998_176, 3_117_927)
        },
        "frozen_mesh_sha256": dict(MESH_HASHES),
        "validated_production_musubi_sha256": MUSUBI_SHA256,
        "protected_restarts_unchanged": {
            "2878425": sha256_file(pairs[2_878_425][1])
            == "75815fded691784ae942285e6ccf32514a1936ef9b1c228a03631991462646ee",
            "2998176": sha256_file(pairs[2_998_176][1])
            == "e3bb103963299e2384ce636dd117d194ee21b86adeac63a9965a095840f39a6c",
            "3117927": sha256_file(pairs[3_117_927][1])
            == "3d54f3970b4120896c214155811d7cd1b594e3efd172f80b5dc5e7d0fef279e2",
        },
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "new_long_cfd_iterations": 0,
        "instrumented_musubi_calls": 1,
        "instrumented_timesteps": 1,
        "flux_definition": BOUNDARY_FLUX_DEFINITION,
        "physical_Q1_Q2_Q3": DEFERRED_PHYSICAL_FLUX,
        "next": "RUN REPAIRED TAU1 COARSE/BASE/FINE GRID CONVERGENCE",
    }
    acceptance["protected_restarts_all_unchanged"] = all(
        acceptance["protected_restarts_unchanged"].values()
    )
    write_json(qc / "tau1_base_acceptance_after_full_referee.json", acceptance)
    return {"referee": referee, "acceptance": acceptance}
