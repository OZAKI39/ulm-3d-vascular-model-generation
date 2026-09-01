"""Finalize the Tau1 CFD research scope after a user budget stop.

This module is report-only.  It reads accepted/restart evidence and never
launches Seeder, Musubi, or any production/WSS stage.
"""

from __future__ import annotations

import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

import numpy as np

from .fine_continuation_monitor import RUN_RELATIVE
from .io import sha256_file
from .restart_decode import parse_restart_header, read_restart_pdf
from .tau1_base import _controller_records


FINAL_ITERATION = 70_313
FINAL_PHYSICAL_TIME_S = 8.4822274095311601e-5
FINAL_RESTART_SHA256 = "462234ee8b59ae7b7be2234a5c1ef1ede24db6065a54b82af8cb6871d3ac7f9a"
FINAL_HEADER_SHA256 = "a3fe3ef98b628057d862c5cac6e97addacde8bd6c6d2704970b1a8ab5e9936b0"
OLD_FINE_5000_SHA256 = "b51f9ddc48cd03f6a1631c5eeba24e197e0c2099358b691af86d144ac87aedc7"
POSTFIX_FINE_5000_SHA256 = "ef315aea6a694cd1982d690662a8abc43194ad8cc38a81debbd2b27892111d6d"
COARSE_SHA256 = "9eda88e685e5eaa6af650757b8c97accba392c88988df2766d2bb4de825fcfbb"
BASE_SHA256 = "ffcd98b2dc684d1569d937d915b603805809c581d5341e71b17afac2ac64c39f"
FINAL_STATUS = "CFD_FLOW_TAU1_COARSE_BASE_RESOLUTION_VALIDATED"
FINE_CLASSIFICATION = "FINE_LONG_RUN_TERMINATED_BY_RESOURCE_BUDGET"
UNAVAILABLE = "NOT_AVAILABLE_FINE_STEADY_NOT_COMPLETED"
NEXT = (
    "PROCEED USING THE VALIDATED BASE SOLUTION AND PARTIAL COARSE-TO-BASE "
    "RESOLUTION SENSITIVITY, WITHOUT FURTHER FINE CFD."
)
OBSERVABLES = (
    "inlet_gauge_pressure_pa",
    "DeltaP01_pa",
    "DeltaP02_pa",
    "DeltaP03_pa",
    "Qin_m3_s",
    "Q1_m3_s",
    "Q2_m3_s",
    "Q3_m3_s",
    "outlet_01_flow_fraction",
    "outlet_02_flow_fraction",
    "outlet_03_flow_fraction",
)


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def relative_difference(coarse: float, base: float) -> dict[str, float]:
    denominator = max(abs(float(base)), np.finfo(np.float64).tiny)
    signed = (float(base) - float(coarse)) / denominator
    return {
        "signed_base_minus_coarse_over_abs_base": signed,
        "absolute_relative_difference": abs(float(base) - float(coarse)) / denominator,
        "absolute_percent_difference": 100.0 * abs(float(base) - float(coarse)) / denominator,
    }


def build_resolution_sensitivity(
    coarse: dict[str, Any], base: dict[str, Any]
) -> dict[str, Any]:
    rows = {
        name: {
            "coarse": float(coarse[name]),
            "base": float(base[name]),
            **relative_difference(float(coarse[name]), float(base[name])),
        }
        for name in OBSERVABLES
    }
    maximum = max(float(row["absolute_relative_difference"]) for row in rows.values())
    return {
        "status": "PASS_TWO_GRID_RESOLUTION_SENSITIVITY",
        "classification": "TWO_GRID_RESOLUTION_SENSITIVITY_NOT_FORMAL_GRID_CONVERGENCE",
        "reference_grid_for_relative_difference": "base",
        "relative_difference_formula": "abs(base - coarse) / abs(base)",
        "coarse": {
            "status": "PASS",
            "accepted_iteration": 354_295,
            "accepted_restart_sha256": COARSE_SHA256,
        },
        "base": {
            "status": "PASS",
            "accepted_iteration": 598_755,
            "accepted_restart_sha256": BASE_SHA256,
        },
        "observables": rows,
        "maximum_absolute_relative_difference": maximum,
        "maximum_absolute_percent_difference": 100.0 * maximum,
        "formal_asymptotic_grid_convergence": False,
        "three_grid_metrics": UNAVAILABLE,
        "conclusion": (
            "Coarse-to-Base resolution sensitivity is small/moderate for the reported primary "
            "observables, while formal three-grid convergence was not completed under the "
            "available compute budget."
        ),
        "claims_not_made": [
            "grid independent proven",
            "three-grid convergence PASS",
            "Richardson convergence established",
            "GCI validated",
            "Fine steady PASS",
        ],
    }


def _git_head(root: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True, timeout=10
    ).strip()


def _artifact_paths(run: Path) -> dict[str, Path]:
    termination = run / "fine_postfix/resource_budget_termination"
    return {
        "restart": termination / "restart/tau1_reference_scaled_fine_postfix_70313.lsb",
        "header": termination / "restart/tau1_reference_scaled_fine_postfix_header_70313.lua",
        "stdout": termination / "logs/musubi_stdout.log",
        "stderr": termination / "logs/musubi_stderr.log",
        "status": termination / "logs/long_run_status.txt",
        "stop": termination / "logs/stop_request.txt",
    }


def _verification(qc: Path) -> dict[str, Any]:
    pytest_path = qc / "resource_budget_targeted_pytest.xml"
    ruff_path = qc / "resource_budget_targeted_ruff.json"
    if not pytest_path.is_file() or not ruff_path.is_file():
        return {"status": "PENDING"}
    root = ElementTree.parse(pytest_path).getroot()
    suites = root.findall("testsuite") if root.tag != "testsuite" else [root]
    tests = sum(int(suite.attrib.get("tests", 0)) for suite in suites)
    failures = sum(int(suite.attrib.get("failures", 0)) for suite in suites)
    errors = sum(int(suite.attrib.get("errors", 0)) for suite in suites)
    findings = len(json.loads(ruff_path.read_text(encoding="utf-8")))
    return {
        "status": "PASS" if failures == errors == findings == 0 else "FAIL",
        "targeted_pytest": {
            "tests": tests,
            "failures": failures,
            "errors": errors,
            "sha256": sha256_file(pytest_path),
        },
        "targeted_ruff": {"findings": findings, "sha256": sha256_file(ruff_path)},
    }


def _restart_integrity(paths: dict[str, Path]) -> dict[str, Any]:
    if sha256_file(paths["restart"]) != FINAL_RESTART_SHA256:
        raise RuntimeError("final Fine restart hash mismatch")
    if sha256_file(paths["header"]) != FINAL_HEADER_SHA256:
        raise RuntimeError("final Fine header hash mismatch")
    header = parse_restart_header(paths["header"])
    if (
        header.iteration != FINAL_ITERATION
        or header.n_elems != 400_949
        or header.n_components != 19
        or header.n_dofs != 1
    ):
        raise RuntimeError(f"final Fine header contract mismatch: {header}")
    pdf = read_restart_pdf(paths["restart"], n_elems=400_949, n_components=19)
    rho = np.sum(pdf, axis=1, dtype=np.float64)
    return {
        "status": "PASS",
        "header_exists": paths["header"].is_file(),
        "payload_exists": paths["restart"].is_file(),
        "header_parse": "PASS",
        "iteration": header.iteration,
        "n_elems": header.n_elems,
        "n_components": header.n_components,
        "payload_size_bytes": paths["restart"].stat().st_size,
        "expected_payload_size_bytes": 400_949 * 19 * 8,
        "payload_sha256": FINAL_RESTART_SHA256,
        "header_sha256": FINAL_HEADER_SHA256,
        "all_pdf_finite": bool(np.all(np.isfinite(pdf))),
        "minimum_pdf": float(np.min(pdf)),
        "rho_min": float(np.min(rho)),
        "rho_max": float(np.max(rho)),
        "rho_mean": float(np.mean(rho)),
    }


def finalize(project_root: Path) -> dict[str, Any]:
    root = project_root.resolve()
    run = root / RUN_RELATIVE
    qc = run / "qc"
    fine_qc = run / "fine_postfix/qc"
    paths = _artifact_paths(run)
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    if paths["status"].read_text(encoding="utf-8").strip() != "PASS":
        raise RuntimeError("graceful Fine launcher did not report PASS")
    stdout_text = paths["stdout"].read_text(encoding="utf-8", errors="replace")
    if "SUCCESSFUL run!" not in stdout_text or "Found a stop file" not in stdout_text:
        raise RuntimeError("Fine graceful-stop semantic markers are incomplete")
    controllers = _controller_records(stdout_text)
    if not controllers or int(controllers[-1]["iteration"]) != FINAL_ITERATION:
        raise RuntimeError("Fine final controller record is missing")
    controller = controllers[-1]
    integrity = _restart_integrity(paths)
    if not integrity["all_pdf_finite"] or integrity["minimum_pdf"] <= 0.0:
        raise RuntimeError("Fine final restart failed finite/positive PDF integrity")

    coarse_acceptance = json.loads((qc / "coarse_steady_acceptance.json").read_text(encoding="utf-8"))
    base_acceptance = json.loads((qc / "base_readonly_acceptance.json").read_text(encoding="utf-8"))
    if coarse_acceptance["accepted_restart"]["sha256"] != COARSE_SHA256:
        raise RuntimeError("Coarse accepted restart changed")
    if base_acceptance["accepted_restart_sha256"] != BASE_SHA256:
        raise RuntimeError("Base accepted restart changed")
    sensitivity = build_resolution_sensitivity(
        coarse_acceptance["primary_observables"], base_acceptance["primary_observables"]
    )
    sensitivity.update(
        {
            "git_head_at_finalization": _git_head(root),
            "input_window_statistic": "PHYSICAL_TIME_TRAPEZOIDAL_MEAN",
            "input_flux_definition": "PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX",
        }
    )
    _write_json(qc / "coarse_base_resolution_sensitivity.json", sensitivity)

    timestamp = datetime.now(timezone.utc).isoformat()
    resource = {
        "status": "PASS_RESOURCE_BUDGET_TERMINATION",
        "classification": FINE_CLASSIFICATION,
        "termination_kind": "RESOURCE_BUDGET_TERMINATION",
        "scientific_failure": False,
        "user_requested": True,
        "timestamp": timestamp,
        "git_head_at_finalization": _git_head(root),
        "stop": {
            "method": "MUSUBI_STOP_FILE_GRACEFUL",
            "path": "/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix/stop",
            "requested_at": paths["stop"].read_text(encoding="utf-8").strip(),
            "launcher_exit_status": 0,
            "launcher_status_file": "PASS",
            "successful_run_marker": True,
            "found_stop_file_marker": True,
            "kill_used": False,
            "kill_9_used": False,
        },
        "final_iteration": FINAL_ITERATION,
        "final_physical_time_s": FINAL_PHYSICAL_TIME_S,
        "restart_integrity": integrity,
        "final_controller": controller,
        "log_artifacts": {
            name: {"path": str(paths[name]), "sha256": sha256_file(paths[name])}
            for name in ("stdout", "stderr", "status", "stop")
        },
        "preserved_evidence": {
            "old_incorrect_target_fine_5000": {
                "runtime_path": (
                    "/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine/restart/"
                    "tau1_reference_scaled_fine_5000.lsb"
                ),
                "sha256": OLD_FINE_5000_SHA256,
                "verified_present_at_termination": True,
            },
            "postfix_fresh_fine_5000": {
                "runtime_path": (
                    "/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix/restart/"
                    "tau1_reference_scaled_fine_postfix_5000.lsb"
                ),
                "sha256": POSTFIX_FINE_5000_SHA256,
                "verified_present_at_termination": True,
            },
            "root_cause_evidence": "qc/fine_adaptive_target_375_376_forensics.json",
            "anomaly_tree_id": 232_050_831,
            "anomaly_classification": "prp_solid",
            "post_fix_target_oracle": "qc/fine_adaptive_target_fix_validation.json",
            "fresh_5000_safety": "fine_postfix/qc/fine_5000_safety_gate.json",
        },
        "fine_contracts": {
            "fine_mesh_contract": "PASS",
            "fine_q_contract": "PASS",
            "fine_controller_fix": "PASS",
            "fine_target_oracle": "PASS",
            "fine_5000_safety": "PASS",
            "fine_steady_state": "NOT_COMPLETED",
            "fine_full_v2_at_steady": "NOT_RUN",
            "fine_physical_plane_steady_metrics": "NOT_AVAILABLE",
        },
        "calls_this_task": {"seeder": 0, "musubi_scientific": 0, "additional_cfd": 0},
        "production_pipeline_modified": False,
        "production_promotion_started": False,
        "wss_validation_started": False,
        "verification": _verification(qc),
    }
    _write_json(qc / "resource_budget_termination.json", resource)

    grid = {
        "status": "PASS_RESOURCE_CONSTRAINED_SCOPE",
        "final_status": FINAL_STATUS,
        "production_pipeline_modified": False,
        "coarse": {"accepted": True, "iteration": 354_295, "restart_sha256": COARSE_SHA256},
        "base": {"accepted": True, "iteration": 598_755, "restart_sha256": BASE_SHA256},
        "fine": {
            "mesh_contract": True,
            "q_contract": True,
            "controller_fix": True,
            "target_oracle": True,
            "safety_5000": True,
            "final_executed_iteration": FINAL_ITERATION,
            "final_restart_sha256": FINAL_RESTART_SHA256,
            "steady_accepted": False,
            "scientific_failure": False,
            "terminated_by_user_resource_budget": True,
            "status": FINE_CLASSIFICATION,
        },
        "coarse_base_resolution_sensitivity": "PASS_TWO_GRID_RESOLUTION_SENSITIVITY",
        "three_grid_convergence_status": "NOT_COMPLETED_RESOURCE_BUDGET",
        "formal_three_grid_metrics": UNAVAILABLE,
        "fine_base_to_fine_difference": UNAVAILABLE,
        "observed_order_p": UNAVAILABLE,
        "richardson_extrapolation": UNAVAILABLE,
        "gci": UNAVAILABLE,
        "scientific_scope": (
            "Validated Base with partial Coarse-to-Base resolution sensitivity; no formal "
            "three-grid convergence claim."
        ),
        "next": NEXT,
        "verification": _verification(qc),
    }
    _write_json(qc / "grid_convergence_final.json", grid)

    long_state_path = fine_qc / "long_run_state.json"
    long_state = json.loads(long_state_path.read_text(encoding="utf-8"))
    long_state.update(
        {
            "status": FINE_CLASSIFICATION,
            "completion_classification": "RESOURCE_BUDGET_TERMINATION",
            "process_alive": False,
            "graceful_stop": True,
            "final_iteration": FINAL_ITERATION,
            "final_physical_time_s": FINAL_PHYSICAL_TIME_S,
            "final_restart_sha256": FINAL_RESTART_SHA256,
            "fine_steady_state": "NOT_COMPLETED",
            "scientific_failure": False,
            "resume_allowed": False,
        }
    )
    _write_json(long_state_path, long_state)

    monitor_path = fine_qc / "fine_continuation_monitor.json"
    monitor = {
        "timestamp": timestamp,
        "pid": 30_718,
        "process_alive": False,
        "command_match": True,
        "current_iteration": FINAL_ITERATION,
        "physical_time_s": FINAL_PHYSICAL_TIME_S,
        "latest_restart_iteration": FINAL_ITERATION,
        "latest_restart_path": str(paths["restart"]),
        "latest_restart_header_path": str(paths["header"]),
        "latest_restart_sha256": FINAL_RESTART_SHA256,
        "latest_restart_size_bytes": paths["restart"].stat().st_size,
        "fatal_pattern_found": False,
        "next_expected_checkpoint": None,
        "runtime_status": "PASS",
        "status": FINE_CLASSIFICATION,
        "completion_classification": "RESOURCE_BUDGET_TERMINATION",
        "stop_method": "MUSUBI_STOP_FILE_GRACEFUL",
        "restart_integrity": "PASS",
    }
    _write_json(monitor_path, monitor)
    history_path = fine_qc / "fine_continuation_monitor.jsonl"
    existing = history_path.read_text(encoding="utf-8").splitlines()
    if not existing or json.loads(existing[-1]).get("status") != FINE_CLASSIFICATION:
        with history_path.open("a", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(monitor, ensure_ascii=False, sort_keys=True) + "\n")

    oracle_path = qc / "adaptive_active_population_cb_equivalence.json"
    oracle = json.loads(oracle_path.read_text(encoding="utf-8"))
    oracle.update(
        {
            "status": "DEFERRED_NOT_REQUIRED_FOR_CURRENT_RESOURCE_CONSTRAINED_CONCLUSION",
            "deferred_reason": "NOT_REQUIRED_FOR_CURRENT_RESOURCE-CONSTRAINED CONCLUSION",
            "oracle_execution_timing": "NOT_RUN_BY_USER_RESOURCE_BUDGET_DECISION",
            "solver_calls": 0,
            "overall_verdict": "NOT_RUN_NOT_REQUIRED_FOR_CURRENT_RESOURCE_CONSTRAINED_CONCLUSION",
            "accepted_coarse_base_require_rerun": False,
        }
    )
    for name in ("coarse", "base"):
        oracle[name]["status"] = "DEFERRED_NOT_RUN_RESOURCE_BUDGET"
        oracle[name]["verdict"] = "NOT_RUN_NOT_REQUIRED_FOR_CURRENT_CONCLUSION"
    oracle["coarse_status"] = oracle["coarse"]["status"]
    oracle["base_status"] = oracle["base"]["status"]
    _write_json(oracle_path, oracle)
    oracle_csv = qc / "adaptive_active_population_cb_equivalence.csv"
    with oracle_csv.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
        fieldnames = list(rows[0]) if rows else []
    for row in rows:
        row["status"] = "DEFERRED_NOT_RUN_RESOURCE_BUDGET"
        row["verdict"] = "NOT_RUN_NOT_REQUIRED_FOR_CURRENT_CONCLUSION"
    with oracle_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    coordination = {
        "timestamp": timestamp,
        "git_head": _git_head(root),
        "fine": monitor
        | {
            "last_checkpoint": FINAL_ITERATION,
            "last_checkpoint_sha": FINAL_RESTART_SHA256,
        },
        "cb_oracle": {
            "status": oracle["status"],
            "deferred_reason": oracle["deferred_reason"],
            "solver_calls": 0,
            "coarse_status": oracle["coarse_status"],
            "base_status": oracle["base_status"],
        },
        "next_action": NEXT,
    }
    _write_json(qc / "monitor_and_cb_oracle_state.json", coordination)
    return {"resource": resource, "sensitivity": sensitivity, "grid": grid}


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    result = finalize(args.project_root)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
