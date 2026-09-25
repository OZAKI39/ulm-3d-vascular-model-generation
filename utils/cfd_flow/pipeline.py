"""Production Tau1 pipeline and resource-bounded promotion regression."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from .apes import (
    compute_lattice_scaling,
    generate_musubi_lua,
    inspect_apes_environment,
    run_wsl_tool,
    windows_to_wsl,
)
from .config import FlowConfig, METHOD
from .io import FlowError, read_json, sha256_file, write_json
from .production import (
    SUCCESS_STATUS,
    ProductionLayout,
    create_production_layout,
    evaluate_smoke_restart,
    git_head,
    replay_accepted_base,
    validate_full_v2,
    validate_local_artifacts,
    validated_scaling_record,
    write_primary_metrics_csv,
)
from .restart_decode import parse_restart_header
from .visualization import create_production_visuals


SUCCESS_NEXT = "USER VISUALLY REVIEW production_review.html"


@dataclass(frozen=True, slots=True)
class FlowResult:
    status: str
    next: str
    run_root: Path
    summary: dict[str, Any]


def _binary_provenance(config: FlowConfig, environment: Any) -> dict[str, Any]:
    sha_checks: dict[str, Any] = {}
    for record in environment.commands:
        if "binary_sha256_checks" in record:
            sha_checks = record["binary_sha256_checks"]
    result = {
        "status": environment.status,
        "execution_environment": environment.execution_environment,
        "wsl_distribution": environment.wsl_distribution,
        "mpi_ranks": environment.mpi_ranks,
        "seeder": {
            **sha_checks.get("seeder", {}),
            "repository": config.apes.seeder_repository,
            "source_commit": config.apes.seeder_commit,
            "patch_sha256": config.apes.seeder_patch_sha256,
            "geometry_q_contract": "DIMENSIONLESS_GEOMETRY_Q_KERNEL_FROZEN",
        },
        "musubi": {
            **sha_checks.get("musubi", {}),
            "repository": config.apes.musubi_repository,
            "upstream_commit": config.apes.musubi_upstream_commit,
            "scheme_commit": config.apes.musubi_scheme_commit,
            "treelm_commit": config.apes.treelm_commit,
            "patch_sha256": config.apes.musubi_patch_sha256,
            "patched_source_sha256": config.apes.musubi_patched_source_sha256,
            "active_population_denominator": "MPI_SUM(active local boundary counts)",
            "fail_fast_invariants": [
                "pntIndex", "finite sample", "active bitmask", "MPI reduction",
            ],
        },
        "harvester": sha_checks.get("mus_harvesting", {}),
    }
    if result["status"] != "PASS":
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "binary discovery/SHA preflight failed")
    return result


def _render_and_check_lua(
    config: FlowConfig,
    layout: ProductionLayout,
    scaling: Any,
    environment: Any,
) -> tuple[str, dict[str, Any]]:
    mesh_wsl = windows_to_wsl(config.paths.frozen_base_mesh, config.apes.wsl_distribution)
    lua = generate_musubi_lua(
        config,
        None,
        None,
        scaling,
        mesh_path=mesh_wsl.rstrip("/") + "/",
        maximum_iterations=config.execution.solver_smoke_iterations,
    )
    lua_path = layout.solver_smoke / "musubi.lua"
    lua_path.write_text(lua, encoding="utf-8")
    (layout.solver_smoke / "restart").mkdir()
    (layout.solver_smoke / "tracking").mkdir()
    run = run_wsl_tool(
        distribution=config.apes.wsl_distribution,
        workdir=layout.solver_smoke,
        command=[str(environment.binaries["lua_compiler"]), "-p", "musubi.lua"],
        stdout_path=layout.logs / "luac_production_stdout.log",
        stderr_path=layout.logs / "luac_production_stderr.log",
        timeout_s=30,
    )
    checks = {
        "lua_syntax": run.returncode == 0,
        "diffusive_scaling": "scaling = 'diffusive'" in lua,
        "tau1_dt": f"dt = {scaling.dt_s:.17g}" in lua,
        "dynamic_pressure_reference": f"pressure_reference_phy = {scaling.pressure_reference_pa:.17g}" in lua,
        "wall_libb": "kind = 'wall_libb'" in lua,
        "adaptive_flux_pressure": "kind = 'adaptive_flux_pressure'" in lua,
        "mfr_eq_absent": "mfr_eq" not in lua,
        "three_pressure_eq": lua.count("kind = 'pressure_eq'") == 3,
        "fresh_initial_pressure": "initial_condition = { pressure = pressure_reference_phy" in lua,
        "fresh_initial_velocity_zero": all(
            token in lua for token in ("velocityX = 0.0", "velocityY = 0.0", "velocityZ = 0.0")
        ),
        "no_restart_read": "read =" not in lua,
        "iterations_within_budget": config.execution.solver_smoke_iterations <= 5000,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "lua_path": str(lua_path), "lua_sha256": sha256_file(lua_path),
        "checks": checks,
    }
    write_json(layout.input / "production_lua_zero_run_preflight.json", result)
    if result["status"] != "PASS":
        raise FlowError("CFD_FLOW_PRODUCTION_CONFIG_INVALID", f"Lua preflight failed: {checks}")
    return lua, result


def _run_smoke(
    config: FlowConfig,
    layout: ProductionLayout,
    scaling: Any,
    environment: Any,
) -> dict[str, Any]:
    run = run_wsl_tool(
        distribution=config.apes.wsl_distribution,
        workdir=layout.solver_smoke,
        command=[
            str(environment.binaries["mpi_launcher"]),
            "-np", str(environment.mpi_ranks),
            str(environment.binaries["musubi"]), "musubi.lua",
        ],
        stdout_path=layout.logs / "production_smoke_musubi_stdout.log",
        stderr_path=layout.logs / "production_smoke_musubi_stderr.log",
        timeout_s=min(config.resources.wallclock_limit_s, 900),
    )
    if run.returncode != 0:
        raise FlowError(
            "CFD_FLOW_PRODUCTION_PROMOTION_SMOKE_FAILED",
            f"Musubi returned {run.returncode}; see {run.stderr_path}",
        )
    headers = sorted((layout.solver_smoke / "restart").glob("*header_*.lua"))
    exact = [
        path for path in headers
        if parse_restart_header(path).iteration == config.execution.solver_smoke_iterations
    ]
    if len(exact) != 1:
        raise FlowError(
            "CFD_FLOW_PRODUCTION_PROMOTION_SMOKE_FAILED",
            f"expected one terminal restart header, found {len(exact)}",
        )
    stdout = run.stdout_path.read_text(encoding="utf-8", errors="replace")
    result = evaluate_smoke_restart(config, scaling, exact[0], stdout)
    result.update(
        {
            "logical_simulation_calls": 1,
            "process_launches": 1,
            "wall_time_s": run.wall_time_s,
            "musubi_stdout": str(run.stdout_path),
            "musubi_stderr": str(run.stderr_path),
            "production_generated_lua": str(layout.solver_smoke / "musubi.lua"),
            "corrected_musubi_sha256": config.apes.musubi_expected_sha256,
        }
    )
    write_json(layout.qc / "production_smoke_qc.json", result)
    return result


def _adaptive_fix_evidence(project_root: Path, config: FlowConfig) -> dict[str, Any]:
    path = (
        project_root
        / "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901/qc/fine_adaptive_target_fix_validation.json"
    )
    value = read_json(path)
    checks = {
        "status": value.get("status") == "PASS",
        "patch_sha256": value.get("candidate_patch_sha256") == config.apes.musubi_patch_sha256,
        "source_sha256": value.get("candidate_source_sha256") == config.apes.musubi_patched_source_sha256,
        "binary_sha256": value.get("candidate_binary_sha256") == config.apes.musubi_expected_sha256,
        "active_population_count": value.get("active_boundary_count") == 375,
        "initial_population_count": value.get("initial_mesh_boundary_count") == 376,
        "controller_target": value.get("controller_target_pass") is True,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "evidence_path": str(path), "evidence_sha256": sha256_file(path),
        "checks": checks,
    }


def _run_promotion_replay(config: FlowConfig, *, project_root: Path) -> FlowResult:
    layout = create_production_layout(config.paths.output_root)
    shutil.copy2(config.source_path, layout.input / "cfd_flow_promotion_regression.yaml")
    started = datetime.now().isoformat()
    local_artifacts = validate_local_artifacts(config)
    plane_contract = read_json(config.paths.physical_plane_contract)
    inlet_area = float(
        plane_contract["ports"]["inlet"]["planes"]["central"]["aperture_physical_area_m2"]
    )
    scaling = compute_lattice_scaling(config, None, inlet_area)
    scaling_record = validated_scaling_record(config, scaling)
    write_json(layout.input / "production_numerical_contract.json", scaling_record)

    environment = inspect_apes_environment(config.apes)
    binary_provenance = _binary_provenance(config, environment)
    write_json(layout.input / "validated_binary_provenance.json", binary_provenance)
    adaptive_fix = _adaptive_fix_evidence(project_root, config)
    if adaptive_fix["status"] != "PASS":
        raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "adaptive fix evidence mismatch")
    write_json(layout.input / "adaptive_active_population_fix_provenance.json", adaptive_fix)
    _, lua_preflight = _render_and_check_lua(config, layout, scaling, environment)

    if not config.execution.run_solver_smoke:
        raise FlowError(
            "CFD_FLOW_PRODUCTION_PROMOTION_SMOKE_FAILED",
            "promotion regression requires run_solver_smoke=true",
        )
    smoke = _run_smoke(config, layout, scaling, environment)

    replay = replay_accepted_base(config, scaling, layout)
    full_v2 = validate_full_v2(config)
    if full_v2["status"] != "PASS":
        raise FlowError(
            "CFD_FLOW_PRODUCTION_REPLAY_REPRODUCTION_FAILED",
            "accepted Base Full timestep V2 lineage failed",
        )
    write_json(layout.qc / "production_full_timestep_v2.json", full_v2)
    coarse_base = read_json(config.paths.coarse_base_grid_evidence)
    if coarse_base.get("status") != "PASS_TWO_GRID_RESOLUTION_SENSITIVITY":
        raise FlowError(
            "CFD_FLOW_PRODUCTION_REPLAY_REPRODUCTION_FAILED",
            "Coarse-to-Base evidence is not the accepted two-grid result",
        )
    write_json(layout.qc / "two_grid_resolution_sensitivity.json", coarse_base)

    visuals = create_production_visuals(
        grid=replay.grid,
        points_m=replay.points_m,
        metrics=replay.metrics,
        steady_qc=replay.field_qc,
        full_v2=full_v2,
        plane_contract=plane_contract,
        coarse_base=coarse_base,
        output=layout.visualization,
        config=config,
        status="PASS",
    )
    visual_manifest_path = layout.visualization / "visual_manifest.json"
    write_primary_metrics_csv(
        layout.qc / "production_primary_metrics.csv",
        replay.metrics,
        replay.field_qc,
        full_v2,
        config,
    )
    write_json(layout.qc / "production_primary_metrics.json", replay.metrics)

    summary = {
        "status": SUCCESS_STATUS,
        "next": SUCCESS_NEXT,
        "started_at": started,
        "completed_at": datetime.now().isoformat(),
        "method": METHOD,
        "production_schema_version": config.schema_version,
        "git_head": git_head(project_root),
        "execution_mode": config.execution.mode,
        "geometry_source": str(config.paths.source_surface_run),
        "steady_solution_source": "VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART",
        "fresh_full_production_steady_solve": False,
        "fresh_full_steady_production_solve": False,
        "accepted_steady_reference": {
            "iteration": 598_755,
            "physical_time_s": 0.001220703363914373,
            "restart_sha256": "ffcd98b2dc684d1569d937d915b603805809c581d5341e71b17afac2ac64c39f",
        },
        "compute_budget": {
            "seeder_calls": 0,
            "coarse_musubi_calls": 0,
            "fine_musubi_calls": 0,
            "fresh_base_long_steady_calls": 0,
            "fresh_base_smoke_logical_simulations": 1,
            "fresh_base_smoke_iterations": config.execution.solver_smoke_iterations,
            "full_v2_new_solver_calls": 0,
        },
        "numerical_contract": scaling_record,
        "mesh_provenance": local_artifacts,
        "binary_provenance": binary_provenance,
        "adaptive_patch_provenance": adaptive_fix,
        "production_lua_preflight": lua_preflight,
        "fresh_production_smoke": smoke,
        "replayed_primary_metrics": replay.metrics,
        "replay_reproduction": replay.reproduction,
        "steady_qc": replay.field_qc,
        "physical_flux": {
            "definition": replay.flux["flux_definition"],
            "algorithm_revision": replay.flux["algorithm_revision"],
            "plane_contract_revision": replay.flux["plane_contract_revision"],
            "plane_contract_sha256": replay.flux["plane_contract_sha256"],
        },
        "full_timestep_v2": full_v2,
        "three_grid_status": "NOT_COMPLETED_RESOURCE_BUDGET",
        "two_grid_status": "PASS_TWO_GRID_RESOLUTION_SENSITIVITY",
        "fine_status": {
            "mesh": "PASS", "controller_fix": "PASS", "safety_5000": "PASS",
            "steady": "NOT_COMPLETED_RESOURCE_BUDGET",
            "classification": "FINE_LONG_RUN_TERMINATED_BY_RESOURCE_BUDGET",
            "scientific_failure": False,
        },
        "WSS_status": "DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION",
        "new_vtu": {
            "path": str(replay.vtu_path), "sha256": sha256_file(replay.vtu_path),
            "size_bytes": replay.vtu_path.stat().st_size,
        },
        "visualization_paths": [
            str(layout.visualization / item["filename"]) for item in visuals["items"]
        ],
        "production_review_html": str(layout.visualization / "production_review.html"),
        "visualization_manifest": {
            "path": str(visual_manifest_path), "sha256": sha256_file(visual_manifest_path),
            "status": visuals["status"],
        },
    }
    write_json(layout.qc / "run_summary.json", summary)
    return FlowResult(SUCCESS_STATUS, SUCCESS_NEXT, layout.root, summary)


def run_cfd_flow(config: FlowConfig, *, project_root: Path) -> FlowResult:
    """Execute either the promotion replay or an explicitly budgeted fresh route."""

    if config.execution.mode == "VALIDATED_BASE_PROMOTION_REPLAY":
        return _run_promotion_replay(config, project_root=Path(project_root).resolve())
    raise FlowError(
        "CFD_FLOW_FRESH_STEADY_REQUIRES_EXPLICIT_COMPUTE_AUTHORIZATION",
        "FRESH_STEADY is the formal production contract. Use it only in a separately "
        "authorized long-compute session; this command will not silently launch a multi-hour solve.",
    )


def print_result(result: FlowResult) -> None:
    print(f"Production CFD output: {result.run_root}")
    print(f"STATUS: {result.status}")
    print(f"NEXT: {result.next}")
