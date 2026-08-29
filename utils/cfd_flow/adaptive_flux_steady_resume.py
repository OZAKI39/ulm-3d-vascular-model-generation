"""One isolated 0.5% adaptive-flux continuation from iteration 154326."""

from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any

from .adaptive_flux_steady import (
    AXIS_MESH_RUN,
    BULK_NU_M2_S,
    CHECKPOINT_INTERVAL,
    CONVERGENCE_INTERVAL,
    EXPECTED_FLUID_CELLS,
    EXPECTED_INLET_GLOBBC,
    EXPECTED_RESTART_BYTES,
    MPI_RANKS,
    NEXT_REVIEW_WSL,
    NU_M2_S,
    PRESSURE_REFERENCE_PA,
    RUNTIME_BASE_WSL,
    SHORT_SIMULATION_NAME,
    STEADY_FAILED,
    STEADY_PENDING_AUDIT,
    WRAPPER_MARGIN_S,
    WSL_RUNTIME_UNSTABLE,
    CheckpointTracker,
    _archive_runtime,
    _git,
    _maximum_runtime_output_path_length,
    _one_wsl_health_check,
    _predicted_maximum_runtime_path_length,
    _prepare_short_tracking_root,
    _run_luac,
    _run_monitored_wsl,
    _wsl_path_to_unc,
    summarize_controller_csv,
)
from .adaptive_flux_validation import (
    BINARY_WSL,
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


RESUME_SOURCE_RUN = "mcclure_adaptive_flux_steady_anchor003274_20260829_153655"
RESUME_ITERATION = 154_326
EXPECTED_FIRST_CONTINUATION_ITERATION = RESUME_ITERATION + 1
RESUME_BINARY_SHA256 = "0edfe80bdf05f5173f185a261975ac43cdeb70719ac775074be8177ce7f7531e"
PRESSURE_THRESHOLD_0P5_PA = 0.729525879482435
VELOCITY_THRESHOLD_0P5_M_S = 4.919279003768086e-7
CONVERGENCE_NVALS_0P5 = 50
EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS = 5_000
RESUME_PREFIX = "mcclure_adaptive_flux_resume_0p5_anchor003274"

CONTINUITY_FAILED = "CFD_FLOW_ADAPTIVE_FLUX_RESTART_CONTINUITY_FAILED"
STEADY_INCOMPLETE = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_STEADY_INCOMPLETE_0P5_PERCENT"
STEADY_PASS_0P5 = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_STEADY_BASELINE_PASS_0P5_PERCENT"
NEXT_INCOMPLETE = "REVIEW 1_PERCENT STEADY ACCEPTANCE OR CONTINUE CURRENT CRITERION"
NEXT_CONTINUITY = "REVIEW RESTART CONFIGURATION BEFORE RETRY"
NEXT_EXACT_AUDIT = "RUN INDEPENDENT FINAL-RESTART PDF-LINK AUDIT"
NEXT_GRID = "RUN ADAPTIVE-FLUX GRID CONVERGENCE"


def generate_adaptive_resume_lua(
    *,
    mesh_wsl: str,
    resume_header_wsl: str,
    maximum_iterations: int,
    wallclock_limit_s: int,
) -> str:
    """Render the accepted physics with only the requested termination change."""

    checkpoint_start = RESUME_ITERATION + CHECKPOINT_INTERVAL
    return f"""-- Isolated adaptive-flux continuation; restart from iteration {RESUME_ITERATION}.
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
      time_control = {{ min = {{ iter = {RESUME_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
      norm = 'average',
      nvals = {CONVERGENCE_NVALS_0P5},
      absolute = true,
      condition = {{
        {{ threshold = {PRESSURE_THRESHOLD_0P5_PA}, operator = '<=' }},
        {{ threshold = {VELOCITY_THRESHOLD_0P5_M_S}, operator = '<=' }}
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
    label = 'p', folder = 'tracking/p/',
    variable = {{ 'pressure_phy' }}, shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {RESUME_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }},
  {{
    label = 'u', folder = 'tracking/u/',
    variable = {{ 'vel_mag_phy' }}, shape = {{ kind = 'all' }},
    reduction = {{ 'average' }},
    time_control = {{ min = {{ iter = {RESUME_ITERATION} }}, max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }},
    output = {{ format = 'ascii' }}
  }}
}}

restart = {{
  read = '{resume_header_wsl}',
  write = 'restart/',
  time_control = {{
    min = {{ iter = {checkpoint_start} }},
    max = {{ iter = maximum_iterations }},
    interval = {{ iter = {CHECKPOINT_INTERVAL} }}
  }}
}}
"""


def adaptive_resume_lua_contract(text: str, *, resume_header_wsl: str) -> dict[str, Any]:
    convergence_time_control = (
        f"time_control = {{ min = {{ iter = {RESUME_ITERATION} }}, "
        f"max = {{ iter = maximum_iterations }}, interval = {{ iter = {CONVERGENCE_INTERVAL} }} }}"
    )
    checks = {
        "restart_read": f"read = '{resume_header_wsl}'" in text,
        "not_fresh_initialization": "read =" in text,
        "adaptive_inlet": "kind = 'adaptive_flux_pressure'" in text,
        "target_mass_flow_frozen": f"mass_flowrate = {TARGET_MASS_FLOW_KG_S:.17g}" in text,
        "pressure_outlets_frozen": text.count("kind = 'pressure_eq'") == 3,
        "wall_libb_frozen": "kind = 'wall_libb'" in text,
        "d3q19_bgk_frozen": "layout = 'd3q19'" in text and "relaxation = 'bgk'" in text,
        "dx_frozen": f"dx = {EXPECTED_DX_M:.17g}" in text,
        "dt_frozen": f"dt = {EXPECTED_DT_S:.17g}" in text,
        "rho_frozen": f"rho0_phy = {REFERENCE_DENSITY_KG_M3:.17g}" in text,
        "nu_frozen": f"nu_phy = {NU_M2_S:.17g}" in text,
        "bulk_nu_frozen": f"bulk_viscosity_phy = {BULK_NU_M2_S:.17g}" in text,
        "criterion_variables_frozen": "variable = { 'pressure_phy', 'vel_mag_phy' }" in text,
        "criterion_reduction_frozen": "reduction = { 'average', 'average' }" in text,
        "criterion_shape_frozen": "shape = { kind = 'all' }" in text,
        "criterion_norm_frozen": "norm = 'average'" in text and "absolute = true" in text,
        "pressure_threshold_0p5": f"threshold = {PRESSURE_THRESHOLD_0P5_PA}" in text,
        "velocity_threshold_0p5": f"threshold = {VELOCITY_THRESHOLD_0P5_M_S}" in text,
        "nvals_50": f"nvals = {CONVERGENCE_NVALS_0P5}" in text,
        "interval_100": f"interval = {{ iter = {CONVERGENCE_INTERVAL} }}" in text,
        "convergence_min_resume": convergence_time_control in text,
        "tracking_min_resume": text.count(convergence_time_control) == 3,
        "checkpoint_interval_20000": f"interval = {{ iter = {CHECKPOINT_INTERVAL} }}" in text,
        "no_harvester": "harvest" not in text.lower(),
        "no_full_field_export": "format = 'vtk'" not in text,
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _archive_and_verify_checkpoints(
    *, tracker: CheckpointTracker, runtime_root_wsl: str, run_root: Path
) -> None:
    for checkpoint in tracker.records:
        header_relative = PurePosixPath(checkpoint["header_path_wsl"]).relative_to(
            PurePosixPath(runtime_root_wsl)
        )
        binary_relative = PurePosixPath(checkpoint["binary_path_wsl"]).relative_to(
            PurePosixPath(runtime_root_wsl)
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
    tracker._write()


def run_adaptive_flux_steady_resume(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    head = _git(root, "rev-parse", "HEAD")
    branch = _git(root, "branch", "--show-current")
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root
    source_root = output_root / RESUME_SOURCE_RUN
    resume_header = source_root / "restart" / f"{SHORT_SIMULATION_NAME}_lastHeader.lua"
    resume_binary = source_root / "restart" / f"{SHORT_SIMULATION_NAME}_3.768E-03.lsb"
    mesh = output_root / AXIS_MESH_RUN / "seeder" / "mesh"
    inlet_rim = output_root / AXIS_MESH_RUN / "qc" / "inlet_rim_audit.json"
    binary_windows = Path(
        r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees\musubi_mcclure_adaptive_flux_20260829_1300\build\musubi_adaptive_flux"
    )

    binary_sha = sha256_file(binary_windows)
    resume_sha = sha256_file(resume_binary)
    resume = parse_restart_header(resume_header)
    resume_contract = restart_binary_size_contract(
        resume_binary,
        n_elems=resume.n_elems,
        n_components=resume.n_components,
        n_dofs=resume.n_dofs,
    )
    if binary_sha != EXPECTED_BINARY_SHA256:
        raise FlowError(STEADY_FAILED, f"Adaptive binary SHA changed: {binary_sha}")
    if resume_sha != RESUME_BINARY_SHA256:
        raise FlowError(STEADY_FAILED, f"Resume binary SHA changed: {resume_sha}")
    if (
        resume.iteration != RESUME_ITERATION
        or resume.n_elems != EXPECTED_FLUID_CELLS
        or resume.n_components != 19
        or resume.n_dofs != 1
        or resume_contract["status"] != "PASS"
    ):
        raise FlowError(STEADY_FAILED, f"Resume checkpoint contract failed: {resume}")
    if parse_mesh_header(mesh)["fluid_element_count"] != EXPECTED_FLUID_CELLS:
        raise FlowError(STEADY_FAILED, "Frozen mesh fluid count changed")
    if config.mesh.dx_target_m != EXPECTED_DX_M:
        raise FlowError(STEADY_FAILED, "dx changed")
    if config.solver.maximum_iterations != 1_000_000:
        raise FlowError(STEADY_FAILED, "maximum_iterations changed")
    if config.solver.wallclock_limit_s != 3_600:
        raise FlowError(STEADY_FAILED, "wallclock_limit_s changed")
    production_diff = _git(
        root,
        "diff",
        "--",
        "cfd_flow.py",
        "configs/cfd_flow.yaml",
        "utils/cfd_flow/pipeline.py",
    )
    if production_diff:
        raise FlowError(STEADY_FAILED, "Production pipeline has local modifications")

    frozen_paths = (
        root / "cfd_flow.py",
        root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
        *sorted(path for path in mesh.iterdir() if path.is_file()),
        inlet_rim,
        resume_header,
        resume_binary,
        binary_windows,
    )
    frozen_before = _file_manifest(frozen_paths)
    local_frozen_before = _file_manifest(frozen_paths[:-1])
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{RESUME_PREFIX}_{stamp}"
    runtime_root_wsl = f"{RUNTIME_BASE_WSL}/c_{stamp[-6:]}"
    runtime_root = _wsl_path_to_unc(config.apes.wsl_distribution, runtime_root_wsl)
    qc_dir = run_root / "qc"
    tracking_dir = run_root / "tracking"
    restart_dir = run_root / "restart"
    qc_dir.mkdir(parents=True, exist_ok=False)
    tracking_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = qc_dir / "adaptive_flux_steady_manifest.json"
    checkpoint_tracker: CheckpointTracker | None = None

    summary: dict[str, Any] = {
        "status": STEADY_FAILED,
        "next": NEXT_INCOMPLETE,
        "run_root": str(run_root),
        "actual_head": head,
        "branch": branch,
        "production_pipeline_modified": False,
        "resume_source_iteration": RESUME_ITERATION,
        "resume_header": str(resume_header),
        "resume_binary": str(resume_binary),
        "resume_checkpoint_sha256": resume_sha,
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
        "fresh_initialization": False,
        "steady_criterion_old": "characteristic-scale 0.1%, nvals=100, interval=100",
        "steady_criterion_new": "characteristic-scale 0.5%, nvals=50, interval=100",
        "pressure_threshold_pa": PRESSURE_THRESHOLD_0P5_PA,
        "velocity_threshold_m_s": VELOCITY_THRESHOLD_0P5_M_S,
        "convergence_nvals": CONVERGENCE_NVALS_0P5,
        "convergence_interval_iterations": CONVERGENCE_INTERVAL,
        "effective_convergence_window_iterations": EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS,
        "maximum_iterations": config.solver.maximum_iterations,
        "wallclock_limit_s": config.solver.wallclock_limit_s,
        "runtime_root_wsl": runtime_root_wsl,
        "synthetic_stream_test": "NOT_RUN_AS_REQUESTED",
        "checkpoint_interval_iterations": CHECKPOINT_INTERVAL,
        "frozen_files_before": frozen_before,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)

    try:
        health = _one_wsl_health_check(config.apes.wsl_distribution)
        write_json(qc_dir / "wsl_health_preflight.json", health)
        if health["status"] != "PASS":
            raise FlowError(WSL_RUNTIME_UNSTABLE, "WSL health preflight failed")
        if runtime_root.exists():
            raise FlowError(STEADY_FAILED, f"Runtime root already exists: {runtime_root_wsl}")
        _prepare_short_tracking_root(
            distribution=config.apes.wsl_distribution, root_wsl=runtime_root_wsl
        )
        # Musubi resolves the header's relative binary_name from its working
        # directory.  Stage the accepted pair without changing the frozen
        # source checkpoint so the official restart reader sees an intact pair.
        staged_resume_header = runtime_root / "restart" / resume_header.name
        staged_resume_binary = runtime_root / "restart" / resume_binary.name
        shutil.copy2(resume_header, staged_resume_header)
        shutil.copy2(resume_binary, staged_resume_binary)
        staged_resume_header_sha = sha256_file(staged_resume_header)
        staged_resume_binary_sha = sha256_file(staged_resume_binary)
        if staged_resume_header_sha != sha256_file(resume_header):
            raise FlowError(STEADY_FAILED, "Staged resume header SHA mismatch")
        if staged_resume_binary_sha != resume_sha:
            raise FlowError(STEADY_FAILED, "Staged resume binary SHA mismatch")
        resume_header_wsl = f"restart/{resume_header.name}"
        lua = generate_adaptive_resume_lua(
            mesh_wsl=windows_to_wsl(mesh, config.apes.wsl_distribution),
            resume_header_wsl=resume_header_wsl,
            maximum_iterations=config.solver.maximum_iterations,
            wallclock_limit_s=config.solver.wallclock_limit_s,
        )
        lua_path = runtime_root / "musubi.lua"
        lua_path.write_text(lua, encoding="utf-8")
        contract = adaptive_resume_lua_contract(lua, resume_header_wsl=resume_header_wsl)
        write_json(qc_dir / "adaptive_flux_resume_lua_contract.json", contract)
        if contract["status"] != "PASS":
            raise FlowError(STEADY_FAILED, f"Resume Lua contract failed: {contract}")
        luac_returncode = _run_luac(
            distribution=config.apes.wsl_distribution,
            workdir_wsl=runtime_root_wsl,
            stdout_path=tracking_dir / "luac_stdout.log",
            stderr_path=tracking_dir / "luac_stderr.log",
        )
        preflight_checks = {
            "wsl_health": health["status"] == "PASS",
            "resume_iteration": resume.iteration == RESUME_ITERATION,
            "resume_dimensions": (
                resume.n_elems == EXPECTED_FLUID_CELLS
                and resume.n_components == 19
                and resume.n_dofs == 1
            ),
            "resume_binary_sha": resume_sha == RESUME_BINARY_SHA256,
            "staged_resume_header_sha": staged_resume_header_sha
            == sha256_file(resume_header),
            "staged_resume_binary_sha": staged_resume_binary_sha == resume_sha,
            "adaptive_binary_sha": binary_sha == EXPECTED_BINARY_SHA256,
            "lua_contract": contract["status"] == "PASS",
            "lua_syntax": luac_returncode == 0,
            "persistent_runtime": runtime_root_wsl.startswith(f"{RUNTIME_BASE_WSL}/"),
            "maximum_path_length_at_most_80": (
                _predicted_maximum_runtime_path_length(runtime_root_wsl) <= 80
            ),
            "production_pipeline_unchanged": not production_diff,
            "synthetic_not_rerun": True,
        }
        preflight = {
            "status": "PASS" if all(preflight_checks.values()) else "FAIL",
            "checks": preflight_checks,
            "runtime_root_wsl": runtime_root_wsl,
            "predicted_maximum_path_length": _predicted_maximum_runtime_path_length(
                runtime_root_wsl
            ),
        }
        write_json(qc_dir / "adaptive_flux_resume_static_preflight.json", preflight)
        if preflight["status"] != "PASS":
            raise FlowError(STEADY_FAILED, f"Resume static preflight failed: {preflight}")

        checkpoint_tracker = CheckpointTracker(
            runtime_root=runtime_root,
            runtime_root_wsl=runtime_root_wsl,
            manifest_path=run_root / "checkpoint_manifest.json",
        )
        summary["musubi_calls"] = 1
        summary["status"] = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_RESUME_0P5_RUNNING"
        write_json(manifest_path, summary)
        run = _run_monitored_wsl(
            distribution=config.apes.wsl_distribution,
            workdir_wsl=runtime_root_wsl,
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
            minimum_controller_iteration_exclusive=RESUME_ITERATION,
        )
        summary.update(
            {
                "musubi_returncode": run.returncode,
                "musubi_wall_time_s": run.wall_time_s,
                "live_safety_failure": run.safety_failure,
                "continuity_failure": run.continuity_failure,
                "wsl_infrastructure_failure": run.infrastructure_failure,
                "wrapper_timeout": run.wrapper_timeout,
                "controller_record_count": run.controller_record_count,
                "last_controller_iteration": run.latest_controller_iteration,
                "checkpoint_iterations": [
                    int(item["iteration"]) for item in checkpoint_tracker.records
                ],
            }
        )
        write_json(manifest_path, summary)
        if run.continuity_failure:
            raise FlowError(CONTINUITY_FAILED, run.continuity_failure)
        if run.infrastructure_failure:
            raise FlowError(WSL_RUNTIME_UNSTABLE, run.infrastructure_failure)
        if run.safety_failure:
            raise FlowError(STEADY_FAILED, run.safety_failure)
        if run.wrapper_timeout:
            raise FlowError(STEADY_FAILED, "Continuation wrapper timeout")
        if run.returncode != 0:
            raise FlowError(STEADY_FAILED, f"Adaptive Musubi returned {run.returncode}")

        maximum_path_length = _maximum_runtime_output_path_length(
            runtime_root_wsl, runtime_root
        )
        if maximum_path_length > 80:
            raise FlowError(STEADY_FAILED, f"Musubi-visible path length {maximum_path_length} > 80")
        _archive_runtime(runtime_root, run_root)
        _archive_and_verify_checkpoints(
            tracker=checkpoint_tracker,
            runtime_root_wsl=runtime_root_wsl,
            run_root=run_root,
        )
        runtime = summarize_controller_csv(tracking_dir / "controller_records.csv")
        if runtime["first_iteration"] != EXPECTED_FIRST_CONTINUATION_ITERATION:
            raise FlowError(
                CONTINUITY_FAILED,
                f"First controller iteration was {runtime['first_iteration']}, expected {EXPECTED_FIRST_CONTINUATION_ITERATION}",
            )
        combined_log = (
            (tracking_dir / "musubi_stdout.log").read_text(encoding="utf-8", errors="replace")
            + "\n"
            + (tracking_dir / "musubi_stderr.log").read_text(encoding="utf-8", errors="replace")
        )
        steady = parse_official_steady_termination(combined_log)
        write_json(qc_dir / "adaptive_flux_steady_termination_qc.json", steady)
        final_iteration = int(runtime["final_iteration"])
        additional_iterations = final_iteration - RESUME_ITERATION
        official_steady = bool(steady["official_steady_termination"])
        if official_steady:
            confirmation = int(steady["confirmation_iteration"])
            if confirmation != final_iteration:
                raise FlowError(
                    STEADY_FAILED,
                    f"Steady confirmation {confirmation} != final controller {final_iteration}",
                )
            if additional_iterations < EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS:
                raise FlowError(
                    STEADY_FAILED,
                    "Official steady termination occurred before the complete 5000-iteration window",
                )

        final_header = restart_dir / f"{SHORT_SIMULATION_NAME}_lastHeader.lua"
        if not final_header.is_file():
            raise FlowError(STEADY_FAILED, "Final restart lastHeader is missing")
        final_restart = parse_restart_header(final_header)
        final_binary = restart_dir / final_restart.binary_path.name
        final_contract = restart_binary_size_contract(
            final_binary,
            n_elems=final_restart.n_elems,
            n_components=final_restart.n_components,
            n_dofs=final_restart.n_dofs,
        )
        if (
            final_restart.iteration != final_iteration
            or final_restart.n_elems != EXPECTED_FLUID_CELLS
            or final_restart.n_components != 19
            or final_restart.n_dofs != 1
            or final_contract["status"] != "PASS"
        ):
            raise FlowError(STEADY_FAILED, f"Final restart contract failed: {final_restart}")
        write_json(
            qc_dir / "adaptive_flux_final_restart_qc.json",
            {
                "status": "PASS",
                "iteration": final_restart.iteration,
                "header": str(final_header),
                "header_sha256": sha256_file(final_header),
                "binary": str(final_binary),
                "binary_sha256": sha256_file(final_binary),
                "binary_size_contract": final_contract,
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
            raise FlowError(STEADY_FAILED, "Frozen source changed during continuation")

        status = STEADY_PENDING_AUDIT if official_steady else STEADY_INCOMPLETE
        next_step = NEXT_EXACT_AUDIT if official_steady else NEXT_INCOMPLETE
        summary.update(
            {
                "status": status,
                "next": next_step,
                "steady_criterion": "CHARACTERISTIC_SCALE_0P5_PERCENT_50_VALUE",
                "steady_criterion_status": "PASS" if official_steady else "INCOMPLETE",
                "official_steady_termination": official_steady,
                "steady_confirmation_iteration": (
                    int(steady["confirmation_iteration"]) if official_steady else None
                ),
                "first_continuation_controller_iteration": runtime["first_iteration"],
                "additional_iterations": additional_iterations,
                "final_absolute_iteration": final_iteration,
                "total_steady_iterations": final_iteration,
                "runtime_controller": runtime,
                "maximum_runtime_output_path_length": maximum_path_length,
                "final_restart_header": str(final_header),
                "final_restart_binary": str(final_binary),
                "runtime_root_preserved_through_audit": True,
                "source_frozen_files_unchanged": True,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        status = error.status if isinstance(error, FlowError) else STEADY_FAILED
        if status == CONTINUITY_FAILED:
            next_step = NEXT_CONTINUITY
        elif status == WSL_RUNTIME_UNSTABLE:
            next_step = NEXT_REVIEW_WSL
        elif status == STEADY_INCOMPLETE:
            next_step = NEXT_INCOMPLETE
        else:
            next_step = NEXT_INCOMPLETE
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
        full_unchanged = False
        try:
            full_unchanged = frozen_before == _file_manifest(frozen_paths)
        except (OSError, ValueError):
            pass
        summary.update(
            {
                "status": status,
                "next": next_step,
                "first_failure": str(error),
                "source_frozen_files_unchanged": (
                    full_unchanged or local_frozen_before == local_frozen_after
                ),
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary


def finalize_resume_with_exact_audit(
    *, steady: dict[str, Any], exact_audit: dict[str, Any] | None
) -> dict[str, Any]:
    """Apply the requested public 0.5% status without changing audit mathematics."""

    if exact_audit is None:
        return {
            "status": steady["status"],
            "next": steady["next"],
            "steady": steady,
            "exact_audit": "NOT_RUN",
        }
    from .adaptive_flux_steady_exact_audit import STEADY_BASELINE_PASS

    passed = exact_audit.get("status") == STEADY_BASELINE_PASS
    return {
        "status": STEADY_PASS_0P5 if passed else exact_audit.get("status"),
        "next": NEXT_GRID if passed else exact_audit.get("next"),
        "steady": steady,
        "exact_audit": exact_audit,
    }
