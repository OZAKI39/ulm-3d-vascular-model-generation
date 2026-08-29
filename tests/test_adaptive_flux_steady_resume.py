from __future__ import annotations

from utils.cfd_flow.adaptive_flux_steady import STEADY_PENDING_AUDIT
from utils.cfd_flow.adaptive_flux_steady_exact_audit import STEADY_BASELINE_PASS
from utils.cfd_flow.adaptive_flux_steady_resume import (
    CONVERGENCE_NVALS_0P5,
    EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS,
    EXPECTED_FIRST_CONTINUATION_ITERATION,
    PRESSURE_THRESHOLD_0P5_PA,
    RESUME_ITERATION,
    STEADY_PASS_0P5,
    VELOCITY_THRESHOLD_0P5_M_S,
    adaptive_resume_lua_contract,
    finalize_resume_with_exact_audit,
    generate_adaptive_resume_lua,
)


def test_resume_lua_changes_only_requested_termination_contract() -> None:
    header = "/mnt/e/project/restart/a3274_lastHeader.lua"
    lua = generate_adaptive_resume_lua(
        mesh_wsl="/mnt/e/project/mesh",
        resume_header_wsl=header,
        maximum_iterations=1_000_000,
        wallclock_limit_s=3_600,
    )
    contract = adaptive_resume_lua_contract(lua, resume_header_wsl=header)
    assert contract["status"] == "PASS"
    assert f"read = '{header}'" in lua
    assert f"min = {{ iter = {RESUME_ITERATION} }}" in lua
    assert f"threshold = {PRESSURE_THRESHOLD_0P5_PA}" in lua
    assert f"threshold = {VELOCITY_THRESHOLD_0P5_M_S}" in lua
    assert f"nvals = {CONVERGENCE_NVALS_0P5}" in lua
    assert "kind = 'adaptive_flux_pressure'" in lua
    assert "kind = 'mfr_eq'" not in lua
    assert "clock = 3600" in lua
    assert EFFECTIVE_CONVERGENCE_WINDOW_ITERATIONS == 5_000
    assert EXPECTED_FIRST_CONTINUATION_ITERATION == 154_327


def test_resume_final_status_requires_steady_and_exact_audit_pass() -> None:
    steady = {"status": STEADY_PENDING_AUDIT, "next": "audit"}
    result = finalize_resume_with_exact_audit(
        steady=steady,
        exact_audit={"status": STEADY_BASELINE_PASS, "next": "grid"},
    )
    assert result["status"] == STEADY_PASS_0P5
    assert result["next"] == "RUN ADAPTIVE-FLUX GRID CONVERGENCE"


def test_resume_incomplete_never_runs_exact_audit() -> None:
    steady = {"status": "INCOMPLETE", "next": "review"}
    result = finalize_resume_with_exact_audit(steady=steady, exact_audit=None)
    assert result["status"] == "INCOMPLETE"
    assert result["exact_audit"] == "NOT_RUN"
