from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.adaptive_flux_steady import (
    CHECKPOINT_INTERVAL,
    EXPECTED_INLET_GLOBBC,
    SHORT_SIMULATION_NAME,
    generate_adaptive_steady_lua,
    steady_lua_contract,
    summarize_controller_csv,
    summarize_runtime_controller,
)
from utils.cfd_flow.adaptive_flux_steady_exact_audit import (
    MASS_BALANCE_FAILED,
    RUNTIME_TARGET_MISMATCH,
    STEADY_BASELINE_PASS,
    classify_steady_exact_audit,
)


def _controller_line(iteration: int, *, speed: float = 1.0e-5, minimum_pdf: float = 0.02) -> str:
    return (
        f"ADAPTIVE_FLUX_PRESSURE iter={iteration} "
        "target_lattice=2.3478724595760924E-003 "
        "controlled_lattice=2.3478724595760924E-003 "
        "relative_error=1.0000000000000000E-014 "
        "rho_boundary=1.0001000000000000E+000 "
        "pressure_pa=2.3624680000000000E+004 "
        f"max_lattice_velocity={speed:.16E} "
        f"minimum_pdf={minimum_pdf:.16E} "
        f"globBC_count={EXPECTED_INLET_GLOBBC}\n"
    )


def test_steady_lua_preserves_accepted_contract() -> None:
    lua = generate_adaptive_steady_lua(
        mesh_wsl="/mesh",
        restart_wsl="restart",
        pressure_tracking_wsl="tracking/p",
        velocity_tracking_wsl="tracking/u",
        maximum_iterations=1_000_000,
        wallclock_limit_s=3600,
    )
    contract = steady_lua_contract(lua)
    assert contract["status"] == "PASS"
    assert "read =" not in lua
    assert "adaptive_flux_pressure" in lua
    assert "mfr_eq" not in lua
    assert "clock = 3600" in lua
    assert "folder = 'tracking/p/'" in lua
    assert "folder = 'tracking/u/'" in lua
    assert f"simulation_name = '{SHORT_SIMULATION_NAME}'" in lua
    assert "write = 'restart/'" in lua
    assert f"interval = {{ iter = {CHECKPOINT_INTERVAL} }}" in lua
    assert "timing_file = 'timing/timing.res'" in lua


def test_runtime_controller_summary_uses_records(tmp_path: Path) -> None:
    stdout = tmp_path / "stdout.log"
    stderr = tmp_path / "stderr.log"
    stdout.write_text(_controller_line(1) + _controller_line(2), encoding="utf-8")
    stderr.write_text("", encoding="utf-8")
    summary = summarize_runtime_controller((stdout, stderr))
    assert summary["status"] == "PASS"
    assert summary["record_count"] == 2
    assert summary["final_iteration"] == 2
    assert summary["globbc_counts"] == [EXPECTED_INLET_GLOBBC]


def test_persisted_controller_csv_summary(tmp_path: Path) -> None:
    path = tmp_path / "controller_records.csv"
    path.write_text(
        "iteration,target_lattice,controlled_lattice,relative_error,rho_boundary,"
        "pressure_pa,max_lattice_velocity,minimum_pdf,globBC_count\n"
        "1,0.1,0.1,1e-4,1.0,23622.0,1e-5,0.02,287\n"
        "2,0.1,0.1,2e-4,1.1,23623.0,2e-5,0.01,287\n",
        encoding="utf-8",
    )
    summary = summarize_controller_csv(path)
    assert summary["record_count"] == 2
    assert summary["first_iteration"] == 1
    assert summary["final_iteration"] == 2
    assert summary["rho_boundary_range"] == [1.0, 1.1]


def test_exact_audit_classification_order() -> None:
    common = {
        "all_pdfs_finite": True,
        "maximum_lattice_speed": 1.0e-4,
        "minimum_pdf": 0.02,
    }
    assert classify_steady_exact_audit(
        inlet_relative_error=1.0e-4,
        mass_balance_relative_error=2.0e-4,
        **common,
    )[0] == STEADY_BASELINE_PASS
    assert classify_steady_exact_audit(
        inlet_relative_error=0.02,
        mass_balance_relative_error=0.0,
        **common,
    )[0] == RUNTIME_TARGET_MISMATCH
    assert classify_steady_exact_audit(
        inlet_relative_error=0.0,
        mass_balance_relative_error=0.02,
        **common,
    )[0] == MASS_BALANCE_FAILED
