from __future__ import annotations

import math

import pytest

from utils.cfd_flow.tau1_grid_convergence import (
    GRID_SPECS,
    PRIMARY_FLUX_DEFINITION,
    PRIMARY_METRICS,
    build_primary_analyses,
    evaluate_repaired_grid_gate,
    three_grid_primary_analysis,
)
from utils.cfd_flow.tau1_reference_scaled_grid import (
    INITIAL_SAFETY_ITERATIONS,
    member_lua_contract,
    plateau_failure,
)


def test_all_grid_runtime_values_are_formula_derived() -> None:
    expected_pressure = {
        "coarse": 2004444.2130177515,
        "base": 3387510.7199999993,
        "fine": 5724893.1168,
    }
    expected_target = {
        "coarse": 5.365234139203659e-4,
        "base": 6.974804380964758e-4,
        "fine": 9.067245695254185e-4,
    }
    for label, spec in GRID_SPECS.items():
        assert spec.dt_s == spec.dx_m**2 / (6.0 * 3.27e-6)
        assert spec.tau == pytest.approx(1.0, abs=1.0e-14)
        assert spec.omega == pytest.approx(1.0, abs=1.0e-14)
        assert spec.pressure_reference_pa == pytest.approx(expected_pressure[label])
        assert spec.target_lattice == pytest.approx(expected_target[label])


def test_grid_insensitive_floor_passes_without_fabricated_order() -> None:
    result = three_grid_primary_analysis(1.005, 1.0, 1.004)
    assert result["classification"] == "GRID_INSENSITIVE_WITHIN_1PCT"
    assert result["observed_order_p"] is None
    assert result["GCI_BF"] is None
    assert result["pass"] is True


def test_oscillatory_primary_is_rejected_even_below_five_percent() -> None:
    result = three_grid_primary_analysis(1.03, 1.0, 1.02)
    assert result["relative_B_F"] < 0.05
    assert result["classification"] == "OSCILLATORY_OR_NONCONVERGENT"
    assert result["pass"] is False


def test_gci_gate_is_independent_of_five_percent_difference_gate() -> None:
    values = {}
    for grid, value in (("coarse", 1.25), ("base", 1.1), ("fine", 1.05)):
        values[grid] = {
            "flux_definition": PRIMARY_FLUX_DEFINITION,
            **{name: value for name in PRIMARY_METRICS},
        }
    gate = evaluate_repaired_grid_gate(build_primary_analyses(values))
    assert gate["status"] == "PASS"
    assert gate["base_fine_primary_within_5_percent"] is True


def test_plateau_requires_four_audits_after_physical_time_floor() -> None:
    audits = [
        {
            "iteration": index,
            "physical_time_s": 0.0015 + index * 1.0e-5,
            "failed_gates": ["R_mass_short"],
            "R_mass_short": 0.011 * (1.0 + index * 0.001),
        }
        for index in range(4)
    ]
    result = plateau_failure(audits)
    assert result is not None
    assert result["failure_mode"] == "SCIENTIFIC_PLATEAU_FAILURE"


def test_member_lua_contract_rejects_cross_grid_restart(tmp_path) -> None:
    # Contract checking is intentionally textual and does not touch a solver.
    spec = GRID_SPECS["coarse"]
    assert INITIAL_SAFETY_ITERATIONS < spec.short_window_iterations
    text = "read='/fine/restart.lua'"
    result = member_lua_contract(
        text,
        tmp_path,
        "coarse",
        maximum_iteration=INITIAL_SAFETY_ITERATIONS,
        restart_header_wsl="/coarse/restart.lua",
    )
    assert result["status"] == "FAIL"
    assert result["checks"]["fresh_or_own_restart"] is False


def test_physical_time_window_rounding_is_within_half_timestep() -> None:
    for spec in GRID_SPECS.values():
        short_error = abs(spec.short_window_iterations * spec.dt_s - 0.0002441406727828746)
        long_error = abs(spec.long_window_iterations * spec.dt_s - 0.0004882813455657492)
        assert short_error <= 0.5 * spec.dt_s
        assert long_error <= 0.5 * spec.dt_s
        assert math.isfinite(spec.pressure_reference_pa)
