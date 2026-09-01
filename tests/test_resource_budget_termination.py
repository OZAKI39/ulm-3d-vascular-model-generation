from __future__ import annotations

import json
from pathlib import Path

import pytest

from utils.cfd_flow.resource_budget_termination import (
    UNAVAILABLE,
    build_resolution_sensitivity,
    relative_difference,
)


def test_relative_difference_uses_base_as_denominator() -> None:
    result = relative_difference(102.0, 100.0)
    assert result["signed_base_minus_coarse_over_abs_base"] == pytest.approx(-0.02)
    assert result["absolute_relative_difference"] == pytest.approx(0.02)
    assert result["absolute_percent_difference"] == pytest.approx(2.0)


def test_two_grid_report_never_claims_formal_convergence() -> None:
    coarse = {
        "inlet_gauge_pressure_pa": 101.0,
        "DeltaP01_pa": 101.0,
        "DeltaP02_pa": 101.0,
        "DeltaP03_pa": 101.0,
        "Qin_m3_s": 1.01,
        "Q1_m3_s": 1.01,
        "Q2_m3_s": 1.01,
        "Q3_m3_s": 1.01,
        "outlet_01_flow_fraction": 0.2,
        "outlet_02_flow_fraction": 0.6,
        "outlet_03_flow_fraction": 0.2,
    }
    base = coarse | {
        "inlet_gauge_pressure_pa": 100.0,
        "DeltaP01_pa": 100.0,
        "DeltaP02_pa": 100.0,
        "DeltaP03_pa": 100.0,
        "Qin_m3_s": 1.0,
        "Q1_m3_s": 1.0,
        "Q2_m3_s": 1.0,
        "Q3_m3_s": 1.0,
    }
    result = build_resolution_sensitivity(coarse, base)
    assert result["status"] == "PASS_TWO_GRID_RESOLUTION_SENSITIVITY"
    assert result["formal_asymptotic_grid_convergence"] is False
    assert result["three_grid_metrics"] == UNAVAILABLE
    assert "grid independent proven" in result["claims_not_made"]


def test_final_resource_budget_qc_contract() -> None:
    root = Path(__file__).parents[1]
    qc = (
        root
        / "outputs/cfd_flow"
        / "healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_anchor003274_20260901"
        / "qc"
    )
    resource = json.loads((qc / "resource_budget_termination.json").read_text(encoding="utf-8"))
    grid = json.loads((qc / "grid_convergence_final.json").read_text(encoding="utf-8"))
    sensitivity = json.loads(
        (qc / "coarse_base_resolution_sensitivity.json").read_text(encoding="utf-8")
    )
    assert resource["classification"] == "FINE_LONG_RUN_TERMINATED_BY_RESOURCE_BUDGET"
    assert resource["scientific_failure"] is False
    assert resource["restart_integrity"]["status"] == "PASS"
    assert grid["coarse"]["accepted"] is True
    assert grid["base"]["accepted"] is True
    assert grid["fine"]["steady_accepted"] is False
    assert grid["fine"]["scientific_failure"] is False
    assert grid["observed_order_p"] == UNAVAILABLE
    assert grid["richardson_extrapolation"] == UNAVAILABLE
    assert grid["gci"] == UNAVAILABLE
    assert sensitivity["formal_asymptotic_grid_convergence"] is False
