"""Targeted tests for the bounded restart diagnostic."""

from __future__ import annotations

import math
import inspect

import numpy as np
import pytest

from utils.cfd_flow.diagnostics import (
    ADDITIONAL_ITERATIONS,
    ASCII_PRESSURE_FOLDER_WSL,
    ASCII_PRESSURE_LABEL,
    ASCII_VELOCITY_FOLDER_WSL,
    ASCII_VELOCITY_LABEL,
    DIAGNOSTIC_WALLCLOCK_S,
    END_ITERATION,
    OFFICIAL_EQUIVALENT_NVALS,
    PROJECT_STEADY_ADDITIONAL_ITERATIONS,
    PROJECT_STEADY_END_ITERATION,
    PROJECT_STEADY_POLICY,
    PROJECT_STEADY_RELATIVE_FRACTION,
    PRESSURE_THRESHOLD_PA,
    START_ITERATION,
    VELOCITY_THRESHOLD_M_S,
    WRAPPER_HARD_TIMEOUT_S,
    _parse_restart_start_iteration,
    _runtime_numerical_errors,
    ascii_path_preflight,
    build_project_criterion_review,
    classify_official_equivalent,
    classify_convergence_probe,
    combine_continuous_tracking,
    diagnostic_lua_contract,
    generate_diagnostic_musubi_lua,
    official_equivalent_residuals,
    parse_total_density_log,
    parse_official_steady_termination,
    run_project_steady_confirmation,
    summarize_convergence,
    summarize_total_density,
)
from utils.cfd_flow.io import FlowError


def _production_lua() -> str:
    return """-- Generated production configuration; official Musubi syntax.
timing_file = 'mus_timing.res'
mesh = '../../musubi_recovery_anchor003274_20260828_162530/seeder/mesh/'
dx = 1.9999999999999999e-07
dt = 2.4414062499999991e-08
rho0_phy = 1056
nu_phy = 3.27e-06
bulk_viscosity_phy = (2.0 / 3.0) * nu_phy
maximum_iterations = 1000000
function outlet_01_pressure(x, y, z, t) return 23636.865106101286 end
function outlet_02_pressure(x, y, z, t) return 23754.524677223188 end
function outlet_03_pressure(x, y, z, t) return 23608.619501326699 end
sim_control = {
  time_control = { max = { iter = maximum_iterations, clock = 3600 } },
  abort_criteria = { convergence = {
    variable = { 'pressure_phy', 'vel_mag_phy' },
    reduction = { 'average', 'average' },
    time_control = { min = { iter = 0 }, max = { iter = maximum_iterations }, interval = { iter = 100 } },
    norm = 'average', nvals = 100, absolute = true,
    condition = {
      { threshold = 0.001, operator = '<=' },
      { threshold = 1.0000000000000001e-09, operator = '<=' }
    }
  } }
}
identify = { label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }
boundary_condition = {
  { label = 'wall', kind = 'wall_libb' },
  { label = 'inlet', kind = 'mfr_eq', mass_flowrate = 8.124344950169123e-13 },
  { label = 'outlet_01', kind = 'pressure_eq', pressure = outlet_01_pressure },
  { label = 'outlet_02', kind = 'pressure_eq', pressure = outlet_02_pressure },
  { label = 'outlet_03', kind = 'pressure_eq', pressure = outlet_03_pressure }
}
restart = { write = 'restart/' }
"""


def test_parse_total_density_and_exact_windows():
    text = """
  iterations: 0
  simTime   : 0.0
  wallClock : 0.0
      field | total density
          1 | 100.0
  iterations: 100
  simTime   : 1.0e-3
  wallClock : 2.0
      field | total density
          1 | 101.0
  iterations: 200
  simTime   : 2.0e-3
  wallClock : 4.0
      field | total density
          1 | 102.0
"""
    records = parse_total_density_log(text)
    assert [item["iteration"] for item in records] == [0, 100, 200]
    summary = summarize_total_density(records)
    assert summary["final_iteration"] == 200
    assert math.isclose(summary["windows"]["100"]["delta"], 1.0)


def test_diagnostic_lua_only_adds_bounded_restart_tracking():
    production = _production_lua()
    diagnostic = generate_diagnostic_musubi_lua(
        production,
        pressure_folder_wsl=ASCII_PRESSURE_FOLDER_WSL,
        velocity_folder_wsl=ASCII_VELOCITY_FOLDER_WSL,
        restart_folder_wsl="/tmp/restart",
        timing_file_wsl="/tmp/u3d/timing.res",
    )
    contract = diagnostic_lua_contract(production, diagnostic)
    assert contract["status"] == "PASS"
    assert f"maximum_iterations = {END_ITERATION}" in diagnostic
    assert ADDITIONAL_ITERATIONS == 5_000
    assert DIAGNOSTIC_WALLCLOCK_S == 300
    assert WRAPPER_HARD_TIMEOUT_S == 330
    assert f"iter = {START_ITERATION}" in diagnostic
    assert "clock =" not in diagnostic
    assert f"label = '{ASCII_PRESSURE_LABEL}'" in diagnostic
    assert f"folder = '{ASCII_PRESSURE_FOLDER_WSL}/'" in diagnostic
    assert f"label = '{ASCII_VELOCITY_LABEL}'" in diagnostic
    assert f"folder = '{ASCII_VELOCITY_FOLDER_WSL}/'" in diagnostic
    assert "read = 'restart/roi003274_steady_lbm_lastHeader.lua'" in diagnostic
    assert "format = 'vtk'" not in diagnostic


def test_ascii_path_preflight_models_full_ranked_filename():
    preflight = ascii_path_preflight()
    assert preflight["status"] == "PASS"
    assert preflight["pinned_label_len"] == 80
    assert preflight["old_predicted_filename_length"] == 193
    assert preflight["max_ascii_filename_length"] == 44
    assert preflight["max_ascii_filename_length"] <= 60
    assert "/tmp/u3d/p/p_p00000.res" in preflight["predicted_filename_lengths"]
    assert "/tmp/u3d/u/u_p00007.res" in preflight["predicted_filename_lengths"]


def test_convergence_summary_reports_fixed_threshold_ratios():
    records = []
    pressure = 23_622.0
    velocity = 1.0e-4
    for index in range(61):
        if index:
            pressure += 2.0e-3
            velocity += 2.0e-9
        records.append(
            {
                "iteration": START_ITERATION + index * 100,
                "physical_time_s": index,
                "avg_pressure_pa": pressure,
                "avg_velocity_m_s": velocity,
                "delta_pressure_vs_previous": None if index == 0 else 2.0e-3,
                "delta_velocity_vs_previous": None if index == 0 else 2.0e-9,
            }
        )
    summary = summarize_convergence(records)
    assert np.isclose(summary["pressure"]["ratio_to_threshold"], 2.0)
    assert np.isclose(summary["velocity"]["ratio_to_threshold"], 2.0)
    assert summary["pressure"]["threshold"] == PRESSURE_THRESHOLD_PA
    assert summary["velocity"]["threshold"] == VELOCITY_THRESHOLD_M_S


def test_probe_classification_uses_only_tracking_and_runtime_error():
    convergence = {
        "pressure": {"ratio_to_threshold": 12.0, "trend_interpretation": "DECAYING"},
        "velocity": {
            "ratio_to_threshold": 2.0,
            "trend_interpretation": "CLEARLY_DECAYING",
        },
    }
    assert classify_convergence_probe(
        convergence,
        numerical_runtime_error=False,
    ) == ("STILL_CONVERGING", "CONTINUE_FROM_RESTART_LONGER")
    assert classify_convergence_probe(
        convergence,
        numerical_runtime_error=True,
    ) == ("NUMERICAL_PROBLEM", "REVIEW_NUMERICAL_OR_BOUNDARY_CONDITION")


def test_failed_ascii_initialization_is_diagnostic_io_not_numerical():
    text = """INFORMATION ON THE EXECUTABLE
Restarting from point in time:
  iterations: 162464
At line 427 of hvs_ascii_module.f90
Fortran runtime error: End of record
the job has been aborted
"""
    assert _parse_restart_start_iteration(text) == START_ITERATION
    scan = _runtime_numerical_errors(text, returncode=2)
    assert scan["inf_detected"] is False
    assert scan["fortran_io_error_detected"] is True
    assert scan["numerical_runtime_error"] is False
    assert scan["failure_classification"] == "DIAGNOSTIC_IO_ERROR"


def test_nonzero_return_and_wallclock_message_are_not_numerical_instability():
    text = "Reached maximal wall clock running time"
    scan = _runtime_numerical_errors(text, returncode=2)
    assert scan["time_control_termination_detected"] is True
    assert scan["returncode_nonzero"] is True
    assert scan["numerical_runtime_error"] is False
    assert scan["failure_classification"] == "TIME_CONTROL_TERMINATION"


def _tracking_record(iteration: int, pressure: float, velocity: float) -> dict[str, float | int | None]:
    return {
        "iteration": iteration,
        "physical_time_s": iteration * 2.44140625e-8,
        "avg_pressure_pa": pressure,
        "avg_velocity_m_s": velocity,
        "delta_pressure_vs_previous": None,
        "delta_velocity_vs_previous": None,
    }


def test_official_history_window_excludes_current_and_101_samples_yield_one():
    records = [
        _tracking_record(index * 100, float(index), float(index) * 1.0e-9)
        for index in range(101)
    ]
    residuals = official_equivalent_residuals(records)
    assert OFFICIAL_EQUIVALENT_NVALS == 100
    assert len(residuals) == 1
    assert np.isclose(residuals[0]["pressure_history_mean"], 49.5)
    assert np.isclose(residuals[0]["pressure_current"], 100.0)
    assert np.isclose(residuals[0]["pressure_residual_pa"], 50.5)


def test_111_samples_yield_eleven_official_equivalent_residuals():
    records = [
        _tracking_record(index * 100, 20_000.0 + index, 1.0e-6 + index * 1.0e-9)
        for index in range(111)
    ]
    assert len(official_equivalent_residuals(records)) == 11


def test_shared_restart_endpoint_is_deduplicated_with_continuity_qc():
    previous = [_tracking_record(100, 10.0, 1.0e-6), _tracking_record(200, 11.0, 2.0e-6)]
    current = [_tracking_record(200, 11.0, 2.0e-6), _tracking_record(300, 12.0, 3.0e-6)]
    combined, qc = combine_continuous_tracking(previous, current)
    assert [item["iteration"] for item in combined] == [100, 200, 300]
    assert qc["shared_endpoint_deduplicated"] is True
    assert qc["combined_unique_sample_count"] == 3


def test_discontinuous_shared_endpoint_is_rejected():
    previous = [_tracking_record(100, 10.0, 1.0e-6), _tracking_record(200, 11.0, 2.0e-6)]
    current = [_tracking_record(200, 12.0, 3.0e-6), _tracking_record(300, 13.0, 4.0e-6)]
    with pytest.raises(FlowError, match="CFD_FLOW_TRACKING_CONTINUITY_FAILED"):
        combine_continuous_tracking(previous, current)


def test_official_ratio_uses_formal_threshold_and_adjacent_delta_cannot_classify():
    records = [
        _tracking_record(index * 100, index * 1.0e-5, index * 1.0e-11)
        for index in range(111)
    ]
    residuals = official_equivalent_residuals(records)
    assert np.isclose(
        residuals[-1]["pressure_ratio_to_threshold"],
        residuals[-1]["pressure_residual_pa"] / PRESSURE_THRESHOLD_PA,
    )
    summary = {
        "pressure": {"ratio_to_threshold": 0.5, "trend_interpretation": "PLATEAU"},
        "velocity": {"ratio_to_threshold": 0.5, "trend_interpretation": "PLATEAU"},
    }
    assert classify_official_equivalent(
        summary,
        numerical_runtime_error=False,
    ) == ("OFFLINE_EQUIVALENT_STEADY_PASS", "RUN_SINGLE_FORMAL_STEADY_CONFIRMATION")


def test_project_pressure_scale_uses_gauge_span_and_ignores_reference_offset():
    pressures = [14.544978101274268, 132.20454922317552, -13.700626673311461]
    shifted = [value + 23_622.320128 for value in pressures]
    kwargs = {
        "inlet_target_mean_velocity_m_s": 9.838558007536173e-5,
        "current_pressure_residual_pa": 0.15522996164509095,
        "current_velocity_residual_m_s": 1.7014817385342933e-7,
    }
    review = build_project_criterion_review(
        outlet_gauge_pressures_pa=pressures,
        **kwargs,
    )
    shifted_review = build_project_criterion_review(
        outlet_gauge_pressures_pa=shifted,
        **kwargs,
    )
    assert np.isclose(
        review["pressure_characteristic_scale"]["value_pa"],
        145.90517589648698,
    )
    assert np.isclose(
        review["pressure_characteristic_scale"]["value_pa"],
        shifted_review["pressure_characteristic_scale"]["value_pa"],
    )


def test_project_fraction_and_fixed_velocity_scale_convert_exactly_once():
    velocity_scale = 9.838558007536173e-5
    review = build_project_criterion_review(
        outlet_gauge_pressures_pa=[14.544978101274268, 132.20454922317552, -13.700626673311461],
        inlet_target_mean_velocity_m_s=velocity_scale,
        current_pressure_residual_pa=0.15522996164509095,
        current_velocity_residual_m_s=1.7014817385342933e-7,
    )
    assert review["criterion_policy_name"] == PROJECT_STEADY_POLICY
    assert PROJECT_STEADY_RELATIVE_FRACTION == 1.0e-3
    assert review["velocity_characteristic_scale"]["value_m_s"] == velocity_scale
    assert np.isclose(review["derived_thresholds"]["pressure_pa"], 0.14590517589648698)
    assert np.isclose(review["derived_thresholds"]["velocity_m_s"], 9.838558007536173e-8)
    contract = review["musubi_convergence_contract"]
    assert contract["absolute"] is True
    assert contract["norm"] == "average"
    assert contract["nvals"] == 100
    assert contract["internal_clock_ceiling"] is False
    assert contract["threshold_sweep"] is False


def test_project_confirmation_has_one_solver_call_and_fixed_30k_ceiling():
    source = inspect.getsource(run_project_steady_confirmation)
    assert source.count("short_run = run_wsl_tool(") == 1
    assert PROJECT_STEADY_ADDITIONAL_ITERATIONS == 30_000
    assert PROJECT_STEADY_END_ITERATION == 203_464


def test_pinned_steady_log_evidence_parser_is_explicit():
    evidence = parse_official_steady_termination("Reached steady state 190464 T")
    assert evidence["official_steady_termination"] is True
    assert evidence["confirmation_iteration"] == 190_464
