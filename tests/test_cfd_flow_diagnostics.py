"""Targeted tests for the bounded restart diagnostic."""

from __future__ import annotations

import math

import numpy as np

from utils.cfd_flow.diagnostics import (
    END_ITERATION,
    PRESSURE_THRESHOLD_PA,
    START_ITERATION,
    VELOCITY_THRESHOLD_M_S,
    diagnostic_lua_contract,
    generate_diagnostic_musubi_lua,
    parse_total_density_log,
    summarize_convergence,
    summarize_total_density,
)


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
sim_control = { time_control = { max = { iter = maximum_iterations, clock = 3600 } } }
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
        tracking_folder_wsl="/tmp/tracking",
        restart_folder_wsl="/tmp/tracking/restart",
        timing_file_wsl="/tmp/tracking/timing.res",
    )
    contract = diagnostic_lua_contract(production, diagnostic)
    assert contract["status"] == "PASS"
    assert f"maximum_iterations = {END_ITERATION}" in diagnostic
    assert f"iter = {START_ITERATION}" in diagnostic
    assert "format = 'vtk'" not in diagnostic


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

