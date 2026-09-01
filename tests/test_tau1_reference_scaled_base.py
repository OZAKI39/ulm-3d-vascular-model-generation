from __future__ import annotations

import math

from utils.cfd_flow.tau1_base import (
    OUTLET_GAUGE_PRESSURE_PA,
    PRESSURE_REFERENCE_PA,
    TAU1_DT_S,
    TARGET_Q_M3_S,
    Tau1BaseRuntimeContract,
)
from utils.cfd_flow.tau1_reference_scaled_base import (
    EARLIEST_AUDIT_ITERATION,
    LONG_WINDOW_ITERATIONS,
    SHORT_WINDOW_ITERATIONS,
    acceptance_transition,
    generate_segment_lua,
    outlet_absolute_pressures,
    plateau_failure,
    restart_compatibility,
    segment_lua_contract,
    steady_window_audit,
)


def _sample(iteration: int, *, gauge: float = 10.0) -> dict:
    q = TARGET_Q_M3_S
    ports = {
        "inlet": {
            "Q_velocity_m3_s": q,
            "Q_rho_u_over_rho0_m3_s": q,
        },
        "outlet_01": {
            "Q_velocity_m3_s": 0.2 * q,
            "Q_rho_u_over_rho0_m3_s": 0.2 * q,
        },
        "outlet_02": {
            "Q_velocity_m3_s": 0.3 * q,
            "Q_rho_u_over_rho0_m3_s": 0.3 * q,
        },
        "outlet_03": {
            "Q_velocity_m3_s": 0.5 * q,
            "Q_rho_u_over_rho0_m3_s": 0.5 * q,
        },
    }
    return {
        "iteration": iteration,
        "physical_time_s": iteration * TAU1_DT_S,
        "rho_lattice": {"mean": 1.0, "median": 1.0, "p1": 1.0, "p99": 1.0},
        "mean_speed_m_s": 1.0e-4,
        "inlet_gauge_pressure_pa": gauge,
        "pressure_drops_pa": {
            label: gauge - value for label, value in OUTLET_GAUGE_PRESSURE_PA.items()
        },
        "ports": ports,
        "Qout_sum_m3_s": q,
        "physical_volume_closure": 0.0,
        "flow_fractions": {
            "outlet_01": 0.2,
            "outlet_02": 0.3,
            "outlet_03": 0.5,
        },
        "Q_density_consistency": {
            label: {"residual": 0.0, "pass": True} for label in ports
        },
        "minimum_pdf": 1.0 / 36.0,
        "maximum_lattice_speed": 1.0e-5,
        "all_finite": True,
        "controller": {
            "target_lattice": Tau1BaseRuntimeContract().target_lattice_flux,
            "relative_error": 0.0,
        },
    }


def _triplet(**end_changes: float) -> list[dict]:
    result = [
        _sample(EARLIEST_AUDIT_ITERATION - LONG_WINDOW_ITERATIONS),
        _sample(EARLIEST_AUDIT_ITERATION - SHORT_WINDOW_ITERATIONS),
        _sample(EARLIEST_AUDIT_ITERATION),
    ]
    result[-1].update(end_changes)
    return result


def test_fresh_base_lua_uses_dynamic_reference_and_no_restart() -> None:
    lua = generate_segment_lua(
        maximum_iteration=2 * SHORT_WINDOW_ITERATIONS,
        restart_header_wsl=None,
        restart_first_iteration=SHORT_WINDOW_ITERATIONS,
    )
    assert f"pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}" in lua
    assert "read='" not in lua
    assert "velocityX=0.0" in lua
    assert segment_lua_contract(
        lua,
        maximum_iteration=2 * SHORT_WINDOW_ITERATIONS,
        restart_header_wsl=None,
        restart_first_iteration=SHORT_WINDOW_ITERATIONS,
    )["status"] == "PASS"


def test_old_misscaled_restart_is_forbidden() -> None:
    lua = generate_segment_lua(
        maximum_iteration=2 * SHORT_WINDOW_ITERATIONS,
        restart_header_wsl=None,
        restart_first_iteration=SHORT_WINDOW_ITERATIONS,
    )
    assert all(value not in lua for value in ("2878425", "2998176", "3117927"))


def test_restart_compatibility_includes_dynamic_reference_pressure() -> None:
    expected = {"pressure_reference_pa": PRESSURE_REFERENCE_PA, "tau": 1.0}
    assert restart_compatibility(expected, expected)["status"] == "PASS"
    changed = dict(expected, pressure_reference_pa=23_622.320128)
    result = restart_compatibility(changed, expected)
    assert result["status"] == "FAIL"
    assert not result["checks"]["pressure_reference_pa"]


def test_physical_time_windows_are_exact() -> None:
    assert SHORT_WINDOW_ITERATIONS == 119_751
    assert LONG_WINDOW_ITERATIONS == 2 * SHORT_WINDOW_ITERATIONS
    assert EARLIEST_AUDIT_ITERATION == 4 * SHORT_WINDOW_ITERATIONS
    assert math.isclose(
        SHORT_WINDOW_ITERATIONS * TAU1_DT_S,
        0.0002441406727828746,
        rel_tol=0.0,
        abs_tol=1.0e-20,
    )


def test_physical_volume_closure_gate() -> None:
    audit = steady_window_audit(_triplet(), all_checkpoint_rho_pass=True)
    assert audit["physical_volume_closure"] == 0.0
    assert audit["gates"]["physical_volume_closure"]
    assert audit["R_mass_short"] < 1.0e-15
    assert audit["R_mass_long"] < 1.0e-15


def test_pressure_residual_uses_gauge_not_absolute_reference() -> None:
    samples = _triplet()
    samples[-2]["inlet_gauge_pressure_pa"] = 1.0
    samples[-1]["inlet_gauge_pressure_pa"] = 2.0
    audit = steady_window_audit(samples, all_checkpoint_rho_pass=True)
    assert math.isclose(
        audit["pressure_residuals"]["inlet_gauge_pressure_pa"], 0.5
    )
    assert audit["R_pressure"] >= 0.5


def test_rho_gate_covers_every_checkpoint() -> None:
    audit = steady_window_audit(_triplet(), all_checkpoint_rho_pass=False)
    assert not audit["gates"]["rho_sanity_all_checkpoints"]
    assert "rho_sanity_all_checkpoints" in audit["failed_gates"]


def test_full_referee_outlet_pressures_are_dynamic_reference_plus_gauge() -> None:
    absolute = outlet_absolute_pressures()
    for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items():
        assert math.isclose(
            absolute[label] - PRESSURE_REFERENCE_PA,
            gauge,
            rel_tol=0.0,
            abs_tol=5.0e-10,
        )


def test_early_stop_requires_one_additional_checkpoint() -> None:
    action, candidate = acceptance_transition(
        candidate_iteration=None,
        current_audit_pass=True,
        iteration=EARLIEST_AUDIT_ITERATION,
    )
    assert (action, candidate) == ("CANDIDATE", EARLIEST_AUDIT_ITERATION)
    action, candidate = acceptance_transition(
        candidate_iteration=candidate,
        current_audit_pass=True,
        iteration=EARLIEST_AUDIT_ITERATION + SHORT_WINDOW_ITERATIONS,
    )
    assert action == "CONFIRMED"


def test_plateau_fail_safe_requires_four_same_single_gate_audits() -> None:
    audits = [
        {
            "iteration": 4_430_787 + index * SHORT_WINDOW_ITERATIONS,
            "physical_time_s": (4_430_787 + index * SHORT_WINDOW_ITERATIONS)
            * TAU1_DT_S,
            "failed_gates": ["R_pressure"],
            "R_pressure": 0.006 + index * 1.0e-6,
        }
        for index in range(4)
    ]
    result = plateau_failure(audits)
    assert result is not None
    assert result["failure_mode"] == "SCIENTIFIC_PLATEAU_FAILURE"
    assert plateau_failure(audits[:3]) is None


def test_decision_tree_resets_candidate_after_failed_confirmation() -> None:
    action, candidate = acceptance_transition(
        candidate_iteration=EARLIEST_AUDIT_ITERATION,
        current_audit_pass=False,
        iteration=EARLIEST_AUDIT_ITERATION + SHORT_WINDOW_ITERATIONS,
    )
    assert (action, candidate) == ("CONTINUE", None)
