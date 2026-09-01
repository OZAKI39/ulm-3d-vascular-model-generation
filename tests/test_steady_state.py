from __future__ import annotations

from copy import deepcopy

from utils.cfd_flow.steady_state import acceptance_transition, audit_steady_window
from utils.cfd_flow.validated_contract import OUTLET_GAUGE_PRESSURES_PA, TARGET_VOLUME_FLOW_M3_S, ValidatedTau1Contract


WINDOW = 119_751


def _sample(iteration: int) -> dict:
    contract = ValidatedTau1Contract()
    q = TARGET_VOLUME_FLOW_M3_S
    flows = {"inlet": q, "outlet_01": 0.2 * q, "outlet_02": 0.3 * q, "outlet_03": 0.5 * q}
    return {
        "iteration": iteration,
        "mean_speed_m_s": 1.0e-4,
        "inlet_gauge_pressure_pa": 10.0,
        "pressure_drops_pa": {label: 10.0 - value for label, value in OUTLET_GAUGE_PRESSURES_PA.items()},
        "ports": {label: {"Q_velocity_m3_s": value} for label, value in flows.items()},
        "physical_volume_closure": 0.0,
        "flow_fractions": {"outlet_01": 0.2, "outlet_02": 0.3, "outlet_03": 0.5},
        "Q_density_consistency": {label: {"residual": 0.0, "pass": True} for label in flows},
        "minimum_pdf": 1.0 / 36.0,
        "maximum_lattice_speed": 1.0e-5,
        "all_finite": True,
        "controller": {"target_lattice": contract.target_lattice_flux, "relative_error": 0.0},
    }


def _triplet() -> list[dict]:
    return [_sample(2 * WINDOW), _sample(3 * WINDOW), _sample(4 * WINDOW)]


def test_generic_steady_auditor_accepts_exact_reference_sample() -> None:
    result = audit_steady_window(
        _triplet(), all_checkpoint_rho_pass=True, expected_short_window_iterations=WINDOW
    )
    assert result["status"] == "PASS_NON_REFEREE"
    assert not result["failed_gates"]
    assert result["R_mass_short"] < 1.0e-15
    assert result["R_mass_long"] < 1.0e-15


def test_pressure_residual_is_gauge_based() -> None:
    samples = _triplet()
    samples[-2]["inlet_gauge_pressure_pa"] = 1.0
    samples[-1]["inlet_gauge_pressure_pa"] = 2.0
    result = audit_steady_window(samples, all_checkpoint_rho_pass=True)
    assert result["pressure_residuals"]["inlet_gauge_pressure_pa"] == 0.5
    assert result["status"] == "FAIL"


def test_all_checkpoint_rho_gate_cannot_be_bypassed() -> None:
    result = audit_steady_window(_triplet(), all_checkpoint_rho_pass=False)
    assert "rho_sanity_all_checkpoints" in result["failed_gates"]


def test_averaged_backflow_is_a_hard_gate() -> None:
    samples = deepcopy(_triplet())
    for sample in samples:
        sample["ports"]["outlet_01"]["Q_velocity_m3_s"] = -0.1 * TARGET_VOLUME_FLOW_M3_S
    result = audit_steady_window(samples, all_checkpoint_rho_pass=True)
    assert result["significant_averaged_backflow"] is True
    assert not result["gates"]["no_significant_averaged_backflow"]


def test_acceptance_requires_two_consecutive_passes() -> None:
    action, candidate = acceptance_transition(candidate_iteration=None, current_audit_pass=True, iteration=4 * WINDOW)
    assert (action, candidate) == ("CANDIDATE", 4 * WINDOW)
    action, candidate = acceptance_transition(candidate_iteration=candidate, current_audit_pass=True, iteration=5 * WINDOW)
    assert (action, candidate) == ("CONFIRMED", 4 * WINDOW)
    assert acceptance_transition(candidate_iteration=candidate, current_audit_pass=False, iteration=5 * WINDOW) == ("CONTINUE", None)

