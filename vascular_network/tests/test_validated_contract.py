from __future__ import annotations

import math

import pytest

from utils.cfd_flow.validated_contract import (
    BASE_DX_M,
    HISTORICAL_DT_S,
    HISTORICAL_FIXED_REFERENCE_PRESSURE_PA,
    KINEMATIC_VISCOSITY_M2_S,
    OUTLET_GAUGE_PRESSURES_PA,
    TARGET_VOLUME_FLOW_M3_S,
    ValidatedTau1Contract,
    gauge_pressure_pa,
    pressure_reference_pa,
    relaxation_from_physical,
    target_lattice_flux,
    tau1_time_step_s,
)


def test_tau1_diffusive_scaling_has_one_canonical_formula() -> None:
    dt_s = tau1_time_step_s(BASE_DX_M, KINEMATIC_VISCOSITY_M2_S)
    assert dt_s == pytest.approx(2.038735983690112e-9, rel=2.0e-15)
    nu_lattice, tau, omega = relaxation_from_physical(BASE_DX_M, dt_s)
    assert nu_lattice == pytest.approx(1.0 / 6.0, rel=2.0e-15)
    assert tau == pytest.approx(1.0, rel=2.0e-15)
    assert omega == pytest.approx(1.0, rel=2.0e-15)


def test_dynamic_reference_pressure_and_gauge_conversion() -> None:
    contract = ValidatedTau1Contract()
    assert contract.pressure_reference_pa == pytest.approx(3387510.7199999993, rel=2.0e-15)
    assert contract.pressure_reference_pa == pressure_reference_pa(contract.dx_m, contract.dt_s)
    for label, gauge in OUTLET_GAUGE_PRESSURES_PA.items():
        absolute = contract.outlet_absolute_pressures_pa[label]
        assert gauge_pressure_pa(absolute, contract.pressure_reference_pa) == pytest.approx(gauge, abs=5.0e-10)


def test_old_23622_reference_is_invalid_for_tau1_and_explains_rho_0007() -> None:
    old_reference = pressure_reference_pa(BASE_DX_M, HISTORICAL_DT_S)
    assert old_reference == pytest.approx(HISTORICAL_FIXED_REFERENCE_PRESSURE_PA, rel=2.0e-15)
    contract = ValidatedTau1Contract()
    bad_density = contract.lattice_density(HISTORICAL_FIXED_REFERENCE_PRESSURE_PA)
    assert bad_density == pytest.approx(0.006973356567857592, rel=2.0e-15)
    assert not 0.9 <= bad_density <= 1.1


def test_controller_target_uses_mass_flow_density_dt_over_dx_cubed() -> None:
    contract = ValidatedTau1Contract()
    expected = TARGET_VOLUME_FLOW_M3_S * contract.dt_s / contract.dx_m**3
    assert contract.target_lattice_flux == pytest.approx(expected, rel=2.0e-15)
    assert contract.target_lattice_flux == target_lattice_flux(
        contract.target_mass_flow_kg_s, contract.dx_m, contract.dt_s
    )


def test_invalid_formula_inputs_fail_fast() -> None:
    with pytest.raises(ValueError):
        tau1_time_step_s(0.0)
    with pytest.raises(ValueError):
        pressure_reference_pa(BASE_DX_M, 0.0)
    with pytest.raises(ValueError):
        target_lattice_flux(0.0, BASE_DX_M, 1.0)
    assert math.isfinite(ValidatedTau1Contract().pressure_reference_pa)

