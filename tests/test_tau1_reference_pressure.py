from __future__ import annotations

import math

from utils.cfd_flow.tau1_base import (
    OLD_DT_S,
    OLD_PRESSURE_REFERENCE_PA,
    OUTLET_GAUGE_PRESSURE_PA,
    TARGET_Q_M3_S,
    Tau1BaseRuntimeContract,
    Tau1ReferencePressureContract,
)
from utils.cfd_flow.tau1_reference_pressure import (
    OLD_BASE_CLASSIFICATION,
    generate_smoke_lua,
    smoke_lua_contract,
)


def test_old_pressure_reference_equals_old_unit_density_pressure() -> None:
    old = Tau1ReferencePressureContract(dt_s=OLD_DT_S)
    assert math.isclose(
        old.unit_density_pressure_pa,
        OLD_PRESSURE_REFERENCE_PA,
        rel_tol=2.0e-15,
        abs_tol=1.0e-10,
    )


def test_tau1_reference_pressure_formula() -> None:
    contract = Tau1ReferencePressureContract()
    expected = contract.rho0_kg_m3 * contract.cs2 * contract.dx_m**2 / contract.dt_s**2
    assert math.isclose(
        contract.pressure_reference_pa, expected, rel_tol=2.0e-15
    )


def test_tau1_reference_density_is_one() -> None:
    contract = Tau1ReferencePressureContract()
    assert contract.lattice_density(contract.pressure_reference_pa) == 1.0


def test_pressure_reference_scales_as_dt_inverse_squared() -> None:
    old = Tau1ReferencePressureContract(dt_s=OLD_DT_S)
    new = Tau1ReferencePressureContract()
    assert math.isclose(
        new.pressure_reference_pa / old.pressure_reference_pa,
        (old.dt_s / new.dt_s) ** 2,
        rel_tol=2.0e-15,
    )


def test_outlet_gauge_pressures_preserved_under_reference_rescaling() -> None:
    contract = Tau1ReferencePressureContract()
    absolute = contract.outlet_absolute_pressures(OUTLET_GAUGE_PRESSURE_PA)
    for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items():
        assert math.isclose(
            absolute[label] - contract.pressure_reference_pa,
            gauge,
            rel_tol=0.0,
            abs_tol=5.0e-10,
        )


def test_tau1_lua_uses_dynamic_reference_pressure() -> None:
    runtime = Tau1BaseRuntimeContract()
    lua = generate_smoke_lua(runtime)
    assert f"pressure_reference_phy = {runtime.pressure_reference_pa:.17g}" in lua
    assert "format='asciiSpatial'" in lua
    assert "format='ascii'" not in lua
    assert smoke_lua_contract(lua, runtime)["status"] == "PASS"


def test_tau1_fresh_smoke_forbids_old_restart() -> None:
    lua = generate_smoke_lua()
    assert "read=" not in lua
    assert "read =" not in lua
    assert "velocityX=0.0" in lua


def test_pressure_output_uses_dynamic_gauge_reference() -> None:
    contract = Tau1ReferencePressureContract()
    assert math.isclose(
        contract.gauge_pressure(contract.pressure_reference_pa + 17.0),
        17.0,
        abs_tol=1.0e-10,
    )


def test_adaptive_flux_target_contract_at_rho_one() -> None:
    runtime = Tau1BaseRuntimeContract()
    expected = TARGET_Q_M3_S * runtime.dt_s / runtime.dx_m**3
    assert math.isclose(runtime.target_lattice_flux, expected, rel_tol=2.0e-15)


def test_old_misscaled_base_classification() -> None:
    contract = Tau1ReferencePressureContract()
    old_density = contract.lattice_density(OLD_PRESSURE_REFERENCE_PA)
    assert OLD_BASE_CLASSIFICATION.endswith("MIS_SCALED_REFERENCE_PRESSURE_OFFSET")
    assert math.isclose(old_density, 0.006973356567857592, rel_tol=2.0e-15)
