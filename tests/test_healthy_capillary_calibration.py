from __future__ import annotations

import math

from utils.cfd_flow.healthy_capillary_calibration import (
    CHECKPOINT_INTERVAL,
    EARLIEST_FINAL_STEADY_ITERATION,
    EXPECTED_RESUME_SHA256,
    PRESSURE_THRESHOLD_PA,
    RESUME_ITERATION,
    SMOOTH_INLET_AREA_M2,
    STATUS_MASS_INCOMPLETE,
    STATUS_OUTLETS_UNRESOLVED,
    STATUS_PASS,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    VELOCITY_THRESHOLD_M_S,
    calculate_healthy_targets,
    classify_healthy_calibration,
    generate_healthy_calibration_lua,
    healthy_literature_reference_contract,
    healthy_lua_contract,
)


def test_literature_contract_has_five_primary_records() -> None:
    contract = healthy_literature_reference_contract()
    assert contract["status"] == "PASS"
    assert len(contract["studies"]) == 5
    assert all(item["doi"] and item["url"] and item["exact_value_used"] for item in contract["studies"])
    assert contract["continuum_mapping"]["status"] == "MODEL_CALIBRATION_ASSUMPTION"


def test_healthy_target_is_recomputed_from_smooth_area() -> None:
    targets = calculate_healthy_targets()
    primary = targets["candidates"]["PRIMARY"]
    assert math.isclose(TARGET_Q_M3_S, 3.5e-4 * SMOOTH_INLET_AREA_M2, rel_tol=0.0, abs_tol=1e-30)
    assert math.isclose(TARGET_Q_M3_S, 2.7369132390905703e-15, rel_tol=1e-15)
    assert math.isclose(TARGET_MASS_FLOW_KG_S, 2.890180380479642e-12, rel_tol=1e-15)
    assert math.isclose(primary["volume_flow_nl_min"], 0.1642147943454342, rel_tol=1e-15)
    assert targets["discrete_cell_area_proxy_is_physical_area"] is False
    assert targets["selected"] == "PRIMARY"


def test_healthy_lua_contract_is_restart_only_and_physics_frozen() -> None:
    lua = generate_healthy_calibration_lua(
        mesh_wsl="/frozen/seeder/mesh",
        resume_header_wsl="restart/a3274_lastHeader.lua",
    )
    contract = healthy_lua_contract(lua, resume_header_wsl="restart/a3274_lastHeader.lua")
    assert contract["status"] == "PASS"
    assert f"min = {{ iter = {RESUME_ITERATION + CHECKPOINT_INTERVAL} }}" in lua
    assert f"nvals = 100" in lua
    assert f"threshold = {PRESSURE_THRESHOLD_PA:.17g}" in lua
    assert f"threshold = {VELOCITY_THRESHOLD_M_S:.17g}" in lua
    assert "adaptive_flux_pressure" in lua
    assert "mfr_eq" not in lua
    assert "read = 'restart/a3274_lastHeader.lua'" in lua
    assert EARLIEST_FINAL_STEADY_ITERATION == 169326
    assert EXPECTED_RESUME_SHA256 == "c911dd8085fea971590758bb0732cc2245f39f9c3720a283215237de0280d3d0"


def _exact(*, balance: float = 0.0, backflow: bool = False) -> dict:
    return {
        "all_finite": True, "minimum_pdf": 0.01, "maximum_lattice_speed": 0.001,
        "inlet": {"target_relative_error": 0.0},
        "healthy_velocity_compatibility": "PASS",
        "healthy_wss_proxy_compatibility": "PASS",
        "instantaneous_mass_balance_relative_error": balance,
        "significant_outlet_02_backflow": backflow,
    }


def _mass(*, ratio: float = 0.0, crosscheck: str = "PASS") -> dict:
    return {"accumulation_to_healthy_target_ratio": ratio, "pdf_pressure_crosscheck_status": crosscheck}


def test_final_classification_requires_all_three_steady_gates() -> None:
    assert classify_healthy_calibration(
        official_steady=True, additional_iterations=10_000,
        exact=_exact(), mass=_mass(), source_unchanged=True,
    )[0] == STATUS_PASS
    assert classify_healthy_calibration(
        official_steady=True, additional_iterations=10_000,
        exact=_exact(), mass=_mass(ratio=0.02), source_unchanged=True,
    )[0] == STATUS_MASS_INCOMPLETE
    assert classify_healthy_calibration(
        official_steady=True, additional_iterations=10_000,
        exact=_exact(backflow=True), mass=_mass(), source_unchanged=True,
    )[0] == STATUS_OUTLETS_UNRESOLVED
