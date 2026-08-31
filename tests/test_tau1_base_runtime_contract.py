from __future__ import annotations

import math

from utils.cfd_flow.tau1_base import (
    TARGET_MASS_FLOW_KG_S,
    Tau1BaseRuntimeContract,
    base_lua_contract,
    generate_base_segment_lua,
)


def test_tau1_runtime_physics_and_boundary_contract() -> None:
    contract = Tau1BaseRuntimeContract()
    assert contract.dx_m == 2.0e-7
    assert math.isclose(contract.dt_s, contract.dx_m**2 / (6.0 * contract.nu_m2_s))
    assert math.isclose(contract.nu_lattice, 1.0 / 6.0)
    assert abs(contract.tau - 1.0) <= 1.0e-12
    assert abs(contract.omega - 1.0) <= 1.0e-12
    assert contract.target_mass_flow_kg_s == TARGET_MASS_FLOW_KG_S
    lua = generate_base_segment_lua(
        maximum_iterations=239_502,
        segment_wsl="/tmp/tau1_segment",
        restart_header_wsl=None,
        contract=contract,
    )
    assert base_lua_contract(
        lua,
        maximum_iterations=239_502,
        restart_header_wsl=None,
        contract=contract,
    )["status"] == "PASS"
    confirmation = generate_base_segment_lua(
        maximum_iterations=contract.hard_max_iterations,
        segment_wsl="/tmp/tau1_confirmation",
        restart_header_wsl="/tmp/tau1_restart.lua",
        restart_first_iteration=598_755,
        restart_interval=119_751,
        contract=contract,
    )
    assert base_lua_contract(
        confirmation,
        maximum_iterations=contract.hard_max_iterations,
        restart_header_wsl="/tmp/tau1_restart.lua",
        restart_first_iteration=598_755,
        restart_interval=119_751,
        contract=contract,
    )["status"] == "PASS"


def test_dense_restart_uses_iteration_filename_key_loaded_by_tem() -> None:
    lua = generate_base_segment_lua(
        maximum_iterations=3_117_959,
        segment_wsl="/tmp/tau1_dense",
        restart_header_wsl="/tmp/tau1_base_header.lua",
        restart_first_iteration=3_117_928,
        restart_interval=1,
        restart_write_wsl="/tmp/tau1_dense/restart/",
        restart_use_iteration_filename=True,
        stop_file_wsl="/tmp/tau1_dense/stop",
        wall_clock_seconds=1_200,
    )
    assert "timeformat={use_iter=true}" in lua
    assert "timeform={use_iter=true}" not in lua
    assert base_lua_contract(
        lua,
        maximum_iterations=3_117_959,
        restart_header_wsl="/tmp/tau1_base_header.lua",
        restart_first_iteration=3_117_928,
        restart_interval=1,
        restart_write_wsl="/tmp/tau1_dense/restart/",
        restart_use_iteration_filename=True,
    )["status"] == "PASS"
