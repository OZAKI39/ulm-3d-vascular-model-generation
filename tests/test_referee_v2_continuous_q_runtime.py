from __future__ import annotations

import math

from utils.cfd_flow.musubi_boundary_mass_referee import (
    REFEREE_REVISION_NEW,
    conservation_identity_residual,
)
from utils.cfd_flow.tau1_base import (
    EXPECTED_CELLS,
    Tau1BaseRuntimeContract,
    generate_referee_lua,
    referee_lua_contract,
)


def test_referee_v2_uses_frozen_mesh_new_dt_and_target_normalization() -> None:
    contract = Tau1BaseRuntimeContract()
    assert REFEREE_REVISION_NEW == "MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2"
    assert contract.expected_cells == EXPECTED_CELLS == 182_320
    lua = generate_referee_lua(contract)
    assert referee_lua_contract(lua, contract)["status"] == "PASS"
    assert f"dt = {contract.dt_s:.17g}" in lua
    target = contract.target_lattice_flux
    assert math.isclose(target, 6.974804380964758e-4, rel_tol=1.0e-15)
    assert conservation_identity_residual(target + 1.0e-12, target, target) < 1.0e-8
