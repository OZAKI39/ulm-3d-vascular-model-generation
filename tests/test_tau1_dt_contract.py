import pytest

from utils.cfd_flow.musubi_pressure_bc_benchmark import NU_M2_S
from utils.cfd_flow.musubi_wall_force_diagnostics import (
    lattice_relaxation_contract,
    tau_one_time_step_s,
)
from utils.cfd_flow.periodic_pipe_force import CASES


def test_tau1_timestep_inverts_musubi_relaxation_exactly() -> None:
    grid = CASES["axis_n27"]
    dt = tau_one_time_step_s(nu_phy_m2_s=NU_M2_S, dx_m=grid.dx_m)
    contract = lattice_relaxation_contract(
        nu_phy_m2_s=NU_M2_S, dx_m=grid.dx_m, dt_s=dt
    )
    assert contract["nu_lattice"] == pytest.approx(1.0 / 6.0)
    assert contract["tau"] == pytest.approx(1.0)
    assert contract["omega"] == pytest.approx(1.0)
