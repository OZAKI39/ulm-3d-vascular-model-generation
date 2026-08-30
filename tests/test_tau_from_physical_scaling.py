import pytest

from utils.cfd_flow.musubi_wall_force_diagnostics import lattice_relaxation_contract
from utils.cfd_flow.periodic_pipe_force import CASES
from utils.cfd_flow.musubi_pressure_bc_benchmark import NU_M2_S


def test_diffusive_scaling_keeps_tau_constant_across_baseline_grids() -> None:
    contracts = [
        lattice_relaxation_contract(
            nu_phy_m2_s=NU_M2_S, dx_m=CASES[name].dx_m, dt_s=CASES[name].dt_s
        )
        for name in ("axis_n16", "axis_n20", "axis_n27")
    ]
    assert [row["nu_lattice"] for row in contracts] == pytest.approx([1.995849609375] * 3)
    assert [row["tau"] for row in contracts] == pytest.approx([6.487548828125] * 3)
    assert [row["omega"] for row in contracts] == pytest.approx([0.15414142174387535] * 3)
