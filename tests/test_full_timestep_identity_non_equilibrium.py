from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.full_timestep_mass_referee import replay_full_timestep
from utils.cfd_flow.musubi_boundary_mass_referee import load_mesh_contract
from utils.cfd_flow.restart_decode import read_restart_pdf
from utils.cfd_flow.tau1_base import (
    EXPECTED_CELLS,
    OUTLET_GAUGE_PRESSURE_PA,
    _mesh_path,
    _run_root,
    historical_tau1_runtime_contract,
)


def test_late_non_equilibrium_3117927_transition_closes() -> None:
    root = Path(__file__).resolve().parents[1]
    run_root = _run_root(root)
    start = read_restart_pdf(
        run_root / "dense_diagnostic" / "restart_endpoint" / "tau1_base_6.357E-03.lsb",
        n_elems=EXPECTED_CELLS,
        n_components=19,
    )
    end = read_restart_pdf(
        run_root
        / "dense_diagnostic_corrected"
        / "attempt_1_8_steps"
        / "restart"
        / "tau1_base_3117928.lsb",
        n_elems=EXPECTED_CELLS,
        n_components=19,
    )
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    contract = historical_tau1_runtime_contract()

    replay = replay_full_timestep(
        start,
        end,
        mesh,
        dx_m=contract.dx_m,
        dt_s=contract.dt_s,
        density_kg_m3=contract.rho_kg_m3,
        target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
        outlet_pressures_pa={
            label: contract.pressure_reference_pa + gauge
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        },
    )

    assert replay["R_full_one_step_identity"] <= 1.0e-8
    assert replay["R_full_one_step_identity"] < replay["R_boundary_only"]
    assert replay["runtime_solid_count"] == 47
    assert replay["compute_target_count"] == EXPECTED_CELLS
