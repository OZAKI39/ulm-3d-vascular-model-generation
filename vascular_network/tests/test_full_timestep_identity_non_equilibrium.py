from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.full_timestep_mass_referee import replay_full_timestep
from utils.cfd_flow.musubi_boundary_mass_referee import load_mesh_contract
from utils.cfd_flow.restart_decode import read_restart_pdf
from utils.cfd_flow.validated_contract import OUTLET_GAUGE_PRESSURES_PA, RHO0_KG_M3, TARGET_MASS_FLOW_KG_S, ValidatedTau1Contract


EXPECTED_CELLS = 182_320


def test_late_non_equilibrium_3117927_transition_closes() -> None:
    root = Path(__file__).resolve().parents[1]
    run_root = root / "outputs/cfd_flow/healthy_mouse_capillary_tau1_base_anchor003274_20260830"
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
    mesh = load_mesh_contract(
        root / "outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh",
        expected_cells=EXPECTED_CELLS,
    )

    replay = replay_full_timestep(
        start,
        end,
        mesh,
        dx_m=2.0e-7,
        dt_s=ValidatedTau1Contract().dt_s,
        density_kg_m3=RHO0_KG_M3,
        target_mass_flow_kg_s=TARGET_MASS_FLOW_KG_S,
        outlet_pressures_pa={
            label: 23622.320128 + gauge
            for label, gauge in OUTLET_GAUGE_PRESSURES_PA.items()
        },
    )

    assert replay["R_full_one_step_identity"] <= 1.0e-8
    assert replay["R_full_one_step_identity"] < replay["R_boundary_only"]
    assert replay["runtime_solid_count"] == 47
    assert replay["compute_target_count"] == EXPECTED_CELLS
