"""Production-only QC gates for validated Tau1 CFD artifacts."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import pyvista as pv

from .physical_port_flux import FLUX_DEFINITION
from .validated_contract import FULL_TIMESTEP_REFEREE_REVISION


REQUIRED_VTU_ARRAYS = {
    "velocity_phy",
    "velocity_magnitude_m_s",
    "velocity_magnitude_mm_s",
    "pressure_gauge_pa",
    "pressure_absolute_solver_pa",
    "rho_lattice",
}


def physical_flow_balance_qc(flux: Mapping[str, Any], *, gate: float = 0.01) -> dict[str, Any]:
    """Gate physical aperture flow; boundary-PDF flux is not accepted as Q."""

    definition = flux.get("flux_definition")
    closure = float(flux["closure"])
    ports = flux["ports"]
    finite = all(
        np.isfinite(float(ports[label]["physical_q_m3_s"]))
        for label in ("inlet", "outlet_01", "outlet_02", "outlet_03")
    )
    checks = {
        "physical_flux_definition": definition == FLUX_DEFINITION,
        "finite": bool(finite),
        "positive_inlet": float(flux["Qin_m3_s"]) > 0.0,
        "positive_outlets": all(
            float(ports[label]["physical_q_m3_s"]) > 0.0
            for label in ("outlet_01", "outlet_02", "outlet_03")
        ),
        "closure": closure <= gate,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "flux_definition": definition,
        "closure": closure,
        "gate": gate,
        "checks": checks,
    }


def accepted_steady_qc(steady: Mapping[str, Any]) -> dict[str, Any]:
    """Require every promoted physical-time gate; no solver-message shortcut."""

    gates = dict(steady["gates"])
    return {
        "status": "PASS" if gates and all(gates.values()) else "FAIL",
        "acceptance_kind": "PHYSICAL_TIME_CANDIDATE_PLUS_CONFIRMATION",
        "built_in_steady_message_is_sufficient": False,
        "gates": gates,
        "failed_gates": [name for name, passed in gates.items() if not passed],
    }


def full_v2_qc(evidence: Mapping[str, Any], *, gate: float = 1.0e-8) -> dict[str, Any]:
    residual = float(evidence["referee"]["R_full_one_step_identity"])
    return {
        "status": "PASS" if evidence.get("status") == "PASS" and residual <= gate else "FAIL",
        "oracle": FULL_TIMESTEP_REFEREE_REVISION,
        "residual": residual,
        "gate": gate,
        "boundary_write_only_is_final_gate": False,
    }


def vtu_contract_qc(grid: pv.DataSet) -> dict[str, Any]:
    arrays = set(grid.cell_data.keys())
    missing = sorted(REQUIRED_VTU_ARRAYS - arrays)
    finite = all(
        bool(np.all(np.isfinite(np.asarray(grid.cell_data[name]))))
        for name in REQUIRED_VTU_ARRAYS & arrays
    )
    velocity = np.asarray(grid.cell_data["velocity_phy"]) if "velocity_phy" in arrays else np.empty((0, 3))
    checks = {
        "required_arrays": not missing,
        "velocity_three_components": velocity.ndim == 2 and velocity.shape[1:] == (3,),
        "all_finite": finite,
        "nonempty": grid.n_cells > 0,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "missing_arrays": missing,
        "cell_arrays": sorted(arrays),
        "cell_count": grid.n_cells,
        "point_count": grid.n_points,
    }
