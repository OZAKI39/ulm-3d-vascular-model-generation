"""Boundary-condition correction for the added numerical straight tubes."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Iterable

import numpy as np

from .io import BoundaryInput, SurfacePrepareError
from .types import BoundarySurfaceResult


CORRECTION_ROLE = "NUMERICAL_ARTIFICIAL_EXTENSION_CORRECTION"


def calculate_pressure_corrections(
    boundaries: Iterable[BoundaryInput],
    results: Iterable[BoundarySurfaceResult],
    *,
    dynamic_viscosity_pa_s: float,
    allow_negative_gauge_pressure: bool,
) -> list[dict[str, Any]]:
    """Apply Poiseuille resistance only to each artificial outlet extension."""

    if not np.isfinite(dynamic_viscosity_pa_s) or dynamic_viscosity_pa_s <= 0:
        raise SurfacePrepareError("INVALID_SAVED_DYNAMIC_VISCOSITY")
    result_by_index = {item.boundary_index: item for item in results}
    rows: list[dict[str, Any]] = []
    for boundary in boundaries:
        result = result_by_index.get(boundary.index)
        if result is None:
            raise SurfacePrepareError(f"Missing surface result for {boundary.port_id}")
        area_m2 = result.actual_cap_area_um2 * 1.0e-12
        radius_m = np.sqrt(area_m2 / np.pi)
        length_m = boundary.extension_length_um * 1.0e-6
        common = {
            "port_id": boundary.port_id,
            "boundary_origin": boundary.boundary_origin,
            "role": boundary.role,
            "P_original_1D_pa": boundary.pressure_original_pa,
            "Q_expected_1D_m3_s": boundary.expected_flow_m3_s,
            "extension_length_um": boundary.extension_length_um,
            "actual_cap_area_um2": result.actual_cap_area_um2,
            "equivalent_radius_um": result.equivalent_radius_um,
            "actual_cap_area_m2": area_m2,
            "extension_pressure_correction_role": CORRECTION_ROLE,
        }
        if boundary.role == "ASSUMED_INLET":
            rows.append(
                {
                    **common,
                    "extension_resistance_pa_s_m3": None,
                    "predicted_extension_pressure_drop_pa": None,
                    "P_solver_boundary_pa": None,
                    "pressure_correction_applied": False,
                    "Q_solver_m3_s": boundary.expected_flow_m3_s,
                }
            )
            continue
        resistance = 8.0 * dynamic_viscosity_pa_s * length_m / (
            np.pi * radius_m**4
        )
        pressure_drop = resistance * boundary.expected_flow_m3_s
        solver_pressure = boundary.pressure_original_pa - pressure_drop
        if solver_pressure < 0 and not allow_negative_gauge_pressure:
            raise SurfacePrepareError("NEGATIVE_GAUGE_PRESSURE_NOT_ALLOWED")
        rows.append(
            {
                **common,
                "extension_resistance_pa_s_m3": float(resistance),
                "predicted_extension_pressure_drop_pa": float(pressure_drop),
                "P_solver_boundary_pa": float(solver_pressure),
                "pressure_correction_applied": True,
                "Q_solver_m3_s": None,
            }
        )
    return rows


def build_extended_boundary_conditions(
    original: dict[str, Any],
    boundaries: Iterable[BoundaryInput],
    corrections: Iterable[dict[str, Any]],
) -> dict[str, Any]:
    """Create traceable solver-plane BCs without changing the original BC object."""

    source = deepcopy(original)
    correction_by_id = {str(row["port_id"]): row for row in corrections}
    boundary_by_id = {item.port_id: item for item in boundaries}
    inlet_source = source["inlet"]
    inlet_boundary = boundary_by_id[str(inlet_source["port_id"])]
    inlet_correction = correction_by_id[inlet_boundary.port_id]
    inlet = {
        "port_id": inlet_boundary.port_id,
        "role": inlet_boundary.role,
        "boundary_origin": inlet_boundary.boundary_origin,
        "type": "VOLUMETRIC_FLOW_RATE",
        "flow_rate_m3_s": inlet_source["flow_rate_m3_s"],
        "profile": inlet_source["profile"],
        "original_plane_center_um": inlet_boundary.center_um.tolist(),
        "solver_plane_center_um": inlet_boundary.extension_end_um.tolist(),
        "actual_solver_cap_area_m2": inlet_correction["actual_cap_area_m2"],
        "pressure_correction_applied": False,
    }
    outlets: list[dict[str, Any]] = []
    for outlet_source in source["outlets"]:
        boundary = boundary_by_id[str(outlet_source["port_id"])]
        correction = correction_by_id[boundary.port_id]
        outlets.append(
            {
                "port_id": boundary.port_id,
                "role": boundary.role,
                "boundary_origin": boundary.boundary_origin,
                "type": "PRESSURE_DIRICHLET",
                "P_original_1D_pa": correction["P_original_1D_pa"],
                "P_solver_boundary_pa": correction["P_solver_boundary_pa"],
                "expected_1D_flow_m3_s": correction["Q_expected_1D_m3_s"],
                "predicted_extension_pressure_drop_pa": correction[
                    "predicted_extension_pressure_drop_pa"
                ],
                "original_plane_center_um": boundary.center_um.tolist(),
                "solver_plane_center_um": boundary.extension_end_um.tolist(),
                "actual_solver_cap_area_m2": correction["actual_cap_area_m2"],
                "pressure_correction_applied": True,
            }
        )
    return {
        "method": source.get("method"),
        "source_boundary_conditions_preserved": True,
        "extension_pressure_correction_role": CORRECTION_ROLE,
        "is_physiological_outlet_model": False,
        "is_resistance_boundary_condition": False,
        "is_windkessel": False,
        "fluid": deepcopy(source.get("fluid")),
        "flow_model": deepcopy(source.get("flow_model")),
        "pressure_reference": source.get("pressure_reference"),
        "inlet": inlet,
        "outlets": outlets,
        "wall": deepcopy(source.get("wall")),
    }
