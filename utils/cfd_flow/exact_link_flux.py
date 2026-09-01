"""Exact, zero-run Musubi boundary-link flux audit from one frozen restart.

The implementation mirrors the pinned D3Q19 PULL boundary path.  It never
launches Seeder, Musubi, or harvesting and treats all solver artefacts as
read-only evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np

from .restart_decode import (
    D3Q19_DIRECTIONS,
    D3Q19_WEIGHTS,
)


SOURCE_BRANCH = "codex/cfd-flow-musubi-recovery-20260828"
CURRENT_SYNCED_BASE_COMMIT = "02660620d9bead5356d931f0c4fa65bd3efa3491"
DIRECT_FIELD_RUN = "musubi_direct_restart_field_anchor003274_20260829_000544"
FROZEN_SEEDER_RUN = "musubi_recovery_anchor003274_20260828_162530"
FROZEN_STEADY_RUN = "musubi_project_steady_confirmation_anchor003274_20260828_225334"
AUDIT_PREFIX = "exact_link_flux_audit_anchor003274"
AUDIT_REVISION = "ZERO_RUN_EXACT_MUSUBI_BOUNDARY_LINK_FLUX_V1"

MUSUBI_SOURCE_COMMIT = "4e8b277b66226277171ef93bf054d21270812793"
MUSUBI_SCHEME_COMMIT = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
TREELM_SOURCE_COMMIT = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
MUSUBI_EXECUTABLE_SHA256 = (
    "a005b4f00bd45df0339adc22460f251c3f300f967ff746c1cd43fa5ad7c07e88"
)

MAXIMUM_RELATIVE_ERROR = 0.01
CS2 = 1.0 / 3.0

NOT_IDENTIFIABLE = "CFD_FLOW_EXACT_LINK_FLUX_NOT_IDENTIFIABLE_FROM_SINGLE_RESTART"
MFR_TARGET_PASS = "MFR_EQ_TARGET_PASS"
MFR_TARGET_MISMATCH = "MFR_EQ_OBLIQUE_PORT_TARGET_MISMATCH_CONFIRMED"
MASS_BALANCE_PASS = "EXACT_LINK_MASS_BALANCE_PASS"
NEXT_NOT_IDENTIFIABLE = "DESIGN MINIMAL ONE-STEP INSTRUMENTED FLUX TRACKING"
NEXT_MISMATCH = "REVIEW MFR_EQ OBLIQUE-PORT AREA AND RAW DISCRETE-NORMAL SEMANTICS"
NEXT_TARGET_PASS = "PROMOTE EXACT PDF-LINK FLUX QC THEN PROCEED TO GRID CONVERGENCE"

INVERSE_DIRECTIONS = np.asarray(
    [
        int(np.flatnonzero(np.all(D3Q19_DIRECTIONS == -direction, axis=1))[0])
        for direction in D3Q19_DIRECTIONS
    ],
    dtype=np.int8,
)
NORMAL_WEIGHTS = np.asarray(
    (*(4 for _ in range(6)), *(2 for _ in range(12))), dtype=np.int8
)
GLOBAL_NORMALS = {
    "outlet_01": np.asarray((0, -1, 0), dtype=np.int8),
    "outlet_02": np.asarray((1, 1, 0), dtype=np.int8),
    "outlet_03": np.asarray((-1, 0, 0), dtype=np.int8),
    "inlet": np.asarray((0, -1, -1), dtype=np.int8),
}
EXPECTED_STENCIL_COUNTS = {
    "inlet": 311,
    "outlet_01": 169,
    "outlet_02": 179,
    "outlet_03": 201,
}
EXPECTED_PRESSURE_VALID_COUNTS = {"outlet_01": 169, "outlet_02": 173, "outlet_03": 201}


@dataclass(frozen=True, slots=True)
class BoundaryReconstruction:
    label: str
    boundary_id: int
    property_rows: np.ndarray
    cell_indices: np.ndarray
    outward_masks: np.ndarray
    incoming_masks: np.ndarray
    raw_normals: np.ndarray
    normal_indices: np.ndarray


def mfr_eq_area_proxy(element_count: int, dx_m: float) -> float:
    """Mirror ``sum(nElems_totalLevel * dxLvl**2)`` for one level."""

    return int(element_count) * float(dx_m) ** 2


def closest_discrete_direction(raw_normal: np.ndarray) -> int:
    """Mirror tem_determine_discreteVector and return a zero-based index."""

    vector = np.asarray(raw_normal, dtype=np.float64).reshape(3)
    length = float(np.linalg.norm(vector))
    if length == 0.0:
        raise ValueError("A boundary normal cannot be zero")
    prevailing = D3Q19_DIRECTIONS[:18].astype(np.float64)
    prevailing /= np.linalg.norm(prevailing, axis=1)[:, None]
    return int(np.argmax(prevailing @ (vector / length)))


def reconstruct_boundary(
    boundary_ids: np.ndarray,
    property_element_indices: np.ndarray,
    *,
    label: str,
    boundary_id: int,
    allow_zero_normals: bool = False,
) -> BoundaryReconstruction:
    """Mirror countnBnds/assignBCList for the active D3Q19 QQN=18."""

    ids = np.asarray(boundary_ids, dtype=np.int64)
    property_indices = np.asarray(property_element_indices, dtype=np.int64).reshape(-1)
    if ids.shape[0] != len(property_indices) or ids.shape[1] < 18:
        raise ValueError("Boundary IDs and property-element mapping are inconsistent")
    active = ids[:, :18] == int(boundary_id)
    rows = np.flatnonzero(np.any(active, axis=1)).astype(np.int64)
    outward = active[rows]
    incoming = np.zeros_like(outward)
    raw = np.zeros((len(rows), 3), dtype=np.int64)
    normal_indices = np.empty(len(rows), dtype=np.int64)
    for local, mask in enumerate(outward):
        outgoing_directions = np.flatnonzero(mask)
        incoming[local, INVERSE_DIRECTIONS[outgoing_directions]] = True
        raw[local] = -np.sum(
            NORMAL_WEIGHTS[outgoing_directions, None]
            * D3Q19_DIRECTIONS[outgoing_directions].astype(np.int64),
            axis=0,
        )
        if np.all(raw[local] == 0) and allow_zero_normals:
            # Adaptive-flux/pressure mesh QC consumes the exact boundary-link
            # masks, not normalInd. A zero aggregate at a coarse Cartesian rim
            # is represented explicitly instead of inventing a direction.
            # The default remains strict for source-proven V2 audits.
            normal_indices[local] = -1
        else:
            normal_indices[local] = closest_discrete_direction(raw[local])
    return BoundaryReconstruction(
        label=label,
        boundary_id=int(boundary_id),
        property_rows=rows,
        cell_indices=property_indices[rows],
        outward_masks=outward,
        incoming_masks=incoming,
        raw_normals=raw,
        normal_indices=normal_indices,
    )


def equilibrium_pdf(rho: np.ndarray | float, velocity: np.ndarray) -> np.ndarray:
    """Strict standard D3Q19 equilibrium used by pinned get_pdfEq_d3q19."""

    density = np.asarray(rho, dtype=np.float64)
    vel = np.asarray(velocity, dtype=np.float64)
    scalar = vel.ndim == 1
    if scalar:
        vel = vel.reshape(1, 3)
    if vel.ndim != 2 or vel.shape[1] != 3:
        raise ValueError("velocity must have shape (3,) or (n, 3)")
    if density.ndim == 0:
        density = np.full(len(vel), float(density), dtype=np.float64)
    else:
        density = density.reshape(-1)
    if len(density) != len(vel):
        raise ValueError("rho and velocity lengths differ")
    cu = vel @ D3Q19_DIRECTIONS.astype(np.float64).T
    u2 = np.sum(vel * vel, axis=1)
    result = (
        density[:, None]
        * D3Q19_WEIGHTS[None, :]
        * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2[:, None])
    )
    return result[0] if scalar else result


def build_coordinate_lookup(cell_ijk: np.ndarray) -> dict[tuple[int, int, int], int]:
    coordinates = np.asarray(cell_ijk, dtype=np.int64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("cell_ijk must have shape (n, 3)")
    lookup = {
        tuple(int(value) for value in row): index
        for index, row in enumerate(coordinates)
    }
    if len(lookup) != len(coordinates):
        raise ValueError("cell_ijk contains duplicates")
    return lookup


def pull_fetch_pdfs(
    pdf: np.ndarray,
    cell_ijk: np.ndarray,
    target_indices: np.ndarray,
    *,
    coordinate_lookup: dict[tuple[int, int, int], int] | None = None,
) -> np.ndarray:
    """Mirror PULL FETCH, including local inverse-direction bounce-back."""

    values = np.asarray(pdf, dtype=np.float64)
    coordinates = np.asarray(cell_ijk, dtype=np.int64)
    targets = np.asarray(target_indices, dtype=np.int64).reshape(-1)
    lookup = coordinate_lookup or build_coordinate_lookup(coordinates)
    fetched = np.empty((len(targets), 19), dtype=np.float64)
    for row, target in enumerate(targets):
        for direction_index, direction in enumerate(D3Q19_DIRECTIONS):
            source_coord = tuple(
                int(value) for value in coordinates[target] - direction
            )
            source = lookup.get(source_coord)
            if source is None:
                fetched[row, direction_index] = values[
                    target, INVERSE_DIRECTIONS[direction_index]
                ]
            else:
                fetched[row, direction_index] = values[source, direction_index]
    return fetched


def velocity_from_pdf(pdf: np.ndarray) -> np.ndarray:
    values = np.asarray(pdf, dtype=np.float64)
    rho = np.sum(values, axis=1)
    return (values @ D3Q19_DIRECTIONS.astype(np.float64)) / rho[:, None]


def signed_mass_balance(
    inlet_kg_s: float, outlets_kg_s: Iterable[float]
) -> dict[str, float]:
    outlet_values = tuple(float(value) for value in outlets_kg_s)
    outlet_sum = float(sum(outlet_values))
    error = abs(float(inlet_kg_s) - outlet_sum) / abs(float(inlet_kg_s))
    return {"outlet_signed_sum_kg_s": outlet_sum, "relative_error": error}


def classify_exact_flux(
    inlet_relative_error: float, mass_balance_error: float, outlet02: float
) -> dict[str, Any]:
    inlet_pass = float(inlet_relative_error) <= MAXIMUM_RELATIVE_ERROR
    balance_pass = float(mass_balance_error) <= MAXIMUM_RELATIVE_ERROR
    status = MFR_TARGET_PASS if inlet_pass else MFR_TARGET_MISMATCH
    return {
        "status": status,
        "next": NEXT_TARGET_PASS if inlet_pass else NEXT_MISMATCH,
        "mfr_eq_target_classification": status,
        "exact_link_mass_balance": MASS_BALANCE_PASS if balance_pass else "FAIL",
        "outlet_02_backflow_confirmed": "YES"
        if balance_pass and outlet02 < 0.0
        else "UNRESOLVED",
        "mfr_eq_oblique_target_mismatch_confirmed": "NO" if inlet_pass else "YES",
    }
