"""Pure Cartesian field reconstruction helpers for decoded Musubi PDFs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv


ALIGNMENT_TOLERANCE_FRACTION = 1.0e-6


@dataclass(frozen=True, slots=True)
class LatticeMapping:
    cell_indices: np.ndarray
    maximum_alignment_error_m: float
    duplicate_cell_count: int
    unique_cell_count: int


_HEX_CORNERS = np.asarray(
    (
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    ),
    dtype=np.int64,
)


def quantize_cell_centers(
    coordinates_m: np.ndarray,
    *,
    origin_m: np.ndarray,
    dx_m: float,
) -> LatticeMapping:
    """Map uniform cell centers to exact integer Cartesian indices."""

    coordinates = np.asarray(coordinates_m, dtype=np.float64)
    origin = np.asarray(origin_m, dtype=np.float64).reshape(3)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or dx_m <= 0.0:
        raise ValueError("coordinates must be Nx3 and dx must be positive")
    scaled = (coordinates - origin) / float(dx_m) - 0.5
    indices = np.rint(scaled).astype(np.int64)
    reconstructed = origin + (indices.astype(np.float64) + 0.5) * float(dx_m)
    errors = np.linalg.norm(coordinates - reconstructed, axis=1)
    maximum_error = float(np.max(errors)) if len(errors) else 0.0
    if maximum_error > ALIGNMENT_TOLERANCE_FRACTION * float(dx_m):
        raise ValueError("cell centers are not aligned to the uniform lattice")
    unique_count = int(len(np.unique(indices, axis=0)))
    return LatticeMapping(
        cell_indices=indices,
        maximum_alignment_error_m=maximum_error,
        duplicate_cell_count=int(len(indices) - unique_count),
        unique_cell_count=unique_count,
    )


def reconstruct_hexahedral_field(
    *,
    mapping: LatticeMapping,
    origin_m: np.ndarray,
    dx_m: float,
    pressure_pa: np.ndarray,
    velocity_m_s: np.ndarray,
    pressure_reference_pa: float,
) -> pv.UnstructuredGrid:
    """Build a shared-vertex VTK hexahedral grid with physical fields."""

    indices = np.asarray(mapping.cell_indices, dtype=np.int64)
    pressure = np.asarray(pressure_pa, dtype=np.float64).reshape(-1)
    velocity = np.asarray(velocity_m_s, dtype=np.float64)
    if velocity.shape != (len(indices), 3) or pressure.shape != (len(indices),):
        raise ValueError("field shape does not match lattice mapping")
    corner_indices = (indices[:, None, :] + _HEX_CORNERS[None, :, :]).reshape(-1, 3)
    points_index, inverse = np.unique(corner_indices, axis=0, return_inverse=True)
    points = np.asarray(origin_m, dtype=np.float64).reshape(3) + points_index * float(dx_m)
    cells = np.column_stack(
        (np.full(len(indices), 8, dtype=np.int64), inverse.reshape(-1, 8))
    ).reshape(-1)
    cell_types = np.full(len(indices), int(pv.CellType.HEXAHEDRON), dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells, cell_types, points)
    grid.cell_data["velocity_phy"] = velocity
    grid.cell_data["pressure_phy"] = pressure
    grid.cell_data["pressure_gauge_pa"] = pressure - float(pressure_reference_pa)
    return grid


def build_proteus_metadata(
    source_flow_vtu: Path,
    *,
    inlet_equivalent_diameter_m: float,
) -> dict[str, Any]:
    """Describe decoded fields using the existing PROTEUS naming contract."""

    return {
        "lengthUnit": 1.0,
        "velocityUnit": 1.0,
        "velocityField": "velocity_phy",
        "pressureField": "pressure_phy",
        "pressureGaugeField": "pressure_gauge_pa",
        "inletEquivalentDiameterM": float(inlet_equivalent_diameter_m),
        "sourceFlowVtu": str(Path(source_flow_vtu)),
    }
