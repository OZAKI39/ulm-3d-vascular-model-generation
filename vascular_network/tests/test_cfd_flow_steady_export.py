"""Tests for solver-free Cartesian field reconstruction helpers."""

from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

from utils.cfd_flow.steady_export import (
    build_proteus_metadata,
    quantize_cell_centers,
    reconstruct_hexahedral_field,
)


def test_cell_center_quantization_and_duplicate_detection():
    origin = np.array([1.0, 2.0, 3.0])
    dx = 0.25
    coordinates = origin + (np.array([[0, 0, 0], [1, 0, 0]]) + 0.5) * dx

    mapping = quantize_cell_centers(coordinates, origin_m=origin, dx_m=dx)
    duplicate = quantize_cell_centers(
        np.vstack((coordinates, coordinates[0])), origin_m=origin, dx_m=dx
    )

    assert mapping.maximum_alignment_error_m == pytest.approx(0.0)
    assert mapping.unique_cell_count == 2
    assert mapping.duplicate_cell_count == 0
    assert duplicate.unique_cell_count == 2
    assert duplicate.duplicate_cell_count == 1


def test_cell_center_quantization_rejects_off_lattice_coordinates():
    with pytest.raises(ValueError, match="not aligned"):
        quantize_cell_centers(
            np.array([[0.6, 0.5, 0.5]]), origin_m=np.zeros(3), dx_m=1.0
        )


def test_hexa_reconstruction_shares_vertices_and_keeps_pressure_fields():
    reference = 100.0
    origin = np.zeros(3)
    dx = 2.0e-7
    centers = origin + (np.array([[0, 0, 0], [1, 0, 0]]) + 0.5) * dx
    mapping = quantize_cell_centers(centers, origin_m=origin, dx_m=dx)
    grid = reconstruct_hexahedral_field(
        mapping=mapping,
        origin_m=origin,
        dx_m=dx,
        pressure_pa=np.array([reference + 1.0, reference + 2.0]),
        velocity_m_s=np.array([[1.0e-4, 0.0, 0.0], [2.0e-4, 0.0, 0.0]]),
        pressure_reference_pa=reference,
    )

    assert grid.n_cells == 2
    assert grid.n_points == 12
    assert np.unique(grid.celltypes).tolist() == [int(pv.CellType.HEXAHEDRON)]
    np.testing.assert_allclose(grid.cell_data["pressure_phy"], [101.0, 102.0])
    np.testing.assert_allclose(grid.cell_data["pressure_gauge_pa"], [1.0, 2.0])
    assert np.asarray(grid.cell_data["velocity_phy"]).shape == (2, 3)


def test_proteus_metadata_contract():
    metadata = build_proteus_metadata(
        Path("flow_field.vtu"), inlet_equivalent_diameter_m=2.0e-6
    )

    assert metadata["velocityField"] == "velocity_phy"
    assert metadata["pressureField"] == "pressure_phy"
    assert metadata["pressureGaugeField"] == "pressure_gauge_pa"
    assert metadata["inletEquivalentDiameterM"] == pytest.approx(2.0e-6)
