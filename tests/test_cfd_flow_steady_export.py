"""Targeted tests for the one-shot frozen steady-field export."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

from utils.cfd_flow.io import FlowError
from utils.cfd_flow.steady_export import (
    PRESSURE_REFERENCE_PA,
    AsciiSpatialField,
    ascii_path_preflight,
    build_proteus_metadata,
    generate_ascii_spatial_harvester_lua,
    harvester_lua_contract,
    parse_ascii_spatial_schema,
    parse_uniform_mesh_lattice,
    quantize_cell_centers,
    read_ascii_spatial_field,
    reconstruct_hexahedral_field,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _write_ascii_spatial(tmp_path: Path) -> tuple[Path, Path]:
    companion = tmp_path / "f.lua"
    companion.write_text(
        """format = 'asciispatial'
nElems = 2
variable = {
  { name = 'pressure_phy', ncomponents = 1 },
  { name = 'velocity_phy', ncomponents = 3 }
}
""",
        encoding="utf-8",
    )
    result = tmp_path / "f_p00000_t1.000E+00.res"
    result.write_text(
        """# Rank of the process:       0
# velocity_phy_02 coordZ pressure_phy coordX velocity_phy_03 coordY velocity_phy_01
2.0 30.5 101.0 10.5 3.0 20.5 1.0
5.0 30.5 102.0 11.5 6.0 20.5 4.0
""",
        encoding="utf-8",
    )
    return result, companion


def test_ascii_spatial_schema_and_vector_parsing_follow_actual_header(tmp_path):
    result, companion = _write_ascii_spatial(tmp_path)
    schema = parse_ascii_spatial_schema(result, companion)
    field = read_ascii_spatial_field(result, companion)

    assert schema["status"] == "PASS"
    np.testing.assert_allclose(field.coordinates_m, [[10.5, 20.5, 30.5], [11.5, 20.5, 30.5]])
    assert field.pressure_pa == pytest.approx([101.0, 102.0])
    np.testing.assert_allclose(field.velocity_m_s, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])


def test_ascii_spatial_ambiguous_schema_stops(tmp_path):
    result, companion = _write_ascii_spatial(tmp_path)
    result.write_text(
        "# coordX coordY coordZ pressure_phy velocity_phy_01 velocity_phy_02\n"
        "0 0 0 1 2 3\n",
        encoding="utf-8",
    )
    with pytest.raises(FlowError, match="CFD_FLOW_STEADY_ASCII_SCHEMA_INVALID"):
        parse_ascii_spatial_schema(result, companion)


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


def test_uniform_mesh_header_parses_nested_origin_table(tmp_path):
    header = tmp_path / "header.lua"
    header.write_text(
        """boundingbox = {
  origin = { 1.0E-06, 2.0E-06, 3.0E-06 },
  length = 102.4E-06
}
nElems = 221109
minLevel = 9
maxLevel = 9
""",
        encoding="utf-8",
    )
    result = parse_uniform_mesh_lattice(header)
    assert result["origin_m"] == pytest.approx([1.0e-6, 2.0e-6, 3.0e-6])
    assert result["dx_m"] == pytest.approx(2.0e-7)
    assert result["n_elems"] == 221109


def test_hexa_reconstruction_shares_vertices_and_assigns_cell_data():
    origin = np.zeros(3)
    dx = 2.0e-7
    centers = origin + (np.array([[0, 0, 0], [1, 0, 0]]) + 0.5) * dx
    mapping = quantize_cell_centers(centers, origin_m=origin, dx_m=dx)
    grid = reconstruct_hexahedral_field(
        mapping=mapping,
        origin_m=origin,
        dx_m=dx,
        pressure_pa=np.array([PRESSURE_REFERENCE_PA + 1.0, PRESSURE_REFERENCE_PA + 2.0]),
        velocity_m_s=np.array([[1.0e-4, 0.0, 0.0], [2.0e-4, 0.0, 0.0]]),
    )

    assert grid.n_cells == 2
    assert grid.n_points == 12
    assert np.unique(grid.celltypes).tolist() == [int(pv.CellType.HEXAHEDRON)]
    assert np.asarray(grid.cell_data["velocity_phy"]).shape == (2, 3)
    assert np.asarray(grid.cell_data["pressure_phy"]).shape == (2,)
    assert grid.cell_data["pressure_gauge_pa"] == pytest.approx([1.0, 2.0])


def test_pressure_gauge_conversion_does_not_modify_pressure_phy():
    field = AsciiSpatialField(
        coordinates_m=np.array([[1.0e-7, 1.0e-7, 1.0e-7]]),
        pressure_pa=np.array([PRESSURE_REFERENCE_PA + 7.5]),
        velocity_m_s=np.zeros((1, 3)),
        columns=(
            "coordX",
            "coordY",
            "coordZ",
            "pressure_phy",
            "velocity_phy_01",
            "velocity_phy_02",
            "velocity_phy_03",
        ),
        result_path=Path("f.res"),
        companion_path=Path("f.lua"),
    )
    mapping = quantize_cell_centers(field.coordinates_m, origin_m=np.zeros(3), dx_m=2.0e-7)
    grid = reconstruct_hexahedral_field(
        mapping=mapping,
        origin_m=np.zeros(3),
        dx_m=2.0e-7,
        pressure_pa=field.pressure_pa,
        velocity_m_s=field.velocity_m_s,
    )

    assert grid.cell_data["pressure_phy"][0] == pytest.approx(PRESSURE_REFERENCE_PA + 7.5)
    assert grid.cell_data["pressure_gauge_pa"][0] == pytest.approx(7.5)


def test_proteus_metadata_contract_uses_pressure_phy_and_null_normal(tmp_path):
    metadata = build_proteus_metadata(
        flow_vtu=tmp_path / "flow_field.vtu",
        inlet_area_m2=4.0e-12,
        dx_m=2.0e-7,
    )

    assert metadata["lengthUnit"] == 1.0
    assert metadata["velocityUnit"] == 1.0
    assert metadata["velocityField"] == "velocity_phy"
    assert metadata["pressureField"] == "pressure_phy"
    assert metadata["coordinateUnit"] == "m"
    assert metadata["velocityUnitName"] == "m/s"
    assert metadata["inletNormal"] is None
    assert metadata["inlet_normal_policy"] == "AUTO_DETECT_BY_BACKPROPAGATION_LATER"


def test_harvester_lua_is_ascii_only_and_uses_exact_solver_config():
    solver = "/mnt/e/project/source/diagnostic_musubi.lua"
    restart = "/mnt/e/project/source/restart/steady_lastHeader.lua"
    text = generate_ascii_spatial_harvester_lua(
        solver_config_wsl=solver,
        restart_header_wsl=restart,
    )
    contract = harvester_lua_contract(
        text,
        solver_config_wsl=solver,
        restart_header_wsl=restart,
    )

    assert contract["status"] == "PASS"
    assert "require 'diagnostic_musubi'" in text
    assert "variable = { 'pressure_phy', 'velocity_phy' }" in text
    assert "format = 'asciiSpatial'" in text
    assert "format = 'vtk'" not in text


def test_ascii_path_preflight_stays_below_project_limit():
    result = ascii_path_preflight("roi003274_steady_lbm", "4.836E-03")
    assert result["status"] == "PASS"
    assert result["maximum_predicted_filename_length"] == 55
    assert result["predicted_paths"]["actual"].endswith(
        "roi003274_steady_lbm_f_p00000_t4.836E-03.res"
    )


def test_formal_entry_contains_one_harvester_and_no_solver_launch():
    source = (PROJECT_ROOT / "utils" / "cfd_flow" / "steady_export.py").read_text(
        encoding="utf-8"
    )
    assert source.count("harvest_run = run_wsl_tool(") == 1
    assert 'environment.binaries["musubi"]' not in source
    assert 'environment.binaries["seeder"]' not in source
    assert "output = { format = 'vtk' }" not in source
