from __future__ import annotations

from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

from utils.cfd_flow.restart_decode import (
    BINARY_SIZE_INVALID,
    D3Q19_DIRECTIONS,
    D3Q19_WEIGHTS,
    EXPECTED_CELL_COUNT,
    EXPECTED_DT_S,
    EXPECTED_DX_M,
    PRESSURE_REFERENCE_PA,
    REFERENCE_DENSITY_KG_M3,
    first_id_at_level,
    parse_d3q19_layout,
    parse_restart_header,
    read_restart_pdf,
    read_treelm_elemlist,
    reconstruct_macroscopic_field,
    reference_reproduction,
    restart_binary_size_contract,
    tree_ids_to_ijk,
    tree_levels,
    unique_cell_gate,
)
from utils.cfd_flow.steady_export import LatticeMapping, reconstruct_hexahedral_field


D3Q19_SOURCE_FIXTURE = """
integer, parameter :: d3q19_cxDir(3,19) &
  = reshape( [ -1,  0,  0, & ! W
                0, -1,  0, & ! S
                0,  0, -1, & ! B
                1,  0,  0, & ! E
                0,  1,  0, & ! N
                0,  0,  1, & ! T
                0, -1, -1, & ! BS
                0, -1,  1, & ! TS
                0,  1, -1, & ! BN
                0,  1,  1, & ! TN
               -1,  0, -1, & ! BW
                1,  0, -1, & ! BE
               -1,  0,  1, & ! TW
                1,  0,  1, & ! TE
               -1, -1,  0, & ! SW
               -1,  1,  0, & ! NW
                1, -1,  0, & ! SE
                1,  1,  0, & ! NE
                0,  0,  0 ], & ! C
             [ 3, 19 ] )
"""


def test_d3q19_source_contract_fixes_every_pdf_column() -> None:
    expected = np.asarray(
        [
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, -1],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [0, -1, -1],
            [0, -1, 1],
            [0, 1, -1],
            [0, 1, 1],
            [-1, 0, -1],
            [1, 0, -1],
            [-1, 0, 1],
            [1, 0, 1],
            [-1, -1, 0],
            [-1, 1, 0],
            [1, -1, 0],
            [1, 1, 0],
            [0, 0, 0],
        ],
        dtype=np.int8,
    )
    assert np.array_equal(parse_d3q19_layout(D3Q19_SOURCE_FIXTURE), expected)
    assert np.array_equal(D3Q19_DIRECTIONS, expected)


def test_current_pinned_d3q19_source_matches_decoder() -> None:
    source = Path(
        r"\\wsl.localhost\Ubuntu\home\lzy\apes-pinned\musubi_official\tem\source\tem_stencil_module.fpp"
    )
    if not source.is_file():
        pytest.skip("Pinned WSL source is not available in this environment")
    assert np.array_equal(parse_d3q19_layout(source.read_text(encoding="utf-8")), D3Q19_DIRECTIONS)


def test_restart_byte_size_and_float64_elementwise_reshape(tmp_path: Path) -> None:
    values = np.arange(3 * 19, dtype="<f8").reshape(3, 19)
    path = tmp_path / "restart.lsb"
    values.tofile(path)
    contract = restart_binary_size_contract(path, n_elems=3, n_components=19, n_dofs=1)
    assert contract["status"] == "PASS"
    decoded = read_restart_pdf(path, n_elems=3, n_components=19)
    assert decoded.dtype == np.dtype("<f8")
    assert decoded.shape == (3, 19)
    assert np.array_equal(decoded[2], values[2])


def test_restart_size_contract_rejects_one_missing_float(tmp_path: Path) -> None:
    path = tmp_path / "short.lsb"
    np.zeros(18, dtype="<f8").tofile(path)
    contract = restart_binary_size_contract(path, n_elems=1, n_components=19, n_dofs=1)
    assert contract["status"] == "FAIL"
    with pytest.raises(Exception, match=BINARY_SIZE_INVALID):
        read_restart_pdf(path, n_elems=1, n_components=19)


def test_synthetic_d3q19_density_and_velocity_recovery() -> None:
    pdf = D3Q19_WEIGHTS[None, :].copy()
    delta = 1.25e-4
    pdf[0, 3] += delta
    pdf[0, 0] -= delta
    field = reconstruct_macroscopic_field(pdf)
    expected_lattice_velocity_x = 2.0 * delta
    assert field.density_lattice[0] == pytest.approx(1.0)
    assert field.velocity_lattice[0] == pytest.approx([expected_lattice_velocity_x, 0.0, 0.0])
    assert field.velocity_phy[0, 0] == pytest.approx(
        expected_lattice_velocity_x * EXPECTED_DX_M / EXPECTED_DT_S
    )


def test_known_equilibrium_pdf_recovers_density_velocity_and_pressure() -> None:
    rho = 1.0025
    velocity = np.asarray([1.0e-5, -2.0e-5, 0.5e-5])
    cu = D3Q19_DIRECTIONS @ velocity
    equilibrium = D3Q19_WEIGHTS * rho * (
        1.0 + 3.0 * cu + 4.5 * cu**2 - 1.5 * float(velocity @ velocity)
    )
    field = reconstruct_macroscopic_field(equilibrium[None, :])
    assert field.density_lattice[0] == pytest.approx(rho, abs=1.0e-15)
    assert field.velocity_lattice[0] == pytest.approx(velocity, abs=1.0e-15)
    pressure_factor = REFERENCE_DENSITY_KG_M3 * EXPECTED_DX_M**2 / EXPECTED_DT_S**2
    assert field.pressure_phy[0] == pytest.approx(rho * pressure_factor / 3.0)


def test_unit_density_pressure_is_frozen_physical_reference() -> None:
    field = reconstruct_macroscopic_field(D3Q19_WEIGHTS[None, :])
    assert field.pressure_phy[0] == pytest.approx(PRESSURE_REFERENCE_PA, abs=1.0e-10)


@pytest.mark.parametrize(
    ("actual", "reference", "absolute_tolerance", "expected"),
    [(1.0 + 1.0e-8, 1.0, 1.0e-6, "PASS"), (1.0 + 1.0e-14, 1.0, 1.0e-12, "PASS"), (2.0, 1.0, 1.0e-12, "FAIL")],
)
def test_musubi_reference_validation_helper(
    actual: float,
    reference: float,
    absolute_tolerance: float,
    expected: str,
) -> None:
    assert (
        reference_reproduction(
            actual,
            reference,
            absolute_tolerance=absolute_tolerance,
        )["status"]
        == expected
    )


def test_treelm_elemlist_parser_is_interleaved_little_endian_int64(tmp_path: Path) -> None:
    records = np.asarray([[1, 8], [9, 256], [73, 264]], dtype="<i8")
    path = tmp_path / "elemlist.lsb"
    records.tofile(path)
    tree_ids, bits, contract = read_treelm_elemlist(path, n_elems=3)
    assert contract["status"] == "PASS"
    assert tree_ids.tolist() == [1, 9, 73]
    assert bits.tolist() == [8, 256, 264]


def test_known_treeid_to_coordinate_cases_follow_treelm_morton_order() -> None:
    level_one = np.asarray([1, 2, 3, 5, 8], dtype=np.int64)
    assert tree_levels(level_one).tolist() == [1, 1, 1, 1, 1]
    assert tree_ids_to_ijk(level_one).tolist() == [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 1],
    ]
    level = 9
    morton = 5 + 2 * 8
    tree_id = first_id_at_level(level) + morton
    assert tree_ids_to_ijk(np.asarray([tree_id])).tolist() == [[1, 2, 1]]


def test_exact_base_unique_cell_gate() -> None:
    indices = np.column_stack(
        (
            np.arange(EXPECTED_CELL_COUNT, dtype=np.int64),
            np.zeros(EXPECTED_CELL_COUNT, dtype=np.int64),
            np.zeros(EXPECTED_CELL_COUNT, dtype=np.int64),
        )
    )
    assert unique_cell_gate(indices)["status"] == "PASS"
    indices[-1] = indices[-2]
    failed = unique_cell_gate(indices)
    assert failed["status"] == "FAIL"
    assert failed["duplicate_cell_count"] == 1


def test_shared_vertex_hexa_reconstruction_and_cell_data() -> None:
    mapping = LatticeMapping(
        cell_indices=np.asarray([[0, 0, 0], [1, 0, 0]], dtype=np.int64),
        maximum_alignment_error_m=0.0,
        duplicate_cell_count=0,
        unique_cell_count=2,
    )
    velocity = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    pressure = np.asarray([10.0, 11.0])
    grid = reconstruct_hexahedral_field(
        mapping=mapping,
        origin_m=np.zeros(3),
        dx_m=1.0,
        pressure_pa=pressure,
        velocity_m_s=velocity,
        pressure_reference_pa=9.0,
    )
    assert grid.n_cells == 2
    assert grid.n_points == 12
    assert np.all(grid.celltypes == int(pv.CellType.HEXAHEDRON))
    assert np.array_equal(grid.cell_data["velocity_phy"], velocity)
    assert np.array_equal(grid.cell_data["pressure_phy"], pressure)
    assert np.array_equal(grid.cell_data["pressure_gauge_pa"], [1.0, 2.0])


def test_restart_header_contract_parser(tmp_path: Path) -> None:
    header = tmp_path / "lastHeader.lua"
    header.write_text(
        """
binary_name = { '/mnt/e/data/restart.lsb' }
solver_configFile = '/mnt/e/data/musubi.lua'
time_point = { iter = 198064 }
nElems = 221109
nDofs = 1
solver = 'Musubi_v2.0.0-4-g4e8b27'
varsys = { variable = { { name = 'pdf', ncomponents = 19 } }, nScalars = 19 }
""",
        encoding="utf-8",
    )
    parsed = parse_restart_header(header)
    assert parsed.iteration == 198064
    assert parsed.n_elems == 221109
    assert parsed.binary_path == Path("E:/data/restart.lsb")
    assert parsed.variable_name == "pdf"
    assert parsed.n_components == 19
