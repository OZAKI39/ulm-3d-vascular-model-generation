from pathlib import Path

import numpy as np

from utils.cfd_flow.periodic_pipe_force import (
    CASES,
    _read_ascii_result,
    generate_periodic_musubi_lua,
    generate_periodic_seeder_lua,
)


def test_periodic_contract_has_no_pressure_boundaries() -> None:
    case = CASES["axis_n27"]
    seeder = generate_periodic_seeder_lua(case)
    musubi = generate_periodic_musubi_lua(case)
    assert "kind = 'periodic'" in seeder
    assert "calc_dist = true" in seeder
    assert "kind = 'wall_libb'" in musubi
    assert "kind = 'fluid_incompressible'" in musubi
    assert "glob_source" in musubi
    assert "adaptive_flux" not in musubi
    assert "pressure_eq" not in musubi
    assert "maximum_iterations = 2000" in musubi


def test_oblique_translation_is_lattice_exact_and_axial() -> None:
    case = CASES["oblique_n27"]
    translation = case.translation_m / case.dx_m
    assert np.allclose(translation, np.rint(translation))
    assert np.allclose(case.translation_m / case.length_m, case.direction)


def test_ascii_result_reader(tmp_path: Path) -> None:
    result = tmp_path / "tracking.res"
    result.write_text("# metadata\n# time value_01\n0 1\n2 3\n", encoding="utf-8")
    header, values = _read_ascii_result(result)
    assert header == ["time", "value_01"]
    assert values.tolist() == [[0.0, 1.0], [2.0, 3.0]]
