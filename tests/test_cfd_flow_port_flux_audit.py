from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from utils.cfd_flow.port_flux_audit import (
    BACKFLOW_CONFIRMED,
    EXPECTED_SIDE_NAMES,
    BoundaryHeader,
    absolute_magnitude_diagnostic,
    boundary_binary_contract,
    boundary_cell_fluxes,
    build_face_neighbor_graph,
    build_port_topology,
    classify_audit,
    extract_boundary_property_indices,
    find_stable_window,
    lattice_internal_cutset_sweep,
    parse_bnd_header,
    parse_boundary_property_header,
    parse_treelm_side_contract,
    read_boundary_ids,
    signed_balance,
)


EXPECTED_OFFSETS = np.asarray(
    [
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (0, -1, -1),
        (0, -1, 1),
        (0, 1, -1),
        (0, 1, 1),
        (-1, 0, -1),
        (1, 0, -1),
        (-1, 0, 1),
        (1, 0, 1),
        (-1, -1, 0),
        (-1, 1, 0),
        (1, -1, 0),
        (1, 1, 0),
        (-1, -1, -1),
        (1, -1, -1),
        (-1, 1, -1),
        (1, 1, -1),
        (-1, -1, 1),
        (1, -1, 1),
        (-1, 1, 1),
        (1, 1, 1),
    ],
    dtype=np.int8,
)


def _side_source_fixture() -> str:
    names = ", ".join(repr(f"{name:>3}") for name in EXPECTED_SIDE_NAMES)
    values = [*EXPECTED_OFFSETS[:, 0], *EXPECTED_OFFSETS[:, 1], *EXPECTED_OFFSETS[:, 2]]
    return (
        f"character(len=3), parameter :: qDirName(qQQQ) = [ {names} ]\n"
        f"integer, dimension(qQQQ,3), parameter :: qOffset = "
        f"reshape((/{', '.join(str(int(value)) for value in values)}/),(/qQQQ,3/))\n"
    )


def test_signed_mass_balance_preserves_negative_outlet() -> None:
    result = signed_balance(10.0, (6.0, -1.0, 5.0))
    assert result["q_out_sum_m3_s"] == 10.0
    assert result["relative_error"] == 0.0


def test_absolute_outlet_formula_is_not_conservation() -> None:
    signed = signed_balance(10.0, (6.0, -1.0, 5.0))["relative_error"]
    magnitude = absolute_magnitude_diagnostic(10.0, (6.0, -1.0, 5.0))
    assert signed == 0.0
    assert magnitude == pytest.approx(0.2)


def test_boundary_property_header_and_bit_extraction() -> None:
    header = parse_boundary_property_header(
        """
property = {
  { label = 'has boundaries', bitpos = 3, nElems = 3 },
  { label = 'has qVal', bitpos = 8, nElems = 1 }
}
"""
    )
    assert header.bit_position == 3
    assert header.element_count == 3
    bits = np.asarray([0, 8, 256, 264, 0], dtype=np.int64)
    assert extract_boundary_property_indices(bits, 3).tolist() == [1, 3]


def test_bnd_binary_size_and_element_major_layout(tmp_path: Path) -> None:
    values = np.arange(3 * 26, dtype="<i8").reshape(3, 26)
    path = tmp_path / "bnd.lsb"
    values.tofile(path)
    contract = boundary_binary_contract(path, element_count=3, side_count=26)
    assert contract["status"] == "PASS"
    assert contract["actual_bytes"] == 3 * 26 * 8
    parsed = read_boundary_ids(path, element_count=3, side_count=26)
    assert np.array_equal(parsed, values)


def test_bnd_label_order_is_parsed_from_lua() -> None:
    header = parse_bnd_header(
        """
nSides = 26
nBCtypes = 5
bclabel = { 'wall', 'outlet_02', 'outlet_01', 'outlet_03', 'inlet' }
"""
    )
    assert header.side_count == 26
    assert header.labels == ("wall", "outlet_02", "outlet_01", "outlet_03", "inlet")


def test_26_side_order_and_qoffset_are_source_parsed() -> None:
    names, offsets = parse_treelm_side_contract(_side_source_fixture())
    assert names == EXPECTED_SIDE_NAMES
    assert np.array_equal(offsets, EXPECTED_OFFSETS)


def test_direct_six_face_normals_are_exact() -> None:
    _, offsets = parse_treelm_side_contract(_side_source_fixture())
    assert offsets[:6].tolist() == [
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]


@dataclass
class _Patch:
    area_um2: float
    outward_normal: np.ndarray


class _Partition:
    def __init__(self, patches: dict[str, _Patch]) -> None:
        self._patches = patches

    def patch(self, label: str) -> _Patch:
        return self._patches[label]


def _four_port_topology() -> tuple[dict[str, object], dict[str, object]]:
    boundary_ids = np.zeros((4, 26), dtype=np.int64)
    boundary_ids[0, 0] = 5  # inlet W
    boundary_ids[1, 3] = 3  # outlet_01 E
    boundary_ids[2, 4] = 2  # outlet_02 N
    boundary_ids[3, 5] = 4  # outlet_03 T
    partition = _Partition(
        {
            "inlet": _Patch(1.0e12, np.asarray([-1.0, 0.0, 0.0])),
            "outlet_01": _Patch(1.0e12, np.asarray([1.0, 0.0, 0.0])),
            "outlet_02": _Patch(1.0e12, np.asarray([0.0, 1.0, 0.0])),
            "outlet_03": _Patch(1.0e12, np.asarray([0.0, 0.0, 1.0])),
        }
    )
    header = BoundaryHeader(26, 5, ("wall", "outlet_02", "outlet_01", "outlet_03", "inlet"))
    return build_port_topology(
        boundary_ids=boundary_ids,
        property_element_indices=np.arange(4),
        boundary_header=header,
        side_names=EXPECTED_SIDE_NAMES,
        side_offsets=EXPECTED_OFFSETS,
        partition=partition,  # type: ignore[arg-type]
        dx_m=1.0,
    )


def test_boundary_face_area_vectors_use_oriented_direct_faces() -> None:
    qc, topology = _four_port_topology()
    assert qc["status"] == "PASS"
    assert topology["inlet"].qc["lattice_area_vector_m2"] == [-1.0, 0.0, 0.0]
    assert topology["outlet_02"].qc["lattice_area_vector_m2"] == [0.0, 1.0, 0.0]
    assert all(item.qc["dot_product"] > 0.0 for item in topology.values())


def test_boundary_cell_estimator_uses_outward_sign_convention() -> None:
    _, topology = _four_port_topology()
    velocity = np.asarray(
        [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    result = boundary_cell_fluxes(topology, velocity, np.ones(4), dx_m=1.0, rho0_kg_m3=2.0)
    assert result["ports"]["inlet"]["q_m3_s"] == 1.0
    assert result["ports"]["inlet"]["mass_flow_kg_s"] == 2.0
    assert result["ports"]["outlet_03"]["q_m3_s"] == 3.0


def test_six_neighbor_graph_contains_only_face_neighbors() -> None:
    coordinates = np.asarray([[0, 0, 0], [1, 0, 0], [1, 1, 0]], dtype=np.int64)
    graph = build_face_neighbor_graph(coordinates)
    assert graph[0].tolist() == [-1, -1, -1, 1, -1, -1]
    assert graph[1].tolist() == [0, -1, -1, -1, 2, -1]


def _pipe_sweep(
    velocity_x: np.ndarray,
    density: np.ndarray,
    *,
    maximum_depth: int,
) -> tuple[list[dict[str, object]], int | None]:
    coordinates = np.column_stack(
        (np.arange(len(velocity_x), dtype=np.int64), np.zeros((len(velocity_x), 2), dtype=np.int64))
    )
    velocity = np.column_stack((velocity_x, np.zeros((len(velocity_x), 2))))
    graph = build_face_neighbor_graph(coordinates)
    return lattice_internal_cutset_sweep(
        graph=graph,
        cell_ijk=coordinates,
        velocity_m_s=velocity,
        density_lattice=density,
        seeds_by_port={"inlet": np.asarray([0]), "outlet_01": np.asarray([len(velocity_x) - 1])},
        dx_m=1.0,
        rho0_kg_m3=1.0,
        maximum_depth=maximum_depth,
    )


def test_uniform_cartesian_pipe_cutset_flux_is_depth_invariant() -> None:
    rows, overlap = _pipe_sweep(np.ones(10), np.ones(10), maximum_depth=4)
    assert overlap is None
    for row in rows:
        assert row["Q_m3_s"] == pytest.approx(1.0)
        assert row["Mdot_kg_s"] == pytest.approx(1.0)


def test_cut_face_orientation_points_from_core_toward_port() -> None:
    rows, _ = _pipe_sweep(np.ones(6), np.ones(6), maximum_depth=1)
    by_port = {row["port"]: row for row in rows}
    assert by_port["inlet"]["q_outward_m3_s"] == pytest.approx(-1.0)
    assert by_port["inlet"]["Q_m3_s"] == pytest.approx(1.0)
    assert by_port["outlet_01"]["q_outward_m3_s"] == pytest.approx(1.0)


def test_cutset_volume_flux_uses_arithmetic_face_velocity() -> None:
    rows, _ = _pipe_sweep(np.asarray([0.0, 2.0, 2.0]), np.ones(3), maximum_depth=1)
    inlet = next(row for row in rows if row["port"] == "inlet")
    assert inlet["Q_m3_s"] == pytest.approx(1.0)


def test_cutset_mass_flux_averages_rho_u_vectors() -> None:
    rows, _ = _pipe_sweep(
        np.asarray([0.0, 2.0, 2.0]),
        np.asarray([1.0, 2.0, 2.0]),
        maximum_depth=1,
    )
    inlet = next(row for row in rows if row["port"] == "inlet")
    assert inlet["Mdot_kg_s"] == pytest.approx(2.0)


def test_cross_port_overlap_depth_is_detected_and_persists_invalid() -> None:
    rows, first_overlap = _pipe_sweep(np.ones(4), np.ones(4), maximum_depth=3)
    assert first_overlap == 3
    depth_three = [row for row in rows if row["depth_cells"] == 3]
    assert all(row["cross_port_overlap"] for row in depth_three)
    assert all(not row["global_balance_valid"] for row in depth_three)


def _stable_depth(depth: int, outlet_two: float = -1.0) -> dict[str, object]:
    outlet_one = 2.0
    outlet_three = 9.0
    inlet = outlet_one + outlet_two + outlet_three
    return {
        "depth_cells": depth,
        "global_balance_valid": True,
        "cross_port_overlap": False,
        "q_in_m3_s": 10.0,
        "q_out_m3_s": {
            "outlet_01": 2.0,
            "outlet_02": -1.0,
            "outlet_03": 9.0,
        },
        "mass_in_kg_s": inlet,
        "mass_out_kg_s": {
            "outlet_01": outlet_one,
            "outlet_02": outlet_two,
            "outlet_03": outlet_three,
        },
        "inlet_mass_flow_relative_error": 0.0,
        "signed_mass_balance_error": 0.0,
    }


def test_stable_window_requires_three_consecutive_depths() -> None:
    assert not find_stable_window([_stable_depth(1), _stable_depth(2)])["found"]
    result = find_stable_window([_stable_depth(1), _stable_depth(2), _stable_depth(3)])
    assert result["found"]
    assert result["depths"] == [1, 2, 3]


def test_backflow_confirmation_rule_requires_stable_negative_outlet_two() -> None:
    stable = find_stable_window([_stable_depth(1), _stable_depth(2), _stable_depth(3)])
    classification = classify_audit(
        topology_pass=True,
        stable_window=stable,
        legacy_outlet_02_m3_s=-2.0,
    )
    assert classification["status"] == BACKFLOW_CONFIRMED
    assert classification["outlet_02_backflow_confirmed"] == "YES"


def test_synthetic_reverse_flow_keeps_negative_outlet_in_balance() -> None:
    result = signed_balance(10.0, (2.0, -1.0, 9.0))
    assert result["relative_error"] == 0.0
