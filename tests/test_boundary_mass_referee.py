from __future__ import annotations

from pathlib import Path

import numpy as np

from utils.cfd_flow.exact_link_flux import INVERSE_DIRECTIONS
from utils.cfd_flow.musubi_boundary_mass_referee import (
    _wall_operations,
    boundary_window_closure,
    conservation_identity_residual,
    detect_slot_overlaps,
    earliest_gate_pass,
    fetch_storage_slot,
    kg_s_per_lattice_population,
    lattice_delta_to_kg_s,
    pull_fetch_pdfs_runtime,
    replay_sequential_writes,
    significant_time_averaged_backflow,
    trapezoidal_integral,
)
from utils.cfd_flow.restart_decode import D3Q19_DIRECTIONS


def test_pull_fetch_memory_slot_and_inverse_mapping() -> None:
    coordinates = np.asarray(((0, 0, 0), (1, 0, 0)), dtype=np.int64)
    lookup = {tuple(row): index for index, row in enumerate(coordinates)}
    # direction -x at target x=1 fetches source x=2, which is missing.
    cell, direction, kind = fetch_storage_slot(
        target_cell=1, direction=0, cell_ijk=coordinates, coordinate_lookup=lookup
    )
    assert (cell, direction) == (1, int(INVERSE_DIRECTIONS[0]))
    assert kind == "MISSING_NEIGHBOR_LOCAL_INVERSE"
    # direction +x at target x=1 pulls from x=0 in the same PDF direction.
    cell, direction, kind = fetch_storage_slot(
        target_cell=1, direction=3, cell_ijk=coordinates, coordinate_lookup=lookup
    )
    assert (cell, direction, kind) == (0, 3, "PULL_SOURCE_SAME_DIRECTION")
    assert np.array_equal(D3Q19_DIRECTIONS[INVERSE_DIRECTIONS], -D3Q19_DIRECTIONS)


def test_runtime_solid_pull_maps_current_and_source_to_local_inverse() -> None:
    from utils.cfd_flow.musubi_boundary_mass_referee import MeshContract

    coordinates = np.asarray(((0, 0, 0), (1, 0, 0), (2, 0, 0)), dtype=np.int64)
    values = np.arange(3 * 19, dtype=np.float64).reshape(3, 19)
    mesh = MeshContract(
        mesh=Path("."),
        tree_ids=np.arange(3),
        cell_ijk=coordinates,
        lookup={tuple(row): index for index, row in enumerate(coordinates)},
        boundaries={},
        boundary_labels=(),
        qvalues_by_cell=np.empty((3, 18)),
    )
    fetched = pull_fetch_pdfs_runtime(values, mesh, np.asarray((1, 2)), {0, 2})
    # target 1, +x pulls from solid source 0 and therefore uses local inverse.
    assert fetched[0, 3] == values[1, INVERSE_DIRECTIONS[3]]
    # solid target 2 uses local inverse even where a non-solid source exists.
    assert fetched[1, 3] == values[2, INVERSE_DIRECTIONS[3]]


def test_wall_libb_neighbor_uses_runtime_solid_pull_connectivity() -> None:
    from utils.cfd_flow.exact_link_flux import BoundaryReconstruction
    from utils.cfd_flow.musubi_boundary_mass_referee import MeshContract

    coordinates = np.asarray(((0, 0, 0), (1, 0, 0)), dtype=np.int64)
    values = np.arange(2 * 19, dtype=np.float64).reshape(2, 19) + 1.0
    incoming = np.zeros((1, 18), dtype=bool)
    incoming[0, 3] = True
    wall = BoundaryReconstruction(
        label="wall",
        boundary_id=1,
        property_rows=np.asarray((0,), dtype=np.int64),
        cell_indices=np.asarray((0,), dtype=np.int64),
        outward_masks=np.zeros((1, 18), dtype=bool),
        incoming_masks=incoming,
        raw_normals=np.zeros((1, 3), dtype=np.int64),
        normal_indices=np.asarray((3,), dtype=np.int64),
    )
    qvalues = np.ones((2, 18), dtype=np.float64)
    inverse = int(INVERSE_DIRECTIONS[3])
    qvalues[0, inverse] = 0.25
    mesh = MeshContract(
        mesh=Path("."),
        tree_ids=np.arange(2),
        cell_ijk=coordinates,
        lookup={tuple(row): index for index, row in enumerate(coordinates)},
        boundaries={"wall": wall},
        boundary_labels=("wall",),
        qvalues_by_cell=qvalues,
    )

    coordinate_only = _wall_operations(values, mesh)[0][3]
    runtime_solid = _wall_operations(values, mesh, {1})[0][3]
    expected = 0.5 * values[0, inverse] + 0.5 * values[0, 3]

    assert coordinate_only != runtime_solid
    assert runtime_solid == expected


def test_runtime_order_overlap_and_sequential_last_writer() -> None:
    overlaps = detect_slot_overlaps({"wall": {(1, 2)}, "inlet": {(1, 2), (3, 4)}})
    assert overlaps[("wall", "inlet")] == {(1, 2)}
    state = np.zeros((5, 19))
    replay, deltas = replay_sequential_writes(
        state, (("wall", 1, 2, 0.25), ("inlet", 1, 2, 0.75))
    )
    assert replay[1, 2] == 0.75
    assert deltas == {"wall": 0.25, "inlet": 0.5}


def test_wall_inlet_and_each_outlet_delta_bookkeeping() -> None:
    state = np.zeros((6, 19))
    operations = (
        ("wall", 0, 1, 0.0),
        ("inlet", 1, 2, 1.0),
        ("outlet_01", 2, 3, -0.1),
        ("outlet_02", 3, 4, -0.2),
        ("outlet_03", 4, 5, -0.3),
    )
    _, deltas = replay_sequential_writes(state, operations)
    assert deltas["wall"] == 0.0
    assert deltas["inlet"] == 1.0
    assert deltas["outlet_01"] == -0.1
    assert deltas["outlet_02"] == -0.2
    assert deltas["outlet_03"] == -0.3


def test_signed_conversion_trapezoid_and_discrete_identity() -> None:
    factor = kg_s_per_lattice_population()
    assert lattice_delta_to_kg_s(-2.0) == -2.0 * factor
    assert trapezoidal_integral((0, 10), (2.0, 4.0)) == 30.0
    assert conservation_identity_residual(1.0, 1.0, 2.0) == 0.0
    assert boundary_window_closure(3.0, 2.0, 10.0) == 0.1


def test_time_averaged_backflow_and_earliest_gate_pass() -> None:
    assert significant_time_averaged_backflow(10.0, (1.0, -0.6, 2.0))
    assert not significant_time_averaged_backflow(10.0, (1.0, -0.4, 2.0))
    base = {
        "R_mass_short": 0.001,
        "R_mass_long": 0.001,
        "R_velocity": 0.001,
        "R_pressure": 0.001,
        "R_inlet": 1e-12,
        "R_conservation_identity": 1e-12,
        "boundary_window_closure": 1e-4,
        "significant_time_averaged_backflow": False,
        "minimum_pdf": 0.01,
        "maximum_lattice_speed": 0.001,
        "inlet_globbc": 287,
    }
    failed = dict(base, iteration=100, R_mass_short=0.02)
    passed = dict(base, iteration=200)
    assert earliest_gate_pass((passed, failed))["iteration"] == 200
