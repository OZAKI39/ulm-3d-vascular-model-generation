from __future__ import annotations

from pathlib import Path

import numpy as np

from utils.cfd_flow.exact_link_flux import BoundaryReconstruction, INVERSE_DIRECTIONS
from utils.cfd_flow.full_timestep_mass_referee import stable_delta
from utils.cfd_flow.musubi_boundary_mass_referee import (
    MeshContract,
    pull_fetch_pdfs_runtime,
)


def test_boundary_neighbor_id_overrides_coordinate_only_pull() -> None:
    coordinates = np.asarray(((0, 0, 0), (1, 0, 0)), dtype=np.int64)
    incoming = np.zeros((1, 18), dtype=bool)
    incoming[0, 3] = True
    boundary = BoundaryReconstruction(
        label="wall",
        boundary_id=1,
        property_rows=np.asarray((0,), dtype=np.int64),
        cell_indices=np.asarray((1,), dtype=np.int64),
        outward_masks=np.zeros((1, 18), dtype=bool),
        incoming_masks=incoming,
        raw_normals=np.zeros((1, 3), dtype=np.int64),
        normal_indices=np.asarray((3,), dtype=np.int64),
    )
    mesh = MeshContract(
        mesh=Path("."),
        tree_ids=np.arange(2),
        cell_ijk=coordinates,
        lookup={tuple(row): index for index, row in enumerate(coordinates)},
        boundaries={"wall": boundary},
        boundary_labels=("wall",),
        qvalues_by_cell=np.ones((2, 18)),
    )
    state = np.arange(38, dtype=np.float64).reshape(2, 19)

    fetched = pull_fetch_pdfs_runtime(state, mesh, np.asarray((1,)), set())

    assert fetched[0, 3] == state[1, INVERSE_DIRECTIONS[3]]
    assert fetched[0, 3] != state[0, 3]
    assert stable_delta(fetched, state[1:2]) == np.sum(
        fetched - state[1:2], dtype=np.float64
    )
