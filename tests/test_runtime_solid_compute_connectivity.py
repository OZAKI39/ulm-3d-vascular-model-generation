from __future__ import annotations

from pathlib import Path

import numpy as np

from utils.cfd_flow.exact_link_flux import INVERSE_DIRECTIONS
from utils.cfd_flow.full_timestep_mass_referee import (
    pull_compute_targets,
    pull_link_counts,
)
from utils.cfd_flow.musubi_boundary_mass_referee import MeshContract


def test_runtime_solid_cells_remain_targets_with_inverse_connectivity() -> None:
    coordinates = np.asarray(((0, 0, 0), (1, 0, 0), (2, 0, 0)), dtype=np.int64)
    mesh = MeshContract(
        mesh=Path("."),
        tree_ids=np.arange(3),
        cell_ijk=coordinates,
        lookup={tuple(row): index for index, row in enumerate(coordinates)},
        boundaries={},
        boundary_labels=(),
        qvalues_by_cell=np.empty((3, 18)),
    )
    state = np.arange(57, dtype=np.float64).reshape(3, 19)

    fetched = pull_compute_targets(state, mesh, {0})
    counts = pull_link_counts(mesh, {0})

    assert fetched.shape == state.shape
    assert fetched[0, 3] == state[0, INVERSE_DIRECTIONS[3]]
    assert fetched[1, 3] == state[1, INVERSE_DIRECTIONS[3]]
    assert counts["current_solid_inverse"] == 19
    assert counts["source_solid_inverse"] == 1
    assert counts["total"] == state.size
