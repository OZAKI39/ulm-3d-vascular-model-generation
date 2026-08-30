from __future__ import annotations

import numpy as np

from utils.cfd_flow.repaired_topology_forensics import sparse_component_labels


def test_face_d3q19_and_full_connectivity_are_distinct() -> None:
    face_diagonal = np.asarray(((0, 0, 0), (1, 1, 0)), dtype=np.int64)
    _, sizes6 = sparse_component_labels(face_diagonal, 6, cells_per_axis=4)
    _, sizes18 = sparse_component_labels(face_diagonal, 18, cells_per_axis=4)
    assert sizes6.tolist() == [1, 1]
    assert sizes18.tolist() == [2]

    corner_diagonal = np.asarray(((0, 0, 0), (1, 1, 1)), dtype=np.int64)
    _, corner18 = sparse_component_labels(corner_diagonal, 18, cells_per_axis=4)
    _, corner26 = sparse_component_labels(corner_diagonal, 26, cells_per_axis=4)
    assert corner18.tolist() == [1, 1]
    assert corner26.tolist() == [2]
