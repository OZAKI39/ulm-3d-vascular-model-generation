from __future__ import annotations

from pathlib import Path

import numpy as np

from utils.cfd_flow.musubi_boundary_mass_referee import MeshContract
from utils.cfd_flow.repaired_topology_forensics import compare_scaled_meshes


def _contract(tree_ids: tuple[int, ...], qvalues: np.ndarray) -> MeshContract:
    count = len(tree_ids)
    return MeshContract(
        mesh=Path("mesh"),
        tree_ids=np.asarray(tree_ids, dtype=np.int64),
        cell_ijk=np.zeros((count, 3), dtype=np.int64),
        lookup={},
        boundaries={},
        boundary_labels=(),
        qvalues_by_cell=np.asarray(qvalues, dtype=np.float64),
    )


def test_q_comparison_aligns_rows_by_tree_id_and_reports_dimensionless_error() -> None:
    patched = _contract((3, 1), np.asarray(((0.75, np.nan), (0.25, np.nan))))
    scaled = _contract((1, 3), np.asarray(((0.25, np.nan), (0.750001, np.nan))))
    patched_bnd = np.asarray(((2, 0), (1, 0)), dtype=np.int64)
    scaled_bnd = np.asarray(((1, 0), (2, 0)), dtype=np.int64)
    result = compare_scaled_meshes(patched, scaled, patched_bnd, scaled_bnd)
    assert result["tree_id_set_exact_match"] is True
    assert result["tree_id_order_exact_match"] is False
    assert result["boundary_id_exact_match"] is True
    assert result["common_q_links"] == 2
    assert np.isclose(result["q_patched_minus_scaled"]["max"], 1.0e-6)
