from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from utils.cfd_flow.exact_link_flux import BoundaryReconstruction, INVERSE_DIRECTIONS
from utils.cfd_flow.musubi_boundary_mass_referee import (
    MeshContract,
    _wall_operations,
)
from utils.cfd_flow.qvalue_repair_validation import (
    libb_reference_replacement,
    q_error_metrics,
    tiny_cylinder_gate,
)


def test_q_error_metrics_report_median_p95_and_bias() -> None:
    metrics = q_error_metrics(np.asarray((-0.02, -0.01, 0.0, 0.01, 0.02)))

    assert metrics["count"] == 5
    assert metrics["bias"] == pytest.approx(0.0, abs=1.0e-18)
    assert metrics["median_absolute"] == 0.01
    assert metrics["p95_absolute"] == 0.02


def test_repaired_tiny_cylinder_passes_analytic_oracle() -> None:
    root = Path(__file__).resolve().parents[1]
    case = (
        root
        / "outputs/cfd_flow/healthy_mouse_capillary_port_grid_sensitivity_"
        "research_anchor003274_20260830/qvalue_repair/"
        "tiny_cylinder_repaired_n16"
    )

    result = tiny_cylinder_gate(case / "mesh", case / "source_mesh_summary.json")

    assert result["status"] == "PASS"
    assert result["q_error"]["median_absolute"] <= 0.01
    assert result["q_error"]["p95_absolute"] <= 0.02
    assert abs(result["q_error"]["bias"]) <= 0.01
    assert result["audit"]["uniform_halfway_fallback_fraction"] == 0.0
    assert result["audit"]["direction_mismatch_count"] == 0


def _single_wall_mesh(qvalue: float) -> tuple[np.ndarray, MeshContract, int]:
    active = 0
    inverse = int(INVERSE_DIRECTIONS[active])
    incoming = np.zeros((1, 18), dtype=bool)
    incoming[0, active] = True
    boundary = BoundaryReconstruction(
        label="wall",
        boundary_id=1,
        property_rows=np.asarray((0,), dtype=np.int64),
        cell_indices=np.asarray((0,), dtype=np.int64),
        outward_masks=np.zeros((1, 18), dtype=bool),
        incoming_masks=incoming,
        raw_normals=np.zeros((1, 3), dtype=np.float64),
        normal_indices=np.asarray((-1,), dtype=np.int64),
    )
    qvalues = np.full((1, 18), np.nan)
    qvalues[0, inverse] = qvalue
    mesh = MeshContract(
        mesh=Path("."),
        tree_ids=np.asarray((1,), dtype=np.int64),
        cell_ijk=np.asarray(((0, 0, 0),), dtype=np.int64),
        lookup={(0, 0, 0): 0},
        boundaries={"wall": boundary},
        boundary_labels=("wall",),
        qvalues_by_cell=qvalues,
    )
    return np.arange(1.0, 20.0).reshape(1, 19), mesh, inverse


@pytest.mark.parametrize("qvalue", (0.25, 0.75))
def test_referee_wall_libb_matches_pinned_equation_for_arbitrary_q(
    qvalue: float,
) -> None:
    pdf, mesh, inverse = _single_wall_mesh(qvalue)

    operation = _wall_operations(pdf, mesh)[0]
    expected = libb_reference_replacement(
        qvalue,
        f_in=float(pdf[0, 0]),
        f_out=float(pdf[0, inverse]),
        f_neigh=float(pdf[0, 0]),
    )

    assert operation[:3] == ("wall", 0, inverse)
    assert operation[3] == pytest.approx(expected)
