from __future__ import annotations

from utils.cfd_flow.repaired_topology_forensics import (
    classify_component_topology,
    unit_scaling_oracle_decision,
)


def _ports(component: int) -> dict:
    result = {
        label: {"component_ids": [component]}
        for label in ("inlet", "outlet_01", "outlet_02", "outlet_03")
    }
    return result


def test_decision_tree_identifies_d3q19_diagonal_connection() -> None:
    assert (
        classify_component_topology((90, 10), (100,), _ports(0))
        == "D3Q19_DIAGONAL_ONLY_CONNECTION"
    )


def test_decision_tree_does_not_delete_portless_pocket() -> None:
    assert (
        classify_component_topology((9999, 1), (9999, 1), _ports(0))
        == "PORTLESS_ISOLATED_POCKET"
    )


def test_decision_tree_marks_large_secondary_as_major_split() -> None:
    assert (
        classify_component_topology((900, 100), (900, 100), _ports(0))
        == "MAJOR_NETWORK_SPLIT"
    )


def test_scale_oracle_requires_exact_topology() -> None:
    comparison = {
        "tree_id_set_exact_match": False,
        "boundary_id_exact_match": True,
        "q_patched_minus_scaled": {"max": 0.0},
    }

    assert (
        unit_scaling_oracle_decision(
            comparison,
            {"6": [10], "18": [10], "26": [10]},
            {"6": [9, 1], "18": [10], "26": [10]},
            derived_q_tolerance=1.0e-5,
        )
        == "FAIL"
    )


def test_scale_oracle_accepts_geometry_bounded_float32_roundoff() -> None:
    comparison = {
        "tree_id_set_exact_match": True,
        "boundary_id_exact_match": True,
        "q_patched_minus_scaled": {"max": 5.0e-6},
    }
    structure = {"6": [10], "18": [10], "26": [10]}

    assert (
        unit_scaling_oracle_decision(
            comparison,
            structure,
            structure,
            derived_q_tolerance=1.0e-5,
        )
        == "PASS_WITH_FLOAT32_STL_ROUNDOFF"
    )
