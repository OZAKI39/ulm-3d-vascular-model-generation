from __future__ import annotations

from utils.cfd_flow.repaired_topology_forensics import classify_component_topology


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
