from __future__ import annotations

from schmid_test_data import write_synthetic_schmid
from utils.schmid_pkl.cleanup import clean_schmid_input
from utils.schmid_pkl.config import SchmidPKLConfig
from utils.schmid_pkl.graph_builder import build_directed_hierarchical_graph
from utils.schmid_pkl.loader import load_schmid_input


def _build(tmp_path, *, equal_pressure_edge: bool = False):
    input_dir = write_synthetic_schmid(
        tmp_path / "NW1_results", equal_pressure_edge=equal_pressure_edge
    )
    config = SchmidPKLConfig(input_dir=input_dir, output_root=tmp_path / "outputs")
    cleanup = clean_schmid_input(load_schmid_input(input_dir), config)
    return build_directed_hierarchical_graph(cleanup, config)


def test_pressure_direction_overrides_tuple_storage_order(tmp_path) -> None:
    result = _build(tmp_path)
    first = next(edge for edge in result.cleanup.edges if edge.edge_id == 0)

    assert (first.node_u, first.node_v) == (1, 0)
    assert (first.upstream_node, first.downstream_node) == (0, 1)
    assert first.direction_status == "known"


def test_degree_two_chains_merge_without_losing_split_and_merge(tmp_path) -> None:
    result = _build(tmp_path)

    assert len(result.branches) == 3
    assert set(result.raw_edge_to_branch) == set(range(5))
    roles = {node.node_id: node.node_role for node in result.nodes}
    assert roles[0] == "source"
    assert roles[1] == "split"
    assert roles[4] == "merge_sink"
    trunk = next(branch for branch in result.branches if branch.upstream_node == 0)
    assert len(trunk.child_branch_ids) == 2
    assert all(trunk.branch_id in result.branches[item].parent_branch_ids for item in trunk.child_branch_ids)


def test_equal_pressure_edge_remains_explicitly_unresolved(tmp_path) -> None:
    result = _build(tmp_path, equal_pressure_edge=True)
    first = next(edge for edge in result.cleanup.edges if edge.edge_id == 0)

    assert first.upstream_node is None
    assert first.direction_status.startswith("unresolved_equal_pressure")
    assert result.raw_unresolved_edge_ids == [0]
