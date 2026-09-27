from __future__ import annotations

import networkx as nx
import numpy as np
from scipy import ndimage

from utils.nne2.centerline import build_centerline_graph
from utils.nne2.hierarchy import build_directed_hierarchy
from utils.nne2.model import AnchorMatch


def _synthetic_graph():
    mask = np.zeros((30, 32, 32), dtype=bool)  # ZYX
    mask[2:25, 16, 16] = True
    for offset in range(11):
        mask[14 + offset, 16, 16 + offset] = True
        mask[14 + offset, 16, 16 - offset] = True
    mask[3:10, 3, 3] = True  # disconnected island
    mask = ndimage.binary_dilation(mask, iterations=1)
    return build_centerline_graph(mask, (1.0, 1.0, 1.0), graph_connectivity=26)


def test_diving_trunk_anchor_builds_acyclic_primary_parent_tree() -> None:
    centerline = _synthetic_graph()
    root_branch = min(
        centerline.graph.branches,
        key=lambda item: float(np.min(item.points_raw_lps_um[:, 2])),
    )
    anchor = AnchorMatch(
        record_id=1,
        subject_id="010101",
        tree_id=1,
        branching_order=0,
        depth_um=20.0,
        expected_stack_index=3,
        matched_stack_index=3,
        seed_xyz_um=tuple(float(v) for v in root_branch.points_raw_lps_um[0]),
        registration_score=0.9,
        registration_status="registered",
        anchor_pixel_method="test",
        matched_branch_id=root_branch.branch_id,
        branch_distance_um=0.0,
        match_status="matched",
    )
    result = build_directed_hierarchy(
        "010101_tree_1", "Stack-Test", centerline.graph, [anchor]
    )

    assert result.root_branch_id == root_branch.branch_id
    assert nx.is_directed_acyclic_graph(result.primary_tree)
    assert result.excluded_branch_ids
    assert next(
        item for item in result.branches if item.branch_id == result.root_branch_id
    ).branching_order == 0
    assert all(
        item.branching_order is not None
        for item in result.branches
        if item.is_primary_tree_edge and item.upstream_node is not None
    )


def test_hierarchy_requires_a_matched_diving_trunk_anchor() -> None:
    centerline = _synthetic_graph()
    anchor = AnchorMatch(
        record_id=1,
        subject_id="010101",
        tree_id=1,
        branching_order=1,
        depth_um=20.0,
        expected_stack_index=3,
        matched_stack_index=3,
        seed_xyz_um=(0.0, 0.0, 0.0),
        registration_score=0.9,
        registration_status="registered",
        anchor_pixel_method="test",
        matched_branch_id=centerline.graph.branches[0].branch_id,
        branch_distance_um=0.0,
        match_status="matched",
    )

    try:
        build_directed_hierarchy("010101_tree_1", "Stack-Test", centerline.graph, [anchor])
    except ValueError as exc:
        assert "Branching Order 0" in str(exc)
    else:
        raise AssertionError("Expected a missing diving-trunk error")
