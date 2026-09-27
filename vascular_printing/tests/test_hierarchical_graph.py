from __future__ import annotations

import json
from pathlib import Path

import networkx as nx
import numpy as np

from utils.config import HierarchicalGraphConfig
from utils.graph.export import export_hierarchical_graph, verify_graph_exports
from utils.graph.extraction import build_hierarchical_graph
from utils.graph.validation import evaluate_hierarchical_graph_acceptance


def _build(
    skeleton: np.ndarray,
    mask: np.ndarray | None = None,
    *,
    connectivity: int = 6,
    spacing: tuple[float, float, float] = (1.0, 1.0, 1.0),
):
    config = HierarchicalGraphConfig(
        neighbor_connectivity=connectivity,
        smoothing_enabled=False,
        short_branch_warning_um=0.0,
    )
    result = build_hierarchical_graph(
        skeleton,
        skeleton if mask is None else mask,
        (10.0, 20.0, 30.0),
        spacing,
        config,
    )
    return result, config


def test_straight_chain_becomes_one_geometry_aware_branch() -> None:
    skeleton = np.zeros((7, 5, 5), dtype=bool)
    skeleton[1:6, 2, 2] = True
    result, _ = _build(skeleton)

    assert len(result.nodes) == 2
    assert len(result.branches) == 1
    assert all(node.node_type == "terminal" for node in result.nodes)
    assert len(result.branches[0].points_raw_lps_um) == 5
    assert result.branches[0].arc_length_raw_um[-1] == 4.0
    assert np.array_equal(result.reconstructed_skeleton, skeleton)
    assert result.missing_voxel_count == 0
    assert result.extra_voxel_count == 0


def test_y_junction_preserves_three_branches_and_branch_as_node_relations() -> None:
    skeleton = np.zeros((7, 7, 7), dtype=bool)
    skeleton[1:6, 3, 3] = True
    skeleton[3, 4:6, 3] = True
    result, _ = _build(skeleton)

    assert len(result.nodes) == 4
    assert len(result.branches) == 3
    assert sum(node.node_type == "junction" for node in result.nodes) == 1
    assert sum(node.node_type == "terminal" for node in result.nodes) == 3
    assert result.branch_as_node_graph.number_of_nodes() == 3
    assert result.branch_as_node_graph.number_of_edges() == 3
    assert len(result.junctions) == 1
    assert len(result.junctions[0].pairwise_angles_deg) == 3
    assert np.array_equal(result.reconstructed_skeleton, skeleton)


def test_diagonal_chain_uses_physical_euclidean_step_length() -> None:
    skeleton = np.zeros((5, 5, 5), dtype=bool)
    skeleton[1, 1, 1] = True
    skeleton[2, 2, 2] = True
    skeleton[3, 3, 3] = True
    result, _ = _build(skeleton, connectivity=26, spacing=(2.0, 2.0, 2.0))

    assert len(result.branches) == 1
    assert np.isclose(result.branches[0].arc_length_raw_um[-1], 4.0 * np.sqrt(3.0))


def test_pure_degree_two_cycle_keeps_cycle_identity() -> None:
    skeleton = np.zeros((7, 7, 7), dtype=bool)
    coordinates = [
        (2, 2, 3),
        (3, 2, 3),
        (4, 2, 3),
        (4, 3, 3),
        (4, 4, 3),
        (3, 4, 3),
        (2, 4, 3),
        (2, 3, 3),
    ]
    for coordinate in coordinates:
        skeleton[coordinate] = True
    result, _ = _build(skeleton)

    assert len(result.nodes) == 1
    assert result.nodes[0].node_type == "cycle_anchor"
    assert len(result.branches) == 1
    assert result.branches[0].node_u == result.branches[0].node_v
    assert result.cycle_rank == 1
    assert len(result.cycles) == 1
    assert np.array_equal(result.reconstructed_skeleton, skeleton)


def test_coarse_radius_sequence_is_retained_and_labeled() -> None:
    skeleton = np.zeros((9, 9, 9), dtype=bool)
    skeleton[2:7, 4, 4] = True
    mask = np.zeros_like(skeleton)
    for x in range(1, 8):
        for y in range(2, 7):
            for z in range(2, 7):
                if (y - 4) ** 2 + (z - 4) ** 2 <= 4:
                    mask[x, y, z] = True
    result, _ = _build(skeleton, mask)

    branch = result.branches[0]
    assert len(branch.coarse_radius_raw_um) == len(branch.points_raw_lps_um)
    assert np.all(branch.coarse_radius_raw_um > 0)
    assert "coarse navigation" in result.radius_source
    assert result.report()["approved_for_cfd"] is False


def test_portable_exports_reopen_with_matching_counts(tmp_path: Path) -> None:
    skeleton = np.zeros((7, 7, 7), dtype=bool)
    skeleton[1:6, 3, 3] = True
    skeleton[3, 4:6, 3] = True
    result, config = _build(skeleton)
    paths = export_hierarchical_graph(
        result,
        tmp_path / "graphs",
        tmp_path / "tables",
        save_graphml=True,
        save_vtp=True,
        save_npz=True,
    )
    errors = verify_graph_exports(paths, len(result.nodes), len(result.branches))

    assert errors == []
    payload = json.loads(
        (tmp_path / "graphs" / "hierarchical_vascular_graph.json").read_text(
            encoding="utf-8"
        )
    )
    assert payload["scale_3_dense_geometry"]["branches"][0]["points_raw_lps_um"]
    graph = nx.read_graphml(tmp_path / "graphs" / "branch_as_node_graph.graphml")
    assert graph.number_of_nodes() == len(result.branches)
    acceptance = evaluate_hierarchical_graph_acceptance(result, config, paths, errors)
    assert acceptance.overall_status == "PASS"
