"""Build graph views used for topology review and later GNN work."""

from __future__ import annotations

from itertools import combinations

import networkx as nx

from .model import BranchRecord, NodeRecord


def build_junction_graph(
    nodes: list[NodeRecord], branches: list[BranchRecord]
) -> nx.MultiGraph:
    graph = nx.MultiGraph(
        representation="junction_as_node",
        directed_for_flow=False,
        relationship="undirected vascular connectivity",
    )
    for node in nodes:
        x, y, z = node.representative_lps_um
        graph.add_node(
            node.node_id,
            node_type=node.node_type,
            graph_degree=node.graph_degree,
            x_lps_um=x,
            y_lps_um=y,
            z_lps_um=z,
            voxel_count=len(node.voxel_indices_xyz),
        )
    for branch in branches:
        summary = branch.summary()
        scalar_summary = {
            key: value
            for key, value in summary.items()
            if key not in {"branch_id", "node_u", "node_v"}
            and isinstance(value, (str, int, float, bool))
            and value is not None
        }
        graph.add_edge(
            branch.node_u,
            branch.node_v,
            key=branch.branch_id,
            branch_id=branch.branch_id,
            **scalar_summary,
        )
    return graph


def build_branch_as_node_graph(
    nodes: list[NodeRecord], branches: list[BranchRecord]
) -> nx.MultiGraph:
    graph = nx.MultiGraph(
        representation="branch_as_node",
        directed_for_flow=False,
        relationship="shares_junction",
    )
    for branch in branches:
        summary = branch.summary()
        attributes = {
            key: value
            for key, value in summary.items()
            if isinstance(value, (str, int, float, bool)) and value is not None
        }
        graph.add_node(branch.branch_id, **attributes)
    for node in nodes:
        endpoint_records: list[tuple[int, str]] = []
        for branch in branches:
            if branch.node_u == node.node_id:
                endpoint_records.append((branch.branch_id, "start"))
            if branch.node_v == node.node_id:
                endpoint_records.append((branch.branch_id, "end"))
        for relation_index, (left, right) in enumerate(combinations(endpoint_records, 2)):
            branch_a, endpoint_a = left
            branch_b, endpoint_b = right
            if branch_a == branch_b:
                continue
            graph.add_edge(
                branch_a,
                branch_b,
                key=f"{node.node_id}_{relation_index}",
                shared_node_id=node.node_id,
                shared_node_type=node.node_type,
                branch_a_endpoint=endpoint_a,
                branch_b_endpoint=endpoint_b,
                relationship="shares_junction",
                flow_direction_known=False,
            )
    return graph
