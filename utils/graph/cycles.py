"""Cycle-preserving analysis for a vascular multigraph."""

from __future__ import annotations

from collections import defaultdict

import networkx as nx

from .model import BranchRecord, CycleRecord, NodeRecord


def identify_cycles(
    graph: nx.MultiGraph,
    nodes: list[NodeRecord],
    branches: list[BranchRecord],
) -> tuple[list[CycleRecord], int]:
    components = nx.number_connected_components(graph) if graph.number_of_nodes() else 0
    cycle_rank = graph.number_of_edges() - graph.number_of_nodes() + components
    records: list[CycleRecord] = []
    by_pair: dict[tuple[int, int], list[int]] = defaultdict(list)
    for branch in branches:
        key = tuple(sorted((branch.node_u, branch.node_v)))
        by_pair[key].append(branch.branch_id)

    for (node_u, node_v), branch_ids in sorted(by_pair.items()):
        if node_u == node_v:
            for branch_id in branch_ids:
                records.append(
                    CycleRecord(len(records), "self_loop", [node_u], [branch_id])
                )
        elif len(branch_ids) > 1:
            base = branch_ids[0]
            for branch_id in branch_ids[1:]:
                records.append(
                    CycleRecord(
                        len(records),
                        "parallel_path",
                        [node_u, node_v],
                        [base, branch_id],
                    )
                )

    simple_graph = nx.Graph()
    simple_graph.add_nodes_from(graph.nodes)
    simple_graph.add_edges_from(pair for pair in by_pair if pair[0] != pair[1])
    for node_cycle in nx.cycle_basis(simple_graph):
        branch_ids: list[int] = []
        for index, node_u in enumerate(node_cycle):
            node_v = node_cycle[(index + 1) % len(node_cycle)]
            branch_ids.append(by_pair[tuple(sorted((node_u, node_v)))][0])
        records.append(
            CycleRecord(len(records), "simple_cycle", list(node_cycle), branch_ids)
        )

    branch_by_id = {item.branch_id: item for item in branches}
    node_by_id = {item.node_id: item for item in nodes}
    for record in records:
        for branch_id in record.branch_ids:
            branch_by_id[branch_id].cycle_ids.append(record.cycle_id)
        for node_id in record.node_ids:
            node_by_id[node_id].cycle_ids.append(record.cycle_id)
    return records, int(cycle_rank)
