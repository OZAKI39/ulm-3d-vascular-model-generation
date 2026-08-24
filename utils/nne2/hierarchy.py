"""Infer an anatomical parent-child hierarchy from NNE2 diving-trunk anchors."""

from __future__ import annotations

import logging
from collections import defaultdict

import networkx as nx
import numpy as np

from ..graph.model import BranchRecord, HierarchicalGraphResult
from .model import (
    AnchorMatch,
    DirectedBranch,
    NNE2HierarchyResult,
    ParentChildRelation,
)


def _length(branch: BranchRecord) -> float:
    return float(branch.arc_length_raw_um[-1]) if len(branch.arc_length_raw_um) else 0.0


def _mean_radius(branch: BranchRecord) -> float:
    return float(np.mean(branch.coarse_radius_raw_um)) if len(branch.coarse_radius_raw_um) else 0.0


def _component_branch_ids(
    graph: nx.MultiGraph, component_nodes: set[int]
) -> set[int]:
    output: set[int] = set()
    for left, right, key, data in graph.edges(keys=True, data=True):
        if left in component_nodes and right in component_nodes:
            output.add(int(data.get("branch_id", key)))
    return output


def build_directed_hierarchy(
    tree_key: str,
    stack_name: str,
    graph_result: HierarchicalGraphResult,
    anchors: list[AnchorMatch],
    logger: logging.Logger | None = None,
) -> NNE2HierarchyResult:
    logger = logger or logging.getLogger("ulm_3d_vascular")
    branch_by_id = {item.branch_id: item for item in graph_result.branches}
    trunk_anchors = [
        item
        for item in anchors
        if item.branching_order == 0 and item.matched_branch_id in branch_by_id
    ]
    if not trunk_anchors:
        raise ValueError(f"{tree_key}/{stack_name} has no matched Branching Order 0 anchor")
    root_anchor = min(
        trunk_anchors,
        key=lambda item: (item.depth_um, -item.registration_score, item.record_id),
    )
    root_branch = branch_by_id[int(root_anchor.matched_branch_id)]
    node_u = graph_result.nodes[root_branch.node_u]
    node_v = graph_result.nodes[root_branch.node_v]
    root_node = min(
        (node_u, node_v),
        key=lambda item: (
            item.representative_lps_um[2],
            item.graph_degree,
            item.node_id,
        ),
    ).node_id

    all_graph = graph_result.junction_graph
    component_nodes = set(nx.node_connected_component(all_graph, root_node))
    retained_ids = _component_branch_ids(all_graph, component_nodes)
    excluded_ids = sorted(set(branch_by_id) - retained_ids)
    distances = nx.single_source_dijkstra_path_length(all_graph, root_node, weight="length_raw_um")

    directed: list[DirectedBranch] = []
    for branch_id in sorted(retained_ids):
        branch = branch_by_id[branch_id]
        left_distance = float(distances.get(branch.node_u, np.inf))
        right_distance = float(distances.get(branch.node_v, np.inf))
        tolerance = max(1.0e-6, 1.0e-6 * max(left_distance, right_distance, 1.0))
        if branch.node_u == branch.node_v or abs(left_distance - right_distance) <= tolerance:
            upstream = downstream = None
            status = "unresolved_equal_root_distance"
            confidence = "unresolved"
        elif left_distance < right_distance:
            upstream, downstream = branch.node_u, branch.node_v
            status = "inferred_away_from_diving_trunk_root"
            confidence = "medium"
        else:
            upstream, downstream = branch.node_v, branch.node_u
            status = "inferred_away_from_diving_trunk_root"
            confidence = "medium"
        directed.append(
            DirectedBranch(
                branch_id=branch_id,
                node_u=branch.node_u,
                node_v=branch.node_v,
                upstream_node=upstream,
                downstream_node=downstream,
                direction_status=status,
                confidence=confidence,
                branching_order=None,
                order_source="unassigned",
            )
        )
    directed_by_id = {item.branch_id: item for item in directed}

    anchor_orders: dict[int, list[int]] = defaultdict(list)
    anchor_records: dict[int, list[int]] = defaultdict(list)
    anchor_registered: dict[int, bool] = defaultdict(lambda: True)
    for anchor in anchors:
        if anchor.matched_branch_id in retained_ids:
            branch_id = int(anchor.matched_branch_id)
            anchor_orders[branch_id].append(anchor.branching_order)
            anchor_records[branch_id].append(anchor.record_id)
            anchor_registered[branch_id] &= anchor.registration_status == "registered"
    order_conflicts: set[int] = set()
    for branch_id, values in anchor_orders.items():
        item = directed_by_id[branch_id]
        item.anchor_record_ids = sorted(anchor_records[branch_id])
        unique = sorted(set(values))
        if len(unique) == 1:
            item.branching_order = unique[0]
            item.order_source = "NNE2_measurement_anchor"
            item.confidence = "high" if anchor_registered[branch_id] else "medium"
        else:
            item.order_source = "conflicting_NNE2_measurement_anchors"
            item.confidence = "low"
            order_conflicts.add(branch_id)
    directed_by_id[root_branch.branch_id].branching_order = 0
    if root_branch.branch_id not in order_conflicts:
        directed_by_id[root_branch.branch_id].order_source = "diving_trunk_root_anchor"

    incoming: dict[int, list[int]] = defaultdict(list)
    outgoing: dict[int, list[int]] = defaultdict(list)
    for item in directed:
        if item.upstream_node is not None and item.downstream_node is not None:
            outgoing[item.upstream_node].append(item.branch_id)
            incoming[item.downstream_node].append(item.branch_id)

    primary_branch_for_node: dict[int, int] = {}
    primary_ids: set[int] = set()
    for node_id in sorted(component_nodes, key=lambda value: distances.get(value, np.inf)):
        if node_id == root_node:
            continue
        candidates = incoming.get(node_id, [])
        if not candidates:
            continue
        best = min(
            candidates,
            key=lambda branch_id: (
                abs(
                    distances[node_id]
                    - distances[directed_by_id[branch_id].upstream_node]  # type: ignore[index]
                    - _length(branch_by_id[branch_id])
                ),
                -_mean_radius(branch_by_id[branch_id]),
                branch_id,
            ),
        )
        primary_branch_for_node[node_id] = best
        primary_ids.add(best)
    primary_ids.add(root_branch.branch_id)
    for item in directed:
        item.is_primary_tree_edge = item.branch_id in primary_ids
        item.is_cross_link = (
            item.upstream_node is not None
            and item.downstream_node is not None
            and item.branch_id not in primary_ids
        )
        if item.is_cross_link:
            item.confidence = "low"

    # Propagate generation order along the primary tree. The largest-radius unlabeled
    # continuation keeps the parent's order; side branches advance one order.
    visited: set[int] = set()
    root_outgoing = [
        branch_id for branch_id in outgoing.get(root_node, []) if branch_id in primary_ids
    ]
    if root_branch.branch_id not in root_outgoing:
        root_outgoing.insert(0, root_branch.branch_id)
    for branch_id in root_outgoing:
        branch = directed_by_id[branch_id]
        if branch.branching_order is None:
            branch.branching_order = 0 if branch_id == root_branch.branch_id else 1
            branch.order_source = (
                "diving_trunk_root_anchor"
                if branch_id == root_branch.branch_id
                else "inferred_side_branch_at_root"
            )
    queue = list(root_outgoing)
    directed_by_id[root_branch.branch_id].branching_order = 0
    while queue:
        parent_id = queue.pop(0)
        if parent_id in visited:
            continue
        visited.add(parent_id)
        parent = directed_by_id[parent_id]
        if parent.downstream_node is None:
            continue
        child_ids = [
            child_id
            for child_id in outgoing.get(parent.downstream_node, [])
            if child_id in primary_ids and child_id != parent_id
        ]
        unlabeled = [
            child_id for child_id in child_ids if directed_by_id[child_id].branching_order is None
        ]
        continuation = (
            max(unlabeled, key=lambda value: (_mean_radius(branch_by_id[value]), -value))
            if unlabeled
            else None
        )
        for child_id in child_ids:
            child = directed_by_id[child_id]
            proposed = min(
                4,
                (parent.branching_order or 0) + (0 if child_id == continuation else 1),
            )
            if child.branching_order is None:
                child.branching_order = proposed
                child.order_source = (
                    "inferred_largest_radius_continuation"
                    if child_id == continuation
                    else "inferred_side_branch_generation"
                )
            elif parent.branching_order is not None and child.branching_order < parent.branching_order:
                order_conflicts.add(child_id)
                child.confidence = "low"
                child.order_source += ";conflicts_with_parent_order"
            queue.append(child_id)

    # A direction-resolved cross-link is not part of the one-parent tree, but it still
    # receives a coarse order from its upstream neighborhood for visualization/export.
    for item in sorted(
        directed,
        key=lambda value: (
            distances.get(value.upstream_node, np.inf)
            if value.upstream_node is not None
            else np.inf,
            value.branch_id,
        ),
    ):
        if item.branching_order is not None or item.upstream_node is None:
            continue
        upstream_orders = [
            directed_by_id[parent_id].branching_order
            for parent_id in incoming.get(item.upstream_node, [])
            if directed_by_id[parent_id].branching_order is not None
        ]
        if upstream_orders:
            item.branching_order = min(4, min(upstream_orders) + 1)
            item.order_source = "inferred_cross_link_from_upstream_neighborhood"

    relations: list[ParentChildRelation] = []
    for child in directed:
        if child.upstream_node is None or child.branch_id == root_branch.branch_id:
            continue
        parents = [
            branch_id
            for branch_id in incoming.get(child.upstream_node, [])
            if branch_id != child.branch_id
        ]
        if not parents:
            continue
        primary_parent = primary_branch_for_node.get(child.upstream_node)
        if primary_parent not in parents:
            primary_parent = min(
                parents,
                key=lambda value: (
                    directed_by_id[value].branching_order
                    if directed_by_id[value].branching_order is not None
                    else 10_000,
                    value,
                ),
            )
        for parent_id in parents:
            relation_type = (
                "primary_parent" if parent_id == primary_parent and child.is_primary_tree_edge
                else "cross_link"
            )
            confidence = (
                "medium"
                if relation_type == "primary_parent"
                and parent_id not in order_conflicts
                and child.branch_id not in order_conflicts
                else "low"
            )
            relations.append(
                ParentChildRelation(
                    parent_branch_id=parent_id,
                    child_branch_id=child.branch_id,
                    shared_node_id=child.upstream_node,
                    relation_type=relation_type,
                    confidence=confidence,
                )
            )
            directed_by_id[parent_id].child_branch_ids.append(child.branch_id)
            child.parent_branch_ids.append(parent_id)

    directed_graph = nx.MultiDiGraph(
        representation="NNE2_directed_junction_branch_graph",
        direction_source="inferred_from_diving_trunk_and_branching_order",
        measured_flow_direction=False,
    )
    for node_id in component_nodes:
        source = graph_result.nodes[node_id]
        x, y, z = source.representative_lps_um
        directed_graph.add_node(node_id, x_um=x, y_um=y, z_um=z, node_type=source.node_type)
    for item in directed:
        attributes = item.to_dict()
        if item.upstream_node is None or item.downstream_node is None:
            directed_graph.add_edge(item.node_u, item.node_v, key=item.branch_id, **attributes)
        else:
            directed_graph.add_edge(
                item.upstream_node, item.downstream_node, key=item.branch_id, **attributes
            )

    primary_tree = nx.DiGraph(
        representation="NNE2_primary_parent_tree",
        measured_flow_direction=False,
    )
    primary_tree.add_nodes_from(directed_graph.nodes(data=True))
    for item in directed:
        if item.is_primary_tree_edge and item.upstream_node is not None and item.downstream_node is not None:
            primary_tree.add_edge(item.upstream_node, item.downstream_node, branch_id=item.branch_id)

    unresolved = sorted(
        item.branch_id for item in directed if item.upstream_node is None or item.downstream_node is None
    )
    cross_links = sorted(item.branch_id for item in directed if item.is_cross_link)
    logger.info(
        "Directed hierarchy %s/%s: root branch=%d, retained=%d, excluded islands=%d, "
        "cross-links=%d, unresolved=%d",
        tree_key,
        stack_name,
        root_branch.branch_id,
        len(directed),
        len(excluded_ids),
        len(cross_links),
        len(unresolved),
    )
    return NNE2HierarchyResult(
        tree_key=tree_key,
        stack_name=stack_name,
        root_node_id=root_node,
        root_branch_id=root_branch.branch_id,
        branches=directed,
        relations=relations,
        anchors=anchors,
        directed_graph=directed_graph,
        primary_tree=primary_tree,
        unresolved_branch_ids=unresolved,
        cross_link_branch_ids=cross_links,
        order_conflict_branch_ids=sorted(order_conflicts),
        excluded_branch_ids=excluded_ids,
    )
