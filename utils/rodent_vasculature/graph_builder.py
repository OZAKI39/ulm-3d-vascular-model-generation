"""Aggregate a parent-to-current SWC forest into directed vascular branches."""

from __future__ import annotations

from collections import defaultdict

import networkx as nx
import numpy as np

from .config import RodentVasculatureConfig
from .geometry import derive_branch_geometry
from .model import DirectedBranch, DirectedVascularGraph
from .swc_io import SWCData


def _node_role(graph: nx.DiGraph, node_id: int) -> str:
    indegree, outdegree = graph.in_degree(node_id), graph.out_degree(node_id)
    if indegree == 0:
        return "inferred_inlet"
    if outdegree == 0:
        return "inferred_outlet"
    if outdegree >= 2:
        return "divergence_junction"
    if indegree >= 2:
        return "convergence_junction"
    return "connector"


def _component_metadata(graph: nx.DiGraph) -> tuple[dict[int, int], dict[int, int]]:
    component_by_node: dict[int, int] = {}
    root_by_node: dict[int, int] = {}
    for component_id, nodes in enumerate(nx.weakly_connected_components(graph)):
        ordered = sorted(int(node) for node in nodes)
        roots = [node for node in ordered if graph.in_degree(node) == 0]
        root = roots[0] if roots else ordered[0]
        for node in ordered:
            component_by_node[node] = component_id
            root_by_node[node] = root
    return component_by_node, root_by_node


def _branch_hierarchy(branches: list[DirectedBranch], graph: nx.DiGraph) -> None:
    by_id = {branch.branch_id: branch for branch in branches}
    for branch in branches:
        branch.parent_branch_ids = sorted(graph.predecessors(branch.branch_id))
        branch.daughter_branch_ids = sorted(graph.successors(branch.branch_id))

    if not nx.is_directed_acyclic_graph(graph):
        raise ValueError("Aggregated branch graph contains a directed cycle")
    for branch_id in nx.topological_sort(graph):
        branch = by_id[branch_id]
        branch.depth = (
            max(by_id[parent].depth for parent in branch.parent_branch_ids) + 1
            if branch.parent_branch_ids
            else 0
        )
    for branch_id in reversed(list(nx.topological_sort(graph))):
        branch = by_id[branch_id]
        daughters = [by_id[value] for value in branch.daughter_branch_ids]
        if not daughters:
            branch.strahler_order = 1
            branch.horsfield_order = 1
            branch.downstream_terminal_count = 1
            branch.downstream_branch_count = 1
            continue
        orders = [daughter.strahler_order for daughter in daughters]
        maximum = max(orders)
        branch.strahler_order = maximum + 1 if orders.count(maximum) >= 2 else maximum
        branch.horsfield_order = 1 + max(daughter.horsfield_order for daughter in daughters)
        branch.downstream_terminal_count = sum(
            daughter.downstream_terminal_count for daughter in daughters
        )
        branch.downstream_branch_count = 1 + sum(
            daughter.downstream_branch_count for daughter in daughters
        )


def build_directed_vascular_graph(
    sample_id: str, swc: SWCData, config: RodentVasculatureConfig
) -> DirectedVascularGraph:
    source = swc.directed_graph()
    if not swc.structurally_valid:
        raise ValueError(f"SWC structural validation failed: {swc.validation}")
    if not nx.is_directed_acyclic_graph(source):
        raise ValueError("SWC parent graph is not acyclic")

    index_by_id = {int(node): index for index, node in enumerate(swc.node_ids.tolist())}
    component_by_node, root_by_node = _component_metadata(source)
    critical = {
        int(node)
        for node in source.nodes
        if source.in_degree(node) != 1 or source.out_degree(node) != 1
    }
    junction_graph = nx.MultiDiGraph()
    for node_id in sorted(critical):
        index = index_by_id[node_id]
        junction_graph.add_node(
            node_id,
            role=_node_role(source, node_id),
            component_id=component_by_node[node_id],
            root_node_id=root_by_node[node_id],
            source_index=index,
            x_um=float(swc.points_um[index, 0]),
            y_um=float(swc.points_um[index, 1]),
            z_um=float(swc.points_um[index, 2]),
            radius_raw_um=float(swc.radius_raw_um[index]),
            indegree=int(source.in_degree(node_id)),
            outdegree=int(source.out_degree(node_id)),
        )

    branches: list[DirectedBranch] = []
    represented_edges: set[tuple[int, int]] = set()
    source_to_branches: dict[int, list[int]] = defaultdict(list)
    shape_xyz = np.asarray(config.expected_shape_zyx[::-1], dtype=float) if config.expected_shape_zyx else None
    for upstream in sorted(critical):
        for successor in sorted(source.successors(upstream)):
            sequence = [upstream, int(successor)]
            while sequence[-1] not in critical:
                next_nodes = list(source.successors(sequence[-1]))
                if len(next_nodes) != 1:
                    raise ValueError(f"Non-critical SWC node has {len(next_nodes)} children")
                sequence.append(int(next_nodes[0]))
            for edge in zip(sequence[:-1], sequence[1:]):
                if edge in represented_edges:
                    raise ValueError(f"SWC edge represented more than once: {edge}")
                represented_edges.add(edge)
            indices = [index_by_id[node] for node in sequence]
            raw_voxel = swc.points_voxel_xyz[indices].copy()
            raw_um = swc.points_um[indices].copy()
            raw_radius = swc.radius_raw_um[indices].copy()
            geometry = derive_branch_geometry(
                raw_um,
                raw_radius,
                smoothing_enabled=config.smoothing_enabled,
                smoothing_window_points=config.smoothing_window_points,
                resample_step_um=config.resample_step_um,
            )
            branch_id = len(branches)
            proximal_boundary = bool(
                shape_xyz is not None
                and np.any((raw_voxel[0] <= 1.0) | (raw_voxel[0] >= shape_xyz - 1.0))
            )
            distal_boundary = bool(
                shape_xyz is not None
                and np.any((raw_voxel[-1] <= 1.0) | (raw_voxel[-1] >= shape_xyz - 1.0))
            )
            branch = DirectedBranch(
                branch_id=branch_id,
                upstream_node_id=sequence[0],
                downstream_node_id=sequence[-1],
                source_node_ids=sequence,
                raw_points_voxel_xyz=raw_voxel,
                raw_points_um=raw_um,
                raw_radius_um=raw_radius,
                derived_points_um=geometry[0],
                derived_radius_um=geometry[1],
                direction_vectors_xyz=geometry[2],
                raw_length_um=geometry[3],
                derived_length_um=geometry[4],
                tortuosity=geometry[5],
                mean_radius_um=geometry[6],
                component_id=component_by_node[upstream],
                root_node_id=root_by_node[upstream],
                proximal_boundary_contact=proximal_boundary,
                distal_boundary_contact=distal_boundary,
            )
            branches.append(branch)
            junction_graph.add_edge(
                branch.upstream_node_id,
                branch.downstream_node_id,
                key=branch_id,
                branch_id=branch_id,
                direction="parent_to_current",
            )
            for node_id in sequence:
                source_to_branches[node_id].append(branch_id)

    missing_edges = set(source.edges) - represented_edges
    if missing_edges:
        raise ValueError(f"SWC edges omitted during aggregation: {sorted(missing_edges)[:10]}")
    branch_graph = nx.DiGraph()
    branch_graph.add_nodes_from(branch.branch_id for branch in branches)
    incoming_by_node: dict[int, list[int]] = defaultdict(list)
    outgoing_by_node: dict[int, list[int]] = defaultdict(list)
    for branch in branches:
        incoming_by_node[branch.downstream_node_id].append(branch.branch_id)
        outgoing_by_node[branch.upstream_node_id].append(branch.branch_id)
    for node_id in critical:
        for parent_branch in incoming_by_node[node_id]:
            for daughter_branch in outgoing_by_node[node_id]:
                branch_graph.add_edge(parent_branch, daughter_branch, junction_node_id=node_id)
    _branch_hierarchy(branches, branch_graph)

    warnings: list[str] = []
    invalid_radius_count = len(swc.validation["nonpositive_radius_node_ids"])
    if invalid_radius_count:
        warnings.append(
            f"{invalid_radius_count} SWC radii are non-positive/non-finite; raw values were preserved "
            "and only derived geometry used interpolation where possible."
        )
    zero_length = [branch.branch_id for branch in branches if branch.raw_length_um <= 1.0e-12]
    if zero_length:
        warnings.append(f"Zero-length branches: {zero_length[:20]}")
    return DirectedVascularGraph(
        sample_id,
        swc,
        source,
        junction_graph,
        branch_graph,
        branches,
        dict(source_to_branches),
        warnings,
    )
