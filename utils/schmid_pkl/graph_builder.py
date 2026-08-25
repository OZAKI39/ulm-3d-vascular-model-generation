"""Step 3: build a directed hierarchical multigraph and parent/child relations."""

from __future__ import annotations

import logging
from collections import defaultdict

import networkx as nx
import numpy as np

from .config import SchmidPKLConfig
from .geometry import concatenate_edge_geometry, measure_branch_geometry
from .model import (
    CleanEdge,
    DirectedBranch,
    DirectedCycle,
    DirectedGraphResult,
    DirectedNode,
    SchmidCleanupResult,
    VESSEL_TYPE_NAMES,
    optional_float,
)


def _assign_edge_directions(cleanup: SchmidCleanupResult, config: SchmidPKLConfig) -> None:
    pressure = cleanup.source.pressure_mmhg
    for edge in cleanup.edges:
        difference = float(pressure[edge.node_u] - pressure[edge.node_v])
        if difference > config.pressure_tolerance_mmhg:
            edge.upstream_node = edge.node_u
            edge.downstream_node = edge.node_v
            edge.pressure_drop_mmhg = difference
            edge.direction_status = "known"
        elif difference < -config.pressure_tolerance_mmhg:
            edge.upstream_node = edge.node_v
            edge.downstream_node = edge.node_u
            edge.pressure_drop_mmhg = -difference
            edge.direction_status = "known"
        else:
            edge.upstream_node = None
            edge.downstream_node = None
            edge.pressure_drop_mmhg = None
            edge.direction_status = (
                "unresolved_equal_pressure_zero_flow"
                if edge.flow_um3_per_ms <= config.zero_flow_tolerance_um3_per_ms
                else "unresolved_equal_pressure_nonzero_flow"
            )


def _raw_graphs(
    cleanup: SchmidCleanupResult,
) -> tuple[nx.MultiGraph, nx.MultiDiGraph, dict[int, CleanEdge]]:
    source = cleanup.source
    undirected = nx.MultiGraph(
        representation="raw_schmid_connectivity",
        directed_for_flow=False,
        coordinate_system="source_XYZ_not_anatomically_labeled",
    )
    directed = nx.MultiDiGraph(
        representation="raw_schmid_pressure_directed",
        directed_for_flow=True,
        direction_method="higher endpoint pressure to lower endpoint pressure",
    )
    for node_id in cleanup.valid_node_ids:
        node_id = int(node_id)
        attributes = {
            "pressure_mmhg": float(source.pressure_mmhg[node_id]),
            "x_um": float(source.coordinates_um[node_id, 0]),
            "y_um": float(source.coordinates_um[node_id, 1]),
            "z_um": float(source.coordinates_um[node_id, 2]),
        }
        undirected.add_node(node_id, **attributes)
        directed.add_node(node_id, **attributes)
    lookup = {edge.edge_id: edge for edge in cleanup.edges}
    for edge in cleanup.edges:
        attributes = {
            "edge_id": edge.edge_id,
            "flow_um3_per_ms": edge.flow_um3_per_ms,
            "source_length_um": edge.source_length_um,
            "mean_diameter_um": edge.mean_diameter_um,
            "vessel_type_code": edge.vessel_type_code,
            "vessel_type": VESSEL_TYPE_NAMES.get(edge.vessel_type_code, "unrecognized"),
            "direction_status": edge.direction_status,
        }
        undirected.add_edge(edge.node_u, edge.node_v, key=edge.edge_id, **attributes)
        if edge.direction_status == "known":
            directed.add_edge(
                edge.upstream_node,
                edge.downstream_node,
                key=edge.edge_id,
                pressure_drop_mmhg=edge.pressure_drop_mmhg,
                **attributes,
            )
    return undirected, directed, lookup


def _critical_nodes(
    graph: nx.MultiGraph,
    directed: nx.MultiDiGraph,
    edges: dict[int, CleanEdge],
    pressure_boundary: np.ndarray,
) -> set[int]:
    critical: set[int] = set()
    for node_id in graph.nodes:
        incident_ids = [int(key) for _, _, key in graph.edges(node_id, keys=True)]
        vessel_types = {edges[edge_id].vessel_type_code for edge_id in incident_ids}
        unresolved = any(edges[edge_id].direction_status != "known" for edge_id in incident_ids)
        if (
            graph.degree(node_id) != 2
            or directed.in_degree(node_id) != 1
            or directed.out_degree(node_id) != 1
            or unresolved
            or len(vessel_types) > 1
            or np.isfinite(pressure_boundary[node_id])
        ):
            critical.add(int(node_id))
    return critical


def _incident_edges(graph: nx.MultiGraph, node_id: int) -> list[int]:
    return sorted(int(key) for _, _, key in graph.edges(node_id, keys=True))


def _other_endpoint(edge: CleanEdge, node_id: int) -> int:
    if edge.node_u == node_id:
        return edge.node_v
    if edge.node_v == node_id:
        return edge.node_u
    raise ValueError(f"Edge {edge.edge_id} is not incident to node {node_id}")


def _trace_from_edge(
    start_node: int,
    first_edge_id: int,
    graph: nx.MultiGraph,
    edges: dict[int, CleanEdge],
    critical: set[int],
    visited: set[int],
) -> tuple[list[int], list[int]]:
    node_path = [start_node]
    edge_path: list[int] = []
    current_node = start_node
    current_edge_id = first_edge_id
    while True:
        if current_edge_id in visited:
            break
        visited.add(current_edge_id)
        edge_path.append(current_edge_id)
        next_node = _other_endpoint(edges[current_edge_id], current_node)
        node_path.append(next_node)
        if next_node in critical and (next_node != start_node or len(edge_path) > 1):
            break
        candidates = [
            edge_id
            for edge_id in _incident_edges(graph, next_node)
            if edge_id != current_edge_id and edge_id not in visited
        ]
        if not candidates:
            break
        current_node = next_node
        current_edge_id = candidates[0]
    return node_path, edge_path


def _branch_paths(
    graph: nx.MultiGraph,
    edges: dict[int, CleanEdge],
    critical: set[int],
) -> list[tuple[list[int], list[int]]]:
    visited: set[int] = set()
    paths: list[tuple[list[int], list[int]]] = []
    for node_id in sorted(critical):
        for edge_id in _incident_edges(graph, node_id):
            if edge_id not in visited:
                paths.append(
                    _trace_from_edge(node_id, edge_id, graph, edges, critical, visited)
                )
    for edge_id in sorted(edges):
        if edge_id in visited:
            continue
        anchor = min(edges[edge_id].node_u, edges[edge_id].node_v)
        critical.add(anchor)
        paths.append(_trace_from_edge(anchor, edge_id, graph, edges, critical, visited))
    return paths


def _orient_path(
    node_path: list[int], edge_path: list[int], edges: dict[int, CleanEdge]
) -> tuple[list[int], list[int], str]:
    signs: list[int] = []
    for left, right, edge_id in zip(node_path[:-1], node_path[1:], edge_path, strict=True):
        edge = edges[edge_id]
        if edge.direction_status != "known":
            signs.append(0)
        elif edge.upstream_node == left and edge.downstream_node == right:
            signs.append(1)
        elif edge.upstream_node == right and edge.downstream_node == left:
            signs.append(-1)
        else:
            signs.append(0)
    if signs and all(value == 1 for value in signs):
        return node_path, edge_path, "known"
    if signs and all(value == -1 for value in signs):
        return list(reversed(node_path)), list(reversed(edge_path)), "known"
    if any(value == 0 for value in signs):
        return node_path, edge_path, "unresolved_source_edge"
    return node_path, edge_path, "unresolved_inconsistent_chain"


def _edge_geometry_for_step(
    edge: CleanEdge, start_node: int, end_node: int
) -> tuple[np.ndarray, np.ndarray]:
    if edge.node_u == start_node and edge.node_v == end_node:
        return edge.points_u_to_v_um, edge.diameter_u_to_v_um
    if edge.node_v == start_node and edge.node_u == end_node:
        return edge.points_u_to_v_um[::-1], edge.diameter_u_to_v_um[::-1]
    raise ValueError(f"Edge {edge.edge_id} does not join {start_node} and {end_node}")


def _make_branch(
    branch_id: int,
    node_path: list[int],
    edge_path: list[int],
    edges: dict[int, CleanEdge],
    config: SchmidPKLConfig,
    pressure: np.ndarray,
) -> DirectedBranch:
    node_path, edge_path, direction_status = _orient_path(node_path, edge_path, edges)
    pieces = [
        _edge_geometry_for_step(edges[edge_id], left, right)
        for left, right, edge_id in zip(node_path[:-1], node_path[1:], edge_path, strict=True)
    ]
    points_raw, radius_raw = concatenate_edge_geometry(pieces)
    measurement = measure_branch_geometry(
        points_raw,
        radius_raw,
        resample_step_um=config.resample_step_um,
        smoothing_enabled=config.smoothing_enabled,
        smoothing_window_points=config.smoothing_window_points,
    )
    source_edges = [edges[edge_id] for edge_id in edge_path]
    statuses = {edge.geometry_status for edge in source_edges}
    if statuses == {"valid"}:
        geometry_status = "valid"
    elif len(source_edges) == 1:
        geometry_status = source_edges[0].geometry_status
    else:
        geometry_status = "mixed_source_geometry_issues"
    flow_values = np.asarray([edge.flow_um3_per_ms for edge in source_edges], dtype=float)
    upstream = node_path[0] if direction_status == "known" else None
    downstream = node_path[-1] if direction_status == "known" else None
    pressure_drop = (
        float(pressure[upstream] - pressure[downstream])
        if upstream is not None and downstream is not None
        else None
    )
    return DirectedBranch(
        branch_id=branch_id,
        raw_edge_ids=edge_path,
        node_a=node_path[0],
        node_b=node_path[-1],
        upstream_node=upstream,
        downstream_node=downstream,
        direction_status=direction_status,
        points_raw_um=points_raw,
        radius_raw_um=radius_raw,
        points_smoothed_um=measurement["points_smoothed_um"],
        radius_smoothed_um=measurement["radius_smoothed_um"],
        arc_length_raw_um=measurement["arc_length_raw_um"],
        arc_length_smoothed_um=measurement["arc_length_smoothed_um"],
        local_direction_smoothed=measurement["local_direction_smoothed"],
        curvature_smoothed_per_um=measurement["curvature_smoothed_per_um"],
        source_length_um=float(sum(edge.source_length_um for edge in source_edges)),
        geometric_length_um=float(measurement["geometric_length_um"]),
        straight_distance_um=float(measurement["straight_distance_um"]),
        tortuosity_raw=measurement["tortuosity_raw"],
        tortuosity_smoothed=measurement["tortuosity_smoothed"],
        smoothing_max_deviation_um=float(measurement["smoothing_max_deviation_um"]),
        flow_um3_per_ms=float(np.median(flow_values)),
        flow_min_um3_per_ms=float(flow_values.min()),
        flow_max_um3_per_ms=float(flow_values.max()),
        pressure_drop_mmhg=pressure_drop,
        vessel_type_codes=sorted({edge.vessel_type_code for edge in source_edges}),
        geometry_status=geometry_status,
    )


def _node_role(node: DirectedNode) -> str:
    incoming = len(node.incoming_branch_ids)
    outgoing = len(node.outgoing_branch_ids)
    unresolved = len(node.unresolved_branch_ids)
    if unresolved:
        return "unresolved_junction" if node.undirected_degree > 1 else "unresolved_terminal"
    if incoming == 0 and outgoing == 0:
        return "isolated"
    if incoming == 0:
        return "source" if outgoing == 1 else "source_split"
    if outgoing == 0:
        return "sink" if incoming == 1 else "merge_sink"
    if incoming == 1 and outgoing == 1:
        return "connector"
    if incoming == 1 and outgoing > 1:
        return "split"
    if incoming > 1 and outgoing == 1:
        return "merge"
    return "mixed_split_merge"


def _build_aggregated_graphs(
    branches: list[DirectedBranch], nodes: list[DirectedNode]
) -> tuple[nx.MultiGraph, nx.MultiDiGraph, nx.MultiDiGraph]:
    all_graph = nx.MultiGraph(
        representation="junction_as_node_all_connectivity",
        directed_for_flow=False,
    )
    directed_graph = nx.MultiDiGraph(
        representation="junction_as_node_pressure_directed",
        directed_for_flow=True,
    )
    for node in nodes:
        attributes = {
            "node_role": node.node_role,
            "pressure_mmhg": node.pressure_mmhg,
            "x_um": node.coordinates_um[0],
            "y_um": node.coordinates_um[1],
            "z_um": node.coordinates_um[2],
        }
        all_graph.add_node(node.node_id, **attributes)
        directed_graph.add_node(node.node_id, **attributes)
    for branch in branches:
        attributes = {
            "branch_id": branch.branch_id,
            "direction_status": branch.direction_status,
            "source_length_um": branch.source_length_um,
            "geometric_length_um": branch.geometric_length_um,
            "flow_um3_per_ms": branch.flow_um3_per_ms,
            "pressure_drop_mmhg": branch.pressure_drop_mmhg or 0.0,
            "mean_radius_um": float(np.mean(branch.radius_raw_um)),
            "geometry_status": branch.geometry_status,
            "vessel_type_codes": ",".join(map(str, branch.vessel_type_codes)),
        }
        all_graph.add_edge(branch.node_a, branch.node_b, key=branch.branch_id, **attributes)
        if branch.direction_status == "known":
            directed_graph.add_edge(
                branch.upstream_node,
                branch.downstream_node,
                key=branch.branch_id,
                **attributes,
            )

    branch_graph = nx.MultiDiGraph(
        representation="branch_as_node_pressure_directed",
        relationship="parent_branch_to_child_branch",
    )
    for branch in branches:
        branch_graph.add_node(
            branch.branch_id,
            direction_status=branch.direction_status,
            source_length_um=branch.source_length_um,
            flow_um3_per_ms=branch.flow_um3_per_ms,
            geometry_status=branch.geometry_status,
        )
    for node in nodes:
        for parent_id in node.incoming_branch_ids:
            parent = branches[parent_id]
            for child_id in node.outgoing_branch_ids:
                if parent_id == child_id:
                    continue
                child = branches[child_id]
                parent.child_branch_ids.append(child_id)
                child.parent_branch_ids.append(parent_id)
                branch_graph.add_edge(
                    parent_id,
                    child_id,
                    key=f"{node.node_id}_{parent_id}_{child_id}",
                    shared_node_id=node.node_id,
                    relationship="parent_to_child",
                )
    for branch in branches:
        branch.parent_branch_ids = sorted(set(branch.parent_branch_ids))
        branch.child_branch_ids = sorted(set(branch.child_branch_ids))
    return all_graph, directed_graph, branch_graph


def _cycle_basis(graph: nx.MultiGraph) -> list[DirectedCycle]:
    simple = nx.Graph()
    simple.add_nodes_from(graph.nodes)
    pair_to_branches: dict[tuple[int, int], list[int]] = defaultdict(list)
    self_loops: list[int] = []
    for left, right, key in graph.edges(keys=True):
        branch_id = int(key)
        if left == right:
            self_loops.append(branch_id)
            continue
        pair = tuple(sorted((int(left), int(right))))
        pair_to_branches[pair].append(branch_id)
        simple.add_edge(*pair)
    records: list[DirectedCycle] = []
    for node_ids in nx.cycle_basis(simple):
        branch_ids = []
        wrapped = node_ids[1:] + node_ids[:1]
        for left, right in zip(node_ids, wrapped, strict=True):
            branch_ids.append(min(pair_to_branches[tuple(sorted((left, right)))]))
        records.append(
            DirectedCycle(len(records), "simple_cycle_basis", list(map(int, node_ids)), branch_ids)
        )
    for pair, branch_ids in sorted(pair_to_branches.items()):
        ordered = sorted(branch_ids)
        for branch_id in ordered[1:]:
            records.append(
                DirectedCycle(len(records), "parallel_branch_cycle", list(pair), [ordered[0], branch_id])
            )
    for branch_id in sorted(self_loops):
        left, _, _ = next((u, v, k) for u, v, k in graph.edges(keys=True) if int(k) == branch_id)
        records.append(DirectedCycle(len(records), "self_loop_cycle", [int(left)], [branch_id]))
    return records


def _flow_conservation(
    cleanup: SchmidCleanupResult, edges: dict[int, CleanEdge]
) -> list[dict[str, object]]:
    incoming: dict[int, float] = defaultdict(float)
    outgoing: dict[int, float] = defaultdict(float)
    unresolved: dict[int, int] = defaultdict(int)
    for edge in edges.values():
        if edge.direction_status == "known":
            outgoing[int(edge.upstream_node)] += edge.flow_um3_per_ms
            incoming[int(edge.downstream_node)] += edge.flow_um3_per_ms
        else:
            unresolved[edge.node_u] += 1
            unresolved[edge.node_v] += 1
    output: list[dict[str, object]] = []
    pbc = cleanup.source.pressure_boundary_mmhg
    for node_id in cleanup.valid_node_ids:
        node_id = int(node_id)
        flow_in = incoming[node_id]
        flow_out = outgoing[node_id]
        absolute = abs(flow_in - flow_out)
        relative = absolute / max(flow_in, flow_out, np.finfo(float).eps)
        output.append(
            {
                "node_id": node_id,
                "is_pressure_boundary": bool(np.isfinite(pbc[node_id])),
                "incoming_flow_um3_per_ms": flow_in,
                "outgoing_flow_um3_per_ms": flow_out,
                "absolute_imbalance_um3_per_ms": absolute,
                "relative_imbalance": relative,
                "unresolved_incident_edge_count": unresolved[node_id],
            }
        )
    return output


def build_directed_hierarchical_graph(
    cleanup: SchmidCleanupResult,
    config: SchmidPKLConfig,
    logger: logging.Logger | None = None,
) -> DirectedGraphResult:
    _assign_edge_directions(cleanup, config)
    raw_undirected, raw_directed, edge_lookup = _raw_graphs(cleanup)
    critical = _critical_nodes(
        raw_undirected,
        raw_directed,
        edge_lookup,
        cleanup.source.pressure_boundary_mmhg,
    )
    paths = _branch_paths(raw_undirected, edge_lookup, critical)
    branches = [
        _make_branch(
            branch_id,
            node_path,
            edge_path,
            edge_lookup,
            config,
            cleanup.source.pressure_mmhg,
        )
        for branch_id, (node_path, edge_path) in enumerate(paths)
    ]
    raw_edge_to_branch = {
        edge_id: branch.branch_id for branch in branches for edge_id in branch.raw_edge_ids
    }
    endpoint_ids = sorted({value for branch in branches for value in (branch.node_a, branch.node_b)})
    source = cleanup.source
    nodes = [
        DirectedNode(
            node_id=node_id,
            coordinates_um=tuple(float(value) for value in source.coordinates_um[node_id]),
            pressure_mmhg=float(source.pressure_mmhg[node_id]),
            pressure_boundary_mmhg=optional_float(float(source.pressure_boundary_mmhg[node_id])),
        )
        for node_id in endpoint_ids
    ]
    node_lookup = {node.node_id: node for node in nodes}
    for branch in branches:
        node_lookup[branch.node_a].undirected_degree += 1
        node_lookup[branch.node_b].undirected_degree += 1
        if branch.direction_status == "known":
            upstream = node_lookup[int(branch.upstream_node)]
            downstream = node_lookup[int(branch.downstream_node)]
            upstream.outgoing_branch_ids.append(branch.branch_id)
            upstream.child_node_ids.append(downstream.node_id)
            downstream.incoming_branch_ids.append(branch.branch_id)
            downstream.parent_node_ids.append(upstream.node_id)
        else:
            node_lookup[branch.node_a].unresolved_branch_ids.append(branch.branch_id)
            node_lookup[branch.node_b].unresolved_branch_ids.append(branch.branch_id)
    for node in nodes:
        node.incoming_branch_ids.sort()
        node.outgoing_branch_ids.sort()
        node.unresolved_branch_ids.sort()
        node.parent_node_ids = sorted(set(node.parent_node_ids))
        node.child_node_ids = sorted(set(node.child_node_ids))
        node.node_role = _node_role(node)

    all_graph, directed_graph, branch_graph = _build_aggregated_graphs(branches, nodes)
    cycles = _cycle_basis(all_graph)
    for cycle in cycles:
        for branch_id in cycle.branch_ids:
            branches[branch_id].cycle_ids.append(cycle.cycle_id)
        for node_id in cycle.node_ids:
            node_lookup[node_id].cycle_ids.append(cycle.cycle_id)
    raw_unresolved = sorted(
        edge.edge_id for edge in cleanup.edges if edge.direction_status != "known"
    )
    weak_components = nx.number_connected_components(all_graph)
    directed_is_acyclic = nx.is_directed_acyclic_graph(directed_graph)
    result = DirectedGraphResult(
        nodes=nodes,
        branches=branches,
        cycles=cycles,
        all_connectivity_graph=all_graph,
        directed_junction_graph=directed_graph,
        directed_branch_graph=branch_graph,
        raw_directed_graph=raw_directed,
        raw_unresolved_edge_ids=raw_unresolved,
        raw_edge_to_branch=raw_edge_to_branch,
        weak_component_count=weak_components,
        directed_is_acyclic=directed_is_acyclic,
        cleanup=cleanup,
        flow_conservation=_flow_conservation(cleanup, edge_lookup),
    )
    if logger is not None:
        report = result.report()
        logger.info(
            "Step 3 graph: %d endpoint/junction nodes, %d branches, %d stored cycles",
            report["node_count"],
            report["branch_count"],
            report["stored_cycle_basis_count"],
        )
        logger.info(
            "Directions: known=%d, unresolved=%d; known-direction graph is acyclic=%s",
            report["known_direction_branch_count"],
            report["unresolved_direction_branch_count"],
            report["known_direction_graph_is_acyclic"],
        )
    return result
