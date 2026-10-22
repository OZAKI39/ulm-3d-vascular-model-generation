"""Bounded directed traversal; exact anatomical subsets, never geometry fitting."""
from dataclasses import dataclass

import networkx as nx

from .mevo_graph import ExactGraph, graph_counts
from .roi_landmarks import LandmarkError
from .swc_export import validate_tree


@dataclass
class ExtractedROI:
    graph: nx.DiGraph
    modelable_graph: nx.DiGraph | None
    context_nodes: set[int]
    coverage: dict
    warnings: list[str]
    anatomical_completeness: str


def downstream(graph, starts):
    seen, stack = set(), list(starts)
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(graph.successors(node))
    return seen


def _path_record(exact, graph, terminal):
    reverse = [terminal]
    while graph.in_degree(reverse[-1]):
        reverse.append(next(graph.predecessors(reverse[-1])))
    nodes = list(reversed(reverse))
    edges = []
    for pair in zip(nodes[:-1], nodes[1:]):
        edge = list(exact.edge_locations[pair])
        if not edges or edges[-1] != edge:
            edges.append(edge)
    return {"terminal_id": terminal, "node_ids": nodes, "topology_edges": edges}


def extract_graph(exact: ExactGraph, roi: dict, *, allow_natural_terminal=False,
                  proximal_context="none") -> ExtractedROI:
    source = exact.graph
    starts, stops, excludes = (set(roi[key]) for key in ("proximal_nodes", "distal_nodes", "exclude_subtree_roots"))
    if not starts or (not stops and not allow_natural_terminal):
        raise LandmarkError("Explicit proximal and distal boundaries are required")
    if roi["proximal_boundary_status"] == "unknown":
        raise LandmarkError("Proximal anatomical boundary is unknown: confirm a landmark or declare truncated_to_available_data")
    if (starts & stops) or (starts & excludes) or (stops & excludes):
        raise LandmarkError("Proximal, distal and excluded node sets must be disjoint")
    internal = {n for values in roi.get("internal_landmarks", {}).values() for n in values}
    missing = sorted((starts | stops | excludes | internal) - set(source))
    if missing:
        raise LandmarkError(f"Landmark nodes do not exist: {missing}")
    reachable = downstream(source, starts)
    if not stops <= reachable:
        raise LandmarkError(f"Configured proximal boundary is not upstream of distal boundary: {sorted(stops-reachable)}")
    if not excludes <= reachable:
        raise LandmarkError(f"Excluded subtree roots are not downstream of proximal boundaries: {sorted(excludes-reachable)}")
    visited, cut, omitted, stack = set(), set(), set(), list(starts)
    while stack:
        node = stack.pop()
        if node in excludes:
            omitted.add(node)
            continue
        if node in visited:
            continue
        visited.add(node)
        if node in stops:
            cut.add(node)
            continue
        stack.extend(source.successors(node))
    graph = source.subgraph(visited).copy()
    bypassed = [n for n in cut if graph.out_degree(n)]
    if bypassed:
        raise LandmarkError(f"Additional proximal nodes bypass distal cuts: {bypassed}")
    unused_cuts, unused_excludes = sorted(stops-cut), sorted(excludes-omitted)
    if unused_cuts or unused_excludes:
        raise LandmarkError("Configured boundaries were not visited before other cuts/exclusions",
                            context={"unvisited_distal_nodes": unused_cuts, "unvisited_exclude_roots": unused_excludes})
    if not internal <= visited:
        raise LandmarkError(f"Optional internal landmarks lie outside the strict ROI: {sorted(internal-visited)}")
    natural = sorted(n for n in graph if source.out_degree(n) == 0 and n not in stops)
    unbounded = [] if allow_natural_terminal else [_path_record(exact, graph, n) for n in natural]
    coverage = {"covered_path_count": len(cut), "excluded_path_count": len(omitted),
                "natural_terminal_count": len(natural), "unbounded_path_count": len(unbounded),
                "visited_distal_nodes": sorted(cut), "visited_exclude_roots": sorted(omitted),
                "allowed_natural_terminal_ids": natural if allow_natural_terminal else [],
                "unbounded_paths": unbounded}
    if unbounded:
        raise LandmarkError("Unbounded distal path: a source terminal was reached without an explicit distal cut",
                            context={"coverage": coverage, "candidate_counts": graph_counts(graph)})
    if len(graph) < 2 or not nx.is_weakly_connected(graph):
        raise LandmarkError("ROI must be one connected directed subgraph with at least two nodes; choose a shared proximal junction",
                            context={"coverage": coverage, "candidate_counts": graph_counts(graph)})
    root = validate_tree(graph)
    if root not in starts:
        raise LandmarkError("ROI root is not a configured proximal boundary")
    warnings = []
    if natural:
        warnings.append(f"Explicit --allow-natural-terminal accepted source terminals {natural}; their anatomical location requires user confirmation")
    completeness = "partial" if roi["proximal_boundary_status"] == "truncated_to_available_data" else "complete_as_annotated"
    if completeness == "partial":
        warnings.append("The anatomical proximal boundary is not directly represented in this source tree.")
    if not roi["manually_verified"]:
        warnings.append("ANATOMICAL_STATUS = UNVERIFIED: candidate ROI only")
    context_nodes, modelable = set(), None
    if proximal_context not in {"none", "upstream-edge"}:
        raise ValueError("proximal_context must be none or upstream-edge")
    if proximal_context == "upstream-edge":
        node = root
        while source.in_degree(node):
            node = next(source.predecessors(node))
            context_nodes.add(node)
            if source.in_degree(node) != 1 or source.out_degree(node) != 1:
                break
        modelable = source.subgraph(visited | context_nodes).copy()
        validate_tree(modelable)
        if not context_nodes:
            warnings.append("No upstream context is available in this source; no proximal geometry was invented")
    return ExtractedROI(graph, modelable, context_nodes, coverage, warnings, completeness)
