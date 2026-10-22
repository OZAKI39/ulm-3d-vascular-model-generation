"""Acceptance checks for parent-to-current graph construction."""

from __future__ import annotations

from pathlib import Path
import networkx as nx

from ..reporting.acceptance import AcceptanceCheck, AcceptanceResult
from .model import DirectedVascularGraph


def _overall(checks: list[AcceptanceCheck]) -> str:
    statuses = {check.status for check in checks}
    return "FAIL" if "FAIL" in statuses else ("WARNING" if "WARNING" in statuses else "PASS")


def evaluate_directed_graph(
    result: DirectedVascularGraph,
    required_files: list[Path],
    *,
    strict_nonpositive_radius: bool,
) -> AcceptanceResult:
    checks: list[AcceptanceCheck] = []
    checks.append(
        AcceptanceCheck(
            "SWC structural integrity",
            "PASS" if result.swc.structurally_valid else "FAIL",
            f"Validation: {result.swc.validation}",
        )
    )
    represented: list[tuple[int, int]] = []
    reversed_or_missing: list[tuple[int, int]] = []
    for branch in result.branches:
        for edge in zip(branch.source_node_ids[:-1], branch.source_node_ids[1:]):
            represented.append(edge)
            if not result.source_graph.has_edge(*edge):
                reversed_or_missing.append(edge)
    source_edges = set(result.source_graph.edges)
    represented_edges = set(represented)
    exact = (
        source_edges == represented_edges
        and len(represented) == len(represented_edges)
        and not reversed_or_missing
    )
    checks.append(
        AcceptanceCheck(
            "Every SWC parent edge is represented once and forward",
            "PASS" if exact else "FAIL",
            f"Source={len(source_edges)}, represented={len(represented)}, "
            f"unique={len(represented_edges)}, invalid={reversed_or_missing[:10]}.",
        )
    )
    endpoint_valid = all(
        branch.source_node_ids[0] == branch.upstream_node_id
        and branch.source_node_ids[-1] == branch.downstream_node_id
        for branch in result.branches
    )
    checks.append(
        AcceptanceCheck(
            "Branch endpoints preserve upstream-to-downstream order",
            "PASS" if endpoint_valid else "FAIL",
            "All branch point arrays start at parent-side and end at child-side."
            if endpoint_valid
            else "At least one branch endpoint disagrees with its source-node sequence.",
        )
    )
    dag = nx.is_directed_acyclic_graph(result.branch_graph)
    checks.append(
        AcceptanceCheck(
            "Directed branch hierarchy is acyclic",
            "PASS" if dag else "FAIL",
            f"Branches={len(result.branches)}, hierarchy edges={result.branch_graph.number_of_edges()}.",
        )
    )
    roots = sum(result.source_graph.in_degree(node) == 0 for node in result.source_graph)
    leaves = sum(result.source_graph.out_degree(node) == 0 for node in result.source_graph)
    checks.append(
        AcceptanceCheck(
            "Inlet/outlet candidates follow topology",
            "PASS" if roots == len(result.swc.root_ids) and leaves > 0 else "FAIL",
            f"Inferred inlet/root candidates={roots}; inferred outlet/leaf candidates={leaves}. "
            "These are structural labels, not measured hemodynamic boundary conditions.",
        )
    )
    invalid_radius = result.swc.validation["nonpositive_radius_node_ids"]
    checks.append(
        AcceptanceCheck(
            "SWC radii are positive and finite",
            ("FAIL" if strict_nonpositive_radius else "WARNING") if invalid_radius else "PASS",
            f"Invalid radius node count={len(invalid_radius)}; raw values are never overwritten.",
        )
    )
    missing = [str(path) for path in required_files if not path.is_file() or path.stat().st_size == 0]
    checks.append(
        AcceptanceCheck(
            "Required graph and visualization artifacts exist",
            "PASS" if not missing else "FAIL",
            "All required artifacts were written." if not missing else f"Missing or empty: {missing}",
        )
    )
    return AcceptanceResult(_overall(checks), checks)
