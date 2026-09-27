"""Explainable acceptance checks for the directed Schmid graph."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np

from ..reporting.acceptance import AcceptanceCheck, AcceptanceResult
from .config import SchmidPKLConfig
from .model import DirectedGraphResult


def _overall(checks: list[AcceptanceCheck]) -> str:
    statuses = {item.status for item in checks}
    if "FAIL" in statuses:
        return "FAIL"
    if "WARNING" in statuses:
        return "WARNING"
    return "PASS"


def evaluate_schmid_acceptance(
    result: DirectedGraphResult,
    config: SchmidPKLConfig,
    required_files: list[Path],
) -> AcceptanceResult:
    cleanup = result.cleanup
    checks: list[AcceptanceCheck] = []
    accounted_edges = len(cleanup.edges) + len(cleanup.removed_edge_ids)
    checks.append(
        AcceptanceCheck(
            "Every source edge is accounted for",
            "PASS" if accounted_edges == cleanup.source.edge_count else "FAIL",
            f"source={cleanup.source.edge_count}, retained={len(cleanup.edges)}, "
            f"removed={len(cleanup.removed_edge_ids)}.",
        )
    )
    represented_edges = sorted(result.raw_edge_to_branch)
    retained_edges = sorted(edge.edge_id for edge in cleanup.edges)
    checks.append(
        AcceptanceCheck(
            "Every retained edge belongs to exactly one branch",
            "PASS" if represented_edges == retained_edges else "FAIL",
            f"represented={len(represented_edges)}, retained={len(retained_edges)}.",
        )
    )
    checks.append(
        AcceptanceCheck(
            "Retained graph is connected",
            "PASS" if result.weak_component_count == 1 else "FAIL",
            f"Weak/undirected components: {result.weak_component_count}.",
        )
    )

    pressure = cleanup.source.pressure_mmhg
    wrong_direction: list[int] = []
    for branch in result.branches:
        if branch.direction_status != "known":
            continue
        if pressure[int(branch.upstream_node)] - pressure[int(branch.downstream_node)] <= config.pressure_tolerance_mmhg:
            wrong_direction.append(branch.branch_id)
    checks.append(
        AcceptanceCheck(
            "Known branches point from higher to lower pressure",
            "PASS" if not wrong_direction else "FAIL",
            "All known branches follow the pressure drop."
            if not wrong_direction
            else f"Wrong branch IDs: {wrong_direction[:20]}",
        )
    )
    unresolved = [item for item in result.branches if item.direction_status != "known"]
    checks.append(
        AcceptanceCheck(
            "Unresolved directions are explicit",
            "PASS" if not unresolved else "WARNING",
            f"Known={len(result.branches) - len(unresolved)}, unresolved={len(unresolved)}. "
            "Unresolved branches were not assigned an arbitrary arrow.",
        )
    )
    checks.append(
        AcceptanceCheck(
            "Known-direction graph contains no directed cycle",
            "PASS" if result.directed_is_acyclic else "FAIL",
            "Strict pressure descent is acyclic."
            if result.directed_is_acyclic
            else "A directed cycle was found despite pressure-based orientation.",
        )
    )
    expected_rank = (
        result.all_connectivity_graph.number_of_edges()
        - result.all_connectivity_graph.number_of_nodes()
        + result.weak_component_count
    )
    checks.append(
        AcceptanceCheck(
            "Undirected cycle basis preserves network redundancy",
            "PASS" if len(result.cycles) == expected_rank else "FAIL",
            f"cycle rank={expected_rank}, stored basis={len(result.cycles)}.",
        )
    )

    inverse_errors: list[int] = []
    branch_lookup = {item.branch_id: item for item in result.branches}
    for branch in result.branches:
        for child_id in branch.child_branch_ids:
            if branch.branch_id not in branch_lookup[child_id].parent_branch_ids:
                inverse_errors.append(branch.branch_id)
        for parent_id in branch.parent_branch_ids:
            if branch.branch_id not in branch_lookup[parent_id].child_branch_ids:
                inverse_errors.append(branch.branch_id)
    checks.append(
        AcceptanceCheck(
            "Parent and child records agree in both directions",
            "PASS" if not inverse_errors else "FAIL",
            "All branch relations are reciprocal."
            if not inverse_errors
            else f"Inconsistent branch IDs: {sorted(set(inverse_errors))[:20]}",
        )
    )

    internal = [
        row
        for row in result.flow_conservation
        if not row["is_pressure_boundary"]
        and row["unresolved_incident_edge_count"] == 0
        and row["incoming_flow_um3_per_ms"] > 0
        and row["outgoing_flow_um3_per_ms"] > 0
    ]
    relative = np.asarray([float(row["relative_imbalance"]) for row in internal], dtype=float)
    median_relative = float(np.median(relative)) if len(relative) else float("nan")
    p99_relative = float(np.quantile(relative, 0.99)) if len(relative) else float("nan")
    if not len(relative):
        balance_status = "WARNING"
    elif median_relative <= 1.0e-6 and p99_relative <= 0.10:
        balance_status = "PASS"
    elif median_relative <= 1.0e-3 and p99_relative <= 0.25:
        balance_status = "WARNING"
    else:
        balance_status = "FAIL"
    checks.append(
        AcceptanceCheck(
            "Internal flow is approximately conserved",
            balance_status,
            f"eligible nodes={len(relative)}, median relative imbalance={median_relative:.3g}, "
            f"99th percentile={p99_relative:.3g}.",
        )
    )

    degenerate = [
        edge.edge_id
        for edge in cleanup.edges
        if edge.geometry_status == "degenerate_source_geometry"
    ]
    checks.append(
        AcceptanceCheck(
            "Source centerline geometry is non-degenerate",
            "PASS" if not degenerate else "WARNING",
            "All retained edges have usable centerlines."
            if not degenerate
            else f"{len(degenerate)} source edges have zero-length point sequences; topology/flow retained.",
        )
    )
    endpoint_bad = [
        edge.edge_id
        for edge in cleanup.edges
        if edge.endpoint_error_um > config.endpoint_tolerance_um
    ]
    checks.append(
        AcceptanceCheck(
            "Centerline endpoints match graph vertices",
            "PASS" if not endpoint_bad else "WARNING",
            "All source paths align with their endpoint vertices."
            if not endpoint_bad
            else f"Endpoint mismatch on {len(endpoint_bad)} edges.",
        )
    )
    missing = [str(path) for path in required_files if not path.is_file() or path.stat().st_size == 0]
    checks.append(
        AcceptanceCheck(
            "Required output files exist",
            "PASS" if not missing else "FAIL",
            "All required files were written." if not missing else f"Missing or empty: {missing}",
        )
    )
    graphml_errors: list[str] = []
    for path in required_files:
        if path.suffix.lower() != ".graphml":
            continue
        try:
            nx.read_graphml(path)
        except Exception as exc:  # pragma: no cover - exact parser message varies
            graphml_errors.append(f"{path.name}: {exc}")
    checks.append(
        AcceptanceCheck(
            "GraphML outputs can be reopened",
            "PASS" if not graphml_errors else "FAIL",
            "All GraphML files reopened successfully."
            if not graphml_errors
            else "; ".join(graphml_errors),
        )
    )
    return AcceptanceResult(overall_status=_overall(checks), checks=checks)
