"""Explainable acceptance checks for hierarchical vascular graphs."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..config import HierarchicalGraphConfig
from ..reporting.acceptance import AcceptanceCheck, AcceptanceResult
from .model import HierarchicalGraphResult


def _overall(checks: list[AcceptanceCheck]) -> str:
    statuses = {item.status for item in checks}
    if "FAIL" in statuses:
        return "FAIL"
    if "WARNING" in statuses:
        return "WARNING"
    return "PASS"


def _check(
    checks: list[AcceptanceCheck], name: str, passed: bool, message: str
) -> None:
    checks.append(AcceptanceCheck(name, "PASS" if passed else "FAIL", message))


def evaluate_hierarchical_graph_acceptance(
    result: HierarchicalGraphResult,
    config: HierarchicalGraphConfig,
    required_files: list[Path],
    export_errors: list[str],
) -> AcceptanceResult:
    checks: list[AcceptanceCheck] = []
    _check(
        checks,
        "Source skeleton is one connected network",
        result.skeleton_component_count == 1,
        f"Detected {result.skeleton_component_count} connected component(s).",
    )
    _check(
        checks,
        "Hierarchical graph is non-empty",
        bool(result.nodes) and bool(result.branches),
        f"Nodes={len(result.nodes):,}, branches={len(result.branches):,}.",
    )
    _check(
        checks,
        "Hierarchical graph remains connected",
        result.graph_component_count == result.skeleton_component_count,
        f"Skeleton components={result.skeleton_component_count}, "
        f"graph components={result.graph_component_count}.",
    )
    _check(
        checks,
        "Graph round-trip preserves every skeleton voxel",
        result.missing_voxel_count == 0 and result.extra_voxel_count == 0,
        f"Missing={result.missing_voxel_count:,}, extra={result.extra_voxel_count:,}.",
    )
    _check(
        checks,
        "Branch interiors are not duplicated",
        result.duplicate_interior_voxel_count == 0,
        f"Duplicated non-node voxels={result.duplicate_interior_voxel_count:,}.",
    )

    node_ids = {item.node_id for item in result.nodes}
    invalid_endpoints = [
        item.branch_id
        for item in result.branches
        if item.node_u not in node_ids or item.node_v not in node_ids
    ]
    _check(
        checks,
        "Every branch references valid endpoint nodes",
        not invalid_endpoints,
        "All endpoints are valid."
        if not invalid_endpoints
        else f"Invalid branch IDs: {invalid_endpoints[:20]}",
    )

    limit = {6: 1, 18: 2, 26: 3}[config.neighbor_connectivity]
    discontinuous: list[int] = []
    nonpositive: list[int] = []
    misaligned: list[int] = []
    nonfinite: list[int] = []
    for branch in result.branches:
        differences = np.abs(np.diff(branch.voxel_indices_xyz, axis=0))
        if len(differences) and (
            np.any(np.max(differences, axis=1) > 1)
            or np.any(np.sum(differences, axis=1) > limit)
            or np.any(np.sum(differences, axis=1) == 0)
        ):
            discontinuous.append(branch.branch_id)
        if branch.arc_length_raw_um[-1] <= 0:
            nonpositive.append(branch.branch_id)
        if not (
            len(branch.voxel_indices_xyz)
            == len(branch.points_raw_lps_um)
            == len(branch.arc_length_raw_um)
            == len(branch.coarse_radius_raw_um)
        ) or not (
            len(branch.points_smoothed_lps_um)
            == len(branch.arc_length_smoothed_um)
            == len(branch.coarse_radius_smoothed_um)
            == len(branch.local_direction_smoothed)
            == len(branch.curvature_smoothed_per_um)
        ):
            misaligned.append(branch.branch_id)
        arrays = (
            branch.points_raw_lps_um,
            branch.arc_length_raw_um,
            branch.coarse_radius_raw_um,
            branch.points_smoothed_lps_um,
            branch.arc_length_smoothed_um,
            branch.coarse_radius_smoothed_um,
            branch.local_direction_smoothed,
            branch.curvature_smoothed_per_um,
        )
        if any(not np.all(np.isfinite(array)) for array in arrays):
            nonfinite.append(branch.branch_id)
    _check(
        checks,
        "Every branch is a continuous voxel path",
        not discontinuous,
        "All paths are continuous."
        if not discontinuous
        else f"Discontinuous branch IDs: {discontinuous[:20]}",
    )
    _check(
        checks,
        "Every branch has positive length",
        not nonpositive,
        "All branch lengths are positive."
        if not nonpositive
        else f"Non-positive branch IDs: {nonpositive[:20]}",
    )
    _check(
        checks,
        "Dense geometry sequences are aligned",
        not misaligned,
        "Coordinates, distances, radii, directions, and curvature have matching lengths."
        if not misaligned
        else f"Misaligned branch IDs: {misaligned[:20]}",
    )
    _check(
        checks,
        "Dense geometry contains finite values",
        not nonfinite,
        "All dense arrays are finite."
        if not nonfinite
        else f"Non-finite branch IDs: {nonfinite[:20]}",
    )

    _check(
        checks,
        "Cycle identity survives degree-2 aggregation",
        len(result.cycles) == result.cycle_rank,
        f"Cycle rank={result.cycle_rank}, stored independent cycles={len(result.cycles)}.",
    )
    _check(
        checks,
        "Branch-as-node representation matches the branch table",
        result.branch_as_node_graph.number_of_nodes() == len(result.branches),
        f"Branch records={len(result.branches)}, branch-as-node nodes="
        f"{result.branch_as_node_graph.number_of_nodes()}.",
    )

    lower = np.asarray(result.origin_lps_um, dtype=float) - 1.0e-6
    upper = lower + np.asarray(result.spacing_um) * (
        np.asarray(result.source_skeleton.shape) - 1
    ) + 2.0e-6
    outside = 0
    for branch in result.branches:
        outside += int(
            np.count_nonzero(
                np.any(
                    (branch.points_raw_lps_um < lower)
                    | (branch.points_raw_lps_um > upper),
                    axis=1,
                )
            )
        )
    _check(
        checks,
        "All raw centerline coordinates stay inside the source grid",
        outside == 0,
        f"Outside raw points={outside:,}.",
    )

    type_mismatches = [
        node.node_id
        for node in result.nodes
        if (node.node_type == "terminal" and node.graph_degree != 1)
        or (node.node_type == "junction" and node.graph_degree != 3)
        or (node.node_type == "complex_junction" and node.graph_degree < 4)
    ]
    _check(
        checks,
        "Node labels agree with graph degree",
        not type_mismatches,
        "All terminal and junction labels agree with graph degree."
        if not type_mismatches
        else f"Mismatched node IDs: {type_mismatches[:20]}",
    )

    missing_files = [
        str(path) for path in required_files if not path.is_file() or path.stat().st_size == 0
    ]
    _check(
        checks,
        "Required Step 3 files exist",
        not missing_files,
        "All required files were written."
        if not missing_files
        else f"Missing or empty files: {missing_files}",
    )
    _check(
        checks,
        "Portable graph exports can be reopened",
        not export_errors,
        "JSON, NPZ, GraphML, and VTP exports reopened successfully."
        if not export_errors
        else "; ".join(export_errors),
    )

    short = [
        item.branch_id
        for item in result.branches
        if item.arc_length_raw_um[-1] < config.short_branch_warning_um
    ]
    checks.append(
        AcceptanceCheck(
            "Very short branches require visual review",
            "WARNING" if short else "PASS",
            f"Detected {len(short)} branch(es) shorter than "
            f"{config.short_branch_warning_um:g} um; none were deleted."
            if short
            else "No branch falls below the configured warning length.",
        )
    )
    large_nodes = [
        item.node_id
        for item in result.nodes
        if len(item.voxel_indices_xyz) >= config.large_junction_warning_voxels
    ]
    checks.append(
        AcceptanceCheck(
            "Large junction regions require visual review",
            "WARNING" if large_nodes else "PASS",
            f"Detected {len(large_nodes)} node region(s) with at least "
            f"{config.large_junction_warning_voxels} voxels; possible nearby-branch fusion."
            if large_nodes
            else "No unusually large junction region was detected.",
        )
    )
    high_degree = [
        item.node_id
        for item in result.nodes
        if item.graph_degree >= config.high_degree_warning
    ]
    checks.append(
        AcceptanceCheck(
            "High-degree nodes require visual review",
            "WARNING" if high_degree else "PASS",
            f"Detected {len(high_degree)} node(s) with degree >= "
            f"{config.high_degree_warning}; none were altered."
            if high_degree
            else "No high-degree node was detected.",
        )
    )
    excessive_smoothing = [
        item.branch_id
        for item in result.branches
        if item.smoothing_max_deviation_um > config.smoothing_deviation_warning_um
    ]
    checks.append(
        AcceptanceCheck(
            "Smoothing stays close to raw geometry",
            "WARNING" if excessive_smoothing else "PASS",
            f"Detected {len(excessive_smoothing)} branch(es) exceeding "
            f"{config.smoothing_deviation_warning_um:g} um; raw geometry remains preserved."
            if excessive_smoothing
            else "All smoothed sequences stay within the configured deviation.",
        )
    )
    checks.append(
        AcceptanceCheck(
            "Geometry quality is explicitly marked as coarse",
            "PASS",
            "Radius and curvature are navigation estimates, not CFD or final training truth.",
        )
    )
    checks.append(
        AcceptanceCheck(
            "Flow direction is not invented",
            "PASS",
            "The graph is undirected; parent, daughter, depth, and downstream subtree are unavailable.",
        )
    )
    return AcceptanceResult(_overall(checks), checks)
