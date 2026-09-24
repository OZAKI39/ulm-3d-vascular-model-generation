"""Acceptance checks for NNE2 stack and hierarchy outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np

from ..graph.model import HierarchicalGraphResult
from .model import NNE2HierarchyResult
from .segmentation import SegmentationResult


def validate_preprocess_result(
    segmentation: SegmentationResult,
    skeleton_zyx: np.ndarray,
    radius_zyx_um: np.ndarray,
    *,
    island_warning_fraction: float,
    island_fail_fraction: float,
) -> dict[str, Any]:
    candidate_count = int(np.count_nonzero(segmentation.candidate_mask_zyx))
    removed_fraction = segmentation.removed_voxel_count / max(1, candidate_count)
    checks = {
        "candidate_mask_nonempty": bool(np.any(segmentation.candidate_mask_zyx)),
        "cleaned_mask_nonempty": bool(np.any(segmentation.mask_zyx)),
        "foreground_fraction_below_35_percent": float(np.mean(segmentation.mask_zyx)) < 0.35,
        "kept_and_removed_masks_do_not_overlap": not bool(
            np.any(segmentation.mask_zyx & segmentation.removed_mask_zyx)
        ),
        "component_decisions_cover_candidate_mask": (
            int(sum(item["voxel_count"] for item in segmentation.component_decisions))
            == candidate_count
        ),
        "skeleton_nonempty": bool(np.any(skeleton_zyx)),
        "skeleton_inside_mask": not bool(np.any(skeleton_zyx & ~segmentation.mask_zyx)),
        "radius_shape_matches_mask": radius_zyx_um.shape == segmentation.mask_zyx.shape,
        "radius_is_finite_and_nonnegative": bool(
            np.all(np.isfinite(radius_zyx_um)) and np.all(radius_zyx_um >= 0)
        ),
    }
    island_status = (
        "FAIL" if removed_fraction >= island_fail_fraction
        else "WARN" if removed_fraction >= island_warning_fraction
        else "PASS"
    )
    status = "FAIL" if not all(checks.values()) or island_status == "FAIL" else (
        "WARN" if island_status == "WARN" else "PASS"
    )
    return {
        "status": status,
        "checks": checks,
        "removed_component_voxel_fraction": removed_fraction,
        "removed_component_safety_status": island_status,
        "warning_threshold": island_warning_fraction,
        "failure_threshold": island_fail_fraction,
    }


def validate_stack_result(
    segmentation: SegmentationResult,
    skeleton_zyx: np.ndarray,
    graph: HierarchicalGraphResult,
) -> dict[str, Any]:
    checks = {
        "mask_nonempty": bool(np.any(segmentation.mask_zyx)),
        "foreground_fraction_below_35_percent": float(np.mean(segmentation.mask_zyx)) < 0.35,
        "skeleton_nonempty": bool(np.any(skeleton_zyx)),
        "skeleton_inside_mask": not bool(np.any(skeleton_zyx & ~segmentation.mask_zyx)),
        "graph_has_branches": len(graph.branches) > 0,
        "graph_round_trip_has_no_missing_skeleton_voxels": graph.missing_voxel_count == 0,
        "graph_round_trip_has_no_extra_skeleton_voxels": graph.extra_voxel_count == 0,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
    }


def validate_hierarchy(result: NNE2HierarchyResult) -> dict[str, Any]:
    root = next(item for item in result.branches if item.branch_id == result.root_branch_id)
    represented = {item.branch_id for item in result.branches}
    classified = represented - set(result.unresolved_branch_ids)
    checks = {
        "has_root": result.root_node_id in result.directed_graph,
        "root_branch_is_diving_trunk_order_zero": root.branching_order == 0,
        "primary_parent_view_is_acyclic": nx.is_directed_acyclic_graph(result.primary_tree),
        "retained_root_component_is_weakly_connected": (
            result.directed_graph.number_of_nodes() > 0
            and nx.is_weakly_connected(result.directed_graph)
        ),
        "all_retained_branches_are_classified_or_unresolved": (
            classified | set(result.unresolved_branch_ids) == represented
        ),
        "has_at_least_one_measurement_anchor": bool(result.anchors),
        "has_matched_diving_trunk_anchor": any(
            item.branching_order == 0 and item.matched_branch_id is not None
            for item in result.anchors
        ),
    }
    warning_conditions = {
        "has_low_confidence_or_unresolved_branches": any(
            item.confidence in {"low", "unresolved"} for item in result.branches
        ),
        "has_order_conflicts": bool(result.order_conflict_branch_ids),
        "has_cross_links": bool(result.cross_link_branch_ids),
    }
    status = "FAIL" if not all(checks.values()) else (
        "WARN" if any(warning_conditions.values()) else "PASS"
    )
    return {
        "status": status,
        "checks": checks,
        "warnings": warning_conditions,
    }


def required_files_status(files: list[Path]) -> dict[str, Any]:
    missing = [str(item) for item in files if not item.is_file() or item.stat().st_size == 0]
    return {
        "status": "PASS" if not missing else "FAIL",
        "required_file_count": len(files),
        "missing_or_empty_files": missing,
    }
