"""Data records for NNE2 landmark matching and directed hierarchy results."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import networkx as nx


@dataclass(frozen=True, slots=True)
class AnchorMatch:
    record_id: int
    subject_id: str
    tree_id: int
    branching_order: int
    depth_um: float
    expected_stack_index: int
    matched_stack_index: int
    seed_xyz_um: tuple[float, float, float]
    registration_score: float
    registration_status: str
    anchor_pixel_method: str
    matched_branch_id: int | None
    branch_distance_um: float | None
    match_status: str
    candidate_branch_ids: tuple[int, ...] = ()
    candidate_distances_um: tuple[float, ...] = ()
    ambiguity_status: str = "not_evaluated"

    def to_dict(self) -> dict[str, Any]:
        return {
            "record_id": self.record_id,
            "subject_id": self.subject_id,
            "tree_id": self.tree_id,
            "branching_order": self.branching_order,
            "depth_um": self.depth_um,
            "expected_stack_index": self.expected_stack_index,
            "matched_stack_index": self.matched_stack_index,
            "seed_xyz_um": self.seed_xyz_um,
            "registration_score": self.registration_score,
            "registration_status": self.registration_status,
            "anchor_pixel_method": self.anchor_pixel_method,
            "matched_branch_id": self.matched_branch_id,
            "branch_distance_um": self.branch_distance_um,
            "match_status": self.match_status,
            "candidate_branch_ids": list(self.candidate_branch_ids),
            "candidate_distances_um": list(self.candidate_distances_um),
            "ambiguity_status": self.ambiguity_status,
        }


@dataclass(slots=True)
class DirectedBranch:
    branch_id: int
    node_u: int
    node_v: int
    upstream_node: int | None
    downstream_node: int | None
    direction_status: str
    confidence: str
    branching_order: int | None
    order_source: str
    anchor_record_ids: list[int] = field(default_factory=list)
    parent_branch_ids: list[int] = field(default_factory=list)
    child_branch_ids: list[int] = field(default_factory=list)
    is_primary_tree_edge: bool = False
    is_cross_link: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "branch_id": self.branch_id,
            "node_u": self.node_u,
            "node_v": self.node_v,
            "upstream_node": self.upstream_node,
            "downstream_node": self.downstream_node,
            "direction_status": self.direction_status,
            "confidence": self.confidence,
            "branching_order": self.branching_order,
            "order_source": self.order_source,
            "anchor_record_ids": self.anchor_record_ids,
            "parent_branch_ids": self.parent_branch_ids,
            "child_branch_ids": self.child_branch_ids,
            "is_primary_tree_edge": self.is_primary_tree_edge,
            "is_cross_link": self.is_cross_link,
        }


@dataclass(frozen=True, slots=True)
class ParentChildRelation:
    parent_branch_id: int
    child_branch_id: int
    shared_node_id: int
    relation_type: str
    confidence: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent_branch_id": self.parent_branch_id,
            "child_branch_id": self.child_branch_id,
            "shared_node_id": self.shared_node_id,
            "relation_type": self.relation_type,
            "confidence": self.confidence,
        }


@dataclass(slots=True)
class NNE2HierarchyResult:
    tree_key: str
    stack_name: str
    root_node_id: int
    root_branch_id: int
    branches: list[DirectedBranch]
    relations: list[ParentChildRelation]
    anchors: list[AnchorMatch]
    directed_graph: nx.MultiDiGraph
    primary_tree: nx.DiGraph
    unresolved_branch_ids: list[int]
    cross_link_branch_ids: list[int]
    order_conflict_branch_ids: list[int]
    excluded_branch_ids: list[int]

    def report(self) -> dict[str, Any]:
        confidence_counts: dict[str, int] = {}
        for item in self.branches:
            confidence_counts[item.confidence] = confidence_counts.get(item.confidence, 0) + 1
        return {
            "tree_key": self.tree_key,
            "stack_name": self.stack_name,
            "root_node_id": self.root_node_id,
            "root_branch_id": self.root_branch_id,
            "branch_count": len(self.branches),
            "parent_child_relation_count": len(self.relations),
            "matched_anchor_count": sum(item.matched_branch_id is not None for item in self.anchors),
            "anchor_count": len(self.anchors),
            "unresolved_branch_count": len(self.unresolved_branch_ids),
            "cross_link_branch_count": len(self.cross_link_branch_ids),
            "order_conflict_branch_count": len(self.order_conflict_branch_ids),
            "excluded_branch_count": len(self.excluded_branch_ids),
            "confidence_counts": confidence_counts,
            "primary_tree_is_directed_acyclic": nx.is_directed_acyclic_graph(self.primary_tree),
            "direction_interpretation": (
                "Anatomical direction inferred away from a Branching Order 0 diving-trunk "
                "anchor; it is not measured blood flow."
            ),
        }
