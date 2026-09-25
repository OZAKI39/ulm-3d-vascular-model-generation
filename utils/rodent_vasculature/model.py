"""Data models for a directed SWC-derived vascular representation."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import networkx as nx
import numpy as np

from .swc_io import SWCData


@dataclass(slots=True)
class DirectedBranch:
    branch_id: int
    upstream_node_id: int
    downstream_node_id: int
    source_node_ids: list[int]
    raw_points_voxel_xyz: np.ndarray
    raw_points_um: np.ndarray
    raw_radius_um: np.ndarray
    derived_points_um: np.ndarray
    derived_radius_um: np.ndarray
    direction_vectors_xyz: np.ndarray
    raw_length_um: float
    derived_length_um: float
    tortuosity: float
    mean_radius_um: float
    component_id: int = -1
    root_node_id: int = -1
    depth: int = 0
    parent_branch_ids: list[int] = field(default_factory=list)
    daughter_branch_ids: list[int] = field(default_factory=list)
    strahler_order: int = 1
    horsfield_order: int = 1
    downstream_terminal_count: int = 1
    downstream_branch_count: int = 1
    proximal_boundary_contact: bool = False
    distal_boundary_contact: bool = False

    def summary(self) -> dict[str, Any]:
        return {
            "branch_id": self.branch_id,
            "upstream_node_id": self.upstream_node_id,
            "downstream_node_id": self.downstream_node_id,
            "source_node_ids": ";".join(map(str, self.source_node_ids)),
            "point_count": len(self.source_node_ids),
            "component_id": self.component_id,
            "root_node_id": self.root_node_id,
            "depth": self.depth,
            "parent_branch_ids": ";".join(map(str, self.parent_branch_ids)),
            "daughter_branch_ids": ";".join(map(str, self.daughter_branch_ids)),
            "strahler_order": self.strahler_order,
            "horsfield_order": self.horsfield_order,
            "downstream_terminal_count": self.downstream_terminal_count,
            "downstream_branch_count": self.downstream_branch_count,
            "raw_length_um": self.raw_length_um,
            "derived_length_um": self.derived_length_um,
            "tortuosity": self.tortuosity,
            "mean_radius_um": self.mean_radius_um,
            "proximal_boundary_contact": self.proximal_boundary_contact,
            "distal_boundary_contact": self.distal_boundary_contact,
            "flow_direction_rule": "upstream_node_id -> downstream_node_id",
        }


@dataclass(slots=True)
class DirectedVascularGraph:
    sample_id: str
    swc: SWCData
    source_graph: nx.DiGraph
    junction_graph: nx.MultiDiGraph
    branch_graph: nx.DiGraph
    branches: list[DirectedBranch]
    source_node_to_branch_ids: dict[int, list[int]]
    warnings: list[str] = field(default_factory=list)

    def report(self) -> dict[str, Any]:
        roles: dict[str, int] = {}
        for _, attributes in self.junction_graph.nodes(data=True):
            role = str(attributes["role"])
            roles[role] = roles.get(role, 0) + 1
        return {
            "sample_id": self.sample_id,
            "source_node_count": self.source_graph.number_of_nodes(),
            "source_edge_count": self.source_graph.number_of_edges(),
            "critical_node_count": self.junction_graph.number_of_nodes(),
            "branch_count": len(self.branches),
            "branch_relation_count": self.branch_graph.number_of_edges(),
            "component_count": self.swc.component_count,
            "root_count": len(self.swc.root_ids),
            "inferred_inlet_count": roles.get("inferred_inlet", 0),
            "inferred_outlet_count": roles.get("inferred_outlet", 0),
            "node_role_counts": roles,
            "flow_direction_rule": "SWC parent_id node -> current node",
            "flow_direction_is_measured": False,
            "warnings": self.warnings,
        }
