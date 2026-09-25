"""In-memory records for a directed, geometry-aware Schmid vascular graph."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import networkx as nx
import numpy as np


VESSEL_TYPE_NAMES = {
    0: "pial_artery",
    1: "pial_vein",
    2: "descending_arteriole",
    3: "ascending_venule",
    4: "capillary",
    5: "unknown",
}


def optional_float(value: float | None) -> float | None:
    if value is None:
        return None
    return float(value) if np.isfinite(value) else None


@dataclass(slots=True)
class SchmidInputData:
    coordinates_um: np.ndarray
    pressure_mmhg: np.ndarray
    pressure_boundary_mmhg: np.ndarray
    edge_tuples: np.ndarray
    mean_diameter_um: np.ndarray
    flow_um3_per_ms: np.ndarray
    hematocrit_boundary: np.ndarray
    vessel_type_code: np.ndarray
    source_length_um: np.ndarray
    hematocrit: np.ndarray
    red_blood_cell_count: np.ndarray
    diameter_profiles_um: list[np.ndarray]
    point_sequences_um: list[np.ndarray]
    source_files: dict[str, dict[str, Any]]

    @property
    def vertex_count(self) -> int:
        return int(len(self.coordinates_um))

    @property
    def edge_count(self) -> int:
        return int(len(self.edge_tuples))


@dataclass(slots=True)
class CleanEdge:
    edge_id: int
    node_u: int
    node_v: int
    points_u_to_v_um: np.ndarray
    diameter_u_to_v_um: np.ndarray
    mean_diameter_um: float
    source_length_um: float
    geometric_length_um: float
    flow_um3_per_ms: float
    vessel_type_code: int
    hematocrit: float | None
    hematocrit_boundary: float | None
    red_blood_cell_count: float | None
    geometry_reversed: bool
    geometry_status: str
    endpoint_error_um: float
    length_relative_error: float | None
    upstream_node: int | None = None
    downstream_node: int | None = None
    pressure_drop_mmhg: float | None = None
    direction_status: str = "unassigned"

    def topology_summary(self) -> dict[str, Any]:
        return {
            "edge_id": self.edge_id,
            "node_u": self.node_u,
            "node_v": self.node_v,
            "upstream_node": self.upstream_node,
            "downstream_node": self.downstream_node,
            "direction_status": self.direction_status,
            "pressure_drop_mmhg": self.pressure_drop_mmhg,
            "flow_um3_per_ms": self.flow_um3_per_ms,
            "source_length_um": self.source_length_um,
            "geometric_length_um": self.geometric_length_um,
            "mean_diameter_um": self.mean_diameter_um,
            "vessel_type_code": self.vessel_type_code,
            "vessel_type": VESSEL_TYPE_NAMES.get(self.vessel_type_code, "unrecognized"),
            "geometry_status": self.geometry_status,
            "geometry_reversed": self.geometry_reversed,
            "endpoint_error_um": self.endpoint_error_um,
            "length_relative_error": self.length_relative_error,
        }


@dataclass(slots=True)
class CleanupDecision:
    item_type: str
    item_id: int
    action: str
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "item_type": self.item_type,
            "item_id": self.item_id,
            "action": self.action,
            "reason": self.reason,
        }


@dataclass(slots=True)
class SchmidCleanupResult:
    source: SchmidInputData
    valid_node_ids: np.ndarray
    edges: list[CleanEdge]
    decisions: list[CleanupDecision]
    original_component_count: int
    retained_component_count: int
    removed_node_ids: list[int]
    removed_edge_ids: list[int]

    def report(self) -> dict[str, Any]:
        action_counts: dict[str, int] = {}
        geometry_counts: dict[str, int] = {}
        for item in self.decisions:
            action_counts[item.action] = action_counts.get(item.action, 0) + 1
        for edge in self.edges:
            geometry_counts[edge.geometry_status] = geometry_counts.get(edge.geometry_status, 0) + 1
        return {
            "source_vertex_count": self.source.vertex_count,
            "source_edge_count": self.source.edge_count,
            "retained_vertex_count": int(len(self.valid_node_ids)),
            "retained_edge_count": len(self.edges),
            "removed_vertex_count": len(self.removed_node_ids),
            "removed_edge_count": len(self.removed_edge_ids),
            "original_component_count": self.original_component_count,
            "retained_component_count": self.retained_component_count,
            "decision_action_counts": action_counts,
            "geometry_status_counts": geometry_counts,
        }


@dataclass(slots=True)
class DirectedNode:
    node_id: int
    coordinates_um: tuple[float, float, float]
    pressure_mmhg: float
    pressure_boundary_mmhg: float | None
    node_role: str = "unclassified"
    undirected_degree: int = 0
    incoming_branch_ids: list[int] = field(default_factory=list)
    outgoing_branch_ids: list[int] = field(default_factory=list)
    unresolved_branch_ids: list[int] = field(default_factory=list)
    parent_node_ids: list[int] = field(default_factory=list)
    child_node_ids: list[int] = field(default_factory=list)
    cycle_ids: list[int] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "node_role": self.node_role,
            "coordinates_source_xyz_um": self.coordinates_um,
            "pressure_mmhg": self.pressure_mmhg,
            "pressure_boundary_mmhg": self.pressure_boundary_mmhg,
            "undirected_degree": self.undirected_degree,
            "directed_indegree": len(self.incoming_branch_ids),
            "directed_outdegree": len(self.outgoing_branch_ids),
            "incoming_branch_ids": self.incoming_branch_ids,
            "outgoing_branch_ids": self.outgoing_branch_ids,
            "unresolved_branch_ids": self.unresolved_branch_ids,
            "parent_node_ids": self.parent_node_ids,
            "child_node_ids": self.child_node_ids,
            "cycle_ids": self.cycle_ids,
        }


@dataclass(slots=True)
class DirectedBranch:
    branch_id: int
    raw_edge_ids: list[int]
    node_a: int
    node_b: int
    upstream_node: int | None
    downstream_node: int | None
    direction_status: str
    points_raw_um: np.ndarray
    radius_raw_um: np.ndarray
    points_smoothed_um: np.ndarray
    radius_smoothed_um: np.ndarray
    arc_length_raw_um: np.ndarray
    arc_length_smoothed_um: np.ndarray
    local_direction_smoothed: np.ndarray
    curvature_smoothed_per_um: np.ndarray
    source_length_um: float
    geometric_length_um: float
    straight_distance_um: float
    tortuosity_raw: float | None
    tortuosity_smoothed: float | None
    smoothing_max_deviation_um: float
    flow_um3_per_ms: float
    flow_min_um3_per_ms: float
    flow_max_um3_per_ms: float
    pressure_drop_mmhg: float | None
    vessel_type_codes: list[int]
    geometry_status: str
    parent_branch_ids: list[int] = field(default_factory=list)
    child_branch_ids: list[int] = field(default_factory=list)
    cycle_ids: list[int] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        radii = self.radius_raw_um[np.isfinite(self.radius_raw_um)]
        radius_mean = float(np.mean(radii)) if len(radii) else None
        radius_median = float(np.median(radii)) if len(radii) else None
        return {
            "branch_id": self.branch_id,
            "raw_edge_ids": self.raw_edge_ids,
            "raw_edge_count": len(self.raw_edge_ids),
            "node_a": self.node_a,
            "node_b": self.node_b,
            "upstream_node": self.upstream_node,
            "downstream_node": self.downstream_node,
            "direction_status": self.direction_status,
            "parent_branch_ids": self.parent_branch_ids,
            "child_branch_ids": self.child_branch_ids,
            "point_count_raw": len(self.points_raw_um),
            "point_count_smoothed": len(self.points_smoothed_um),
            "source_length_um": self.source_length_um,
            "geometric_length_um": self.geometric_length_um,
            "smoothed_length_um": float(self.arc_length_smoothed_um[-1]) if len(self.arc_length_smoothed_um) else 0.0,
            "straight_distance_um": self.straight_distance_um,
            "tortuosity_raw": self.tortuosity_raw,
            "tortuosity_smoothed": self.tortuosity_smoothed,
            "mean_radius_um": radius_mean,
            "median_radius_um": radius_median,
            "proximal_radius_um": float(radii[0]) if len(radii) and self.direction_status == "known" else None,
            "distal_radius_um": float(radii[-1]) if len(radii) and self.direction_status == "known" else None,
            "flow_um3_per_ms": self.flow_um3_per_ms,
            "flow_min_um3_per_ms": self.flow_min_um3_per_ms,
            "flow_max_um3_per_ms": self.flow_max_um3_per_ms,
            "pressure_drop_mmhg": self.pressure_drop_mmhg,
            "vessel_type_codes": self.vessel_type_codes,
            "vessel_types": [VESSEL_TYPE_NAMES.get(code, "unrecognized") for code in self.vessel_type_codes],
            "curvature_mean_per_um": float(np.mean(self.curvature_smoothed_per_um)) if len(self.curvature_smoothed_per_um) else 0.0,
            "curvature_max_per_um": float(np.max(self.curvature_smoothed_per_um)) if len(self.curvature_smoothed_per_um) else 0.0,
            "smoothing_max_deviation_um": self.smoothing_max_deviation_um,
            "geometry_status": self.geometry_status,
            "cycle_ids": self.cycle_ids,
        }


@dataclass(slots=True)
class DirectedCycle:
    cycle_id: int
    cycle_type: str
    node_ids: list[int]
    branch_ids: list[int]

    def to_dict(self) -> dict[str, Any]:
        return {
            "cycle_id": self.cycle_id,
            "cycle_type": self.cycle_type,
            "node_ids": self.node_ids,
            "branch_ids": self.branch_ids,
        }


@dataclass(slots=True)
class DirectedGraphResult:
    nodes: list[DirectedNode]
    branches: list[DirectedBranch]
    cycles: list[DirectedCycle]
    all_connectivity_graph: nx.MultiGraph
    directed_junction_graph: nx.MultiDiGraph
    directed_branch_graph: nx.MultiDiGraph
    raw_directed_graph: nx.MultiDiGraph
    raw_unresolved_edge_ids: list[int]
    raw_edge_to_branch: dict[int, int]
    weak_component_count: int
    directed_is_acyclic: bool
    cleanup: SchmidCleanupResult
    flow_conservation: list[dict[str, Any]]

    def report(self) -> dict[str, Any]:
        role_counts: dict[str, int] = {}
        for node in self.nodes:
            role_counts[node.node_role] = role_counts.get(node.node_role, 0) + 1
        known = sum(branch.direction_status == "known" for branch in self.branches)
        unresolved = len(self.branches) - known
        branch_lengths = np.asarray([item.source_length_um for item in self.branches], dtype=float)
        return {
            "representation_name": "Directed Schmid Hierarchical Vascular Graph",
            "source_dataset": "Schmid NW1_results",
            "coordinate_system": "source_XYZ_not_anatomically_labeled",
            "physical_unit": "micrometer",
            "pressure_unit": "mmHg",
            "flow_unit": "micrometer^3/millisecond",
            "directed": True,
            "direction_method": "higher endpoint pressure to lower endpoint pressure",
            "is_tree": False,
            "allows_multiple_parents_and_children": True,
            "node_count": len(self.nodes),
            "branch_count": len(self.branches),
            "known_direction_branch_count": known,
            "unresolved_direction_branch_count": unresolved,
            "raw_unresolved_edge_count": len(self.raw_unresolved_edge_ids),
            "weak_component_count": self.weak_component_count,
            "undirected_cycle_rank": (
                self.all_connectivity_graph.number_of_edges()
                - self.all_connectivity_graph.number_of_nodes()
                + self.weak_component_count
            ),
            "stored_cycle_basis_count": len(self.cycles),
            "known_direction_graph_is_acyclic": self.directed_is_acyclic,
            "node_role_counts": role_counts,
            "branch_length_um": {
                "minimum": float(branch_lengths.min()) if len(branch_lengths) else 0.0,
                "median": float(np.median(branch_lengths)) if len(branch_lengths) else 0.0,
                "mean": float(branch_lengths.mean()) if len(branch_lengths) else 0.0,
                "maximum": float(branch_lengths.max()) if len(branch_lengths) else 0.0,
                "total": float(branch_lengths.sum()),
            },
        }
