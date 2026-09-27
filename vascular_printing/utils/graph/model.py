"""Data records for a geometry-aware hierarchical vascular graph."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import networkx as nx
import numpy as np


def _optional_float(value: float) -> float | None:
    return float(value) if np.isfinite(value) else None


@dataclass(slots=True)
class NodeRecord:
    node_id: int
    node_type: str
    graph_degree: int
    voxel_indices_xyz: np.ndarray
    representative_index_xyz: tuple[int, int, int]
    representative_lps_um: tuple[float, float, float]
    incident_branch_ids: list[int] = field(default_factory=list)
    cycle_ids: list[int] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "node_type": self.node_type,
            "graph_degree": self.graph_degree,
            "voxel_count": int(len(self.voxel_indices_xyz)),
            "representative_index_xyz": self.representative_index_xyz,
            "representative_lps_um": self.representative_lps_um,
            "incident_branch_ids": self.incident_branch_ids,
            "cycle_ids": self.cycle_ids,
        }


@dataclass(slots=True)
class BranchRecord:
    branch_id: int
    node_u: int
    node_v: int
    voxel_indices_xyz: np.ndarray
    points_raw_lps_um: np.ndarray
    arc_length_raw_um: np.ndarray
    coarse_radius_raw_um: np.ndarray
    points_smoothed_lps_um: np.ndarray
    arc_length_smoothed_um: np.ndarray
    coarse_radius_smoothed_um: np.ndarray
    local_direction_smoothed: np.ndarray
    curvature_smoothed_per_um: np.ndarray
    straight_distance_um: float
    tortuosity_raw: float
    tortuosity_smoothed: float
    smoothing_max_deviation_um: float
    cycle_ids: list[int] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        start_direction = (
            self.local_direction_smoothed[0].tolist()
            if len(self.local_direction_smoothed)
            else [0.0, 0.0, 0.0]
        )
        end_direction = (
            self.local_direction_smoothed[-1].tolist()
            if len(self.local_direction_smoothed)
            else [0.0, 0.0, 0.0]
        )
        return {
            "branch_id": self.branch_id,
            "node_u": self.node_u,
            "node_v": self.node_v,
            "point_count_raw": int(len(self.points_raw_lps_um)),
            "point_count_smoothed": int(len(self.points_smoothed_lps_um)),
            "length_raw_um": float(self.arc_length_raw_um[-1]),
            "length_smoothed_um": float(self.arc_length_smoothed_um[-1]),
            "straight_distance_um": float(self.straight_distance_um),
            "tortuosity_raw": _optional_float(self.tortuosity_raw),
            "tortuosity_smoothed": _optional_float(self.tortuosity_smoothed),
            "coarse_radius_mean_um": float(np.mean(self.coarse_radius_raw_um)),
            "coarse_radius_median_um": float(np.median(self.coarse_radius_raw_um)),
            "coarse_radius_proximal_um": float(self.coarse_radius_raw_um[0]),
            "coarse_radius_distal_um": float(self.coarse_radius_raw_um[-1]),
            "coarse_radius_cv": (
                float(np.std(self.coarse_radius_raw_um) / np.mean(self.coarse_radius_raw_um))
                if np.mean(self.coarse_radius_raw_um) > 0
                else None
            ),
            "coarse_taper_um_per_um": (
                float(
                    (self.coarse_radius_raw_um[-1] - self.coarse_radius_raw_um[0])
                    / self.arc_length_raw_um[-1]
                )
                if self.arc_length_raw_um[-1] > 0
                else None
            ),
            "curvature_mean_per_um": float(np.mean(self.curvature_smoothed_per_um)),
            "curvature_max_per_um": float(np.max(self.curvature_smoothed_per_um)),
            "initial_direction_lps": start_direction,
            "end_direction_lps": end_direction,
            "smoothing_max_deviation_um": float(self.smoothing_max_deviation_um),
            "cycle_ids": self.cycle_ids,
        }

    def dense_geometry(self) -> dict[str, Any]:
        return {
            "branch_id": self.branch_id,
            "voxel_indices_xyz": self.voxel_indices_xyz,
            "points_raw_lps_um": self.points_raw_lps_um,
            "arc_length_raw_um": self.arc_length_raw_um,
            "coarse_radius_raw_um": self.coarse_radius_raw_um,
            "points_smoothed_lps_um": self.points_smoothed_lps_um,
            "arc_length_smoothed_um": self.arc_length_smoothed_um,
            "coarse_radius_smoothed_um": self.coarse_radius_smoothed_um,
            "local_direction_smoothed": self.local_direction_smoothed,
            "curvature_smoothed_per_um": self.curvature_smoothed_per_um,
        }


@dataclass(slots=True)
class JunctionGeometry:
    node_id: int
    incident_directions: list[dict[str, Any]]
    pairwise_angles_deg: list[dict[str, Any]]

    def to_dict(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "incident_directions": self.incident_directions,
            "pairwise_angles_deg": self.pairwise_angles_deg,
        }


@dataclass(slots=True)
class CycleRecord:
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
class HierarchicalGraphResult:
    nodes: list[NodeRecord]
    branches: list[BranchRecord]
    junctions: list[JunctionGeometry]
    cycles: list[CycleRecord]
    junction_graph: nx.MultiGraph
    branch_as_node_graph: nx.MultiGraph
    source_skeleton: np.ndarray
    reconstructed_skeleton: np.ndarray
    skeleton_component_count: int
    graph_component_count: int
    skeleton_voxel_count: int
    represented_voxel_count: int
    missing_voxel_count: int
    extra_voxel_count: int
    duplicate_interior_voxel_count: int
    cycle_rank: int
    coordinate_system: str
    physical_unit: str
    origin_lps_um: tuple[float, float, float]
    spacing_um: tuple[float, float, float]
    radius_source: str

    def report(self) -> dict[str, Any]:
        branch_lengths = np.asarray(
            [item.arc_length_raw_um[-1] for item in self.branches], dtype=float
        )
        degrees = np.asarray([item.graph_degree for item in self.nodes], dtype=int)
        return {
            "representation_name": "Hierarchical Vascular Representation",
            "quality_level": "whole-network coarse navigation graph",
            "approved_for_cfd": False,
            "approved_as_final_geometry_training_truth": False,
            "coordinate_system": self.coordinate_system,
            "physical_unit": self.physical_unit,
            "origin_lps_um": self.origin_lps_um,
            "spacing_um": self.spacing_um,
            "radius_source": self.radius_source,
            "node_count": len(self.nodes),
            "terminal_node_count": int(sum(item.node_type == "terminal" for item in self.nodes)),
            "junction_node_count": int(sum(item.node_type == "junction" for item in self.nodes)),
            "complex_junction_node_count": int(
                sum(item.node_type == "complex_junction" for item in self.nodes)
            ),
            "connector_node_count": int(sum(item.node_type == "connector" for item in self.nodes)),
            "cycle_anchor_node_count": int(
                sum(item.node_type == "cycle_anchor" for item in self.nodes)
            ),
            "branch_count": len(self.branches),
            "cycle_basis_count": len(self.cycles),
            "cycle_rank": self.cycle_rank,
            "skeleton_component_count": self.skeleton_component_count,
            "graph_component_count": self.graph_component_count,
            "skeleton_voxel_count": self.skeleton_voxel_count,
            "represented_voxel_count": self.represented_voxel_count,
            "missing_voxel_count": self.missing_voxel_count,
            "extra_voxel_count": self.extra_voxel_count,
            "duplicate_interior_voxel_count": self.duplicate_interior_voxel_count,
            "branch_as_node_count": self.branch_as_node_graph.number_of_nodes(),
            "branch_relation_count": self.branch_as_node_graph.number_of_edges(),
            "node_degree_max": int(degrees.max()) if len(degrees) else 0,
            "branch_length_um": {
                "minimum": float(branch_lengths.min()) if len(branch_lengths) else 0.0,
                "median": float(np.median(branch_lengths)) if len(branch_lengths) else 0.0,
                "mean": float(branch_lengths.mean()) if len(branch_lengths) else 0.0,
                "maximum": float(branch_lengths.max()) if len(branch_lengths) else 0.0,
                "total": float(branch_lengths.sum()),
            },
        }
