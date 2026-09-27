"""Typed data objects shared by the real vascular ROI sampling pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True, slots=True)
class GlobalEdge:
    """One stable, traceable SWC parent-to-current edge."""

    edge_id: int
    upstream_node_id: int
    downstream_node_id: int
    upstream_position_um: np.ndarray
    downstream_position_um: np.ndarray
    upstream_radius_um: float
    downstream_radius_um: float


@dataclass(slots=True)
class GlobalVascularModel:
    """Read-only source arrays and graph metadata used by spatial extraction."""

    source_model_id: str
    source_mouse_id: str
    node_ids: np.ndarray
    node_positions_um: np.ndarray
    node_radius_um: np.ndarray
    parent_ids: np.ndarray
    edges: tuple[GlobalEdge, ...]
    model_bounds_xyz_um: tuple[float, float, float, float, float, float]
    node_index_by_id: dict[int, int]
    incident_edge_ids_by_node: dict[int, tuple[int, ...]]
    global_degree_by_node: dict[int, int]

    @property
    def node_count(self) -> int:
        return int(len(self.node_ids))

    @property
    def edge_count(self) -> int:
        return len(self.edges)


@dataclass(frozen=True, slots=True)
class CutPort:
    """An exact ROI-boundary intersection traced to one global SWC edge."""

    cut_port_id: str
    local_node_id: int
    global_edge_id: int
    intersection_position_um: tuple[float, float, float]
    radius_at_cut_um: float
    boundary_face: str
    boundary_role: str = "CORE_ROI_BOUNDARY"

    def report(self) -> dict[str, Any]:
        return {
            "cut_port_id": self.cut_port_id,
            "local_node_id": self.local_node_id,
            "global_edge_id": self.global_edge_id,
            "intersection_x_um": self.intersection_position_um[0],
            "intersection_y_um": self.intersection_position_um[1],
            "intersection_z_um": self.intersection_position_um[2],
            "radius_at_cut_um": self.radius_at_cut_um,
            "boundary_face": self.boundary_face,
            "boundary_role": self.boundary_role,
        }


@dataclass(slots=True)
class ROIRecord:
    """A connected, source-traceable, spatially clipped vascular subgraph."""

    roi_id: str
    source_model_id: str
    source_mouse_id: str
    anchor_id: int
    anchor_position_um: tuple[float, float, float]
    bbox_min_um: tuple[float, float, float]
    bbox_max_um: tuple[float, float, float]
    bbox_center_um: tuple[float, float, float]
    bbox_size_um: tuple[float, float, float]
    global_node_ids: tuple[int, ...]
    global_edge_ids: tuple[int, ...]
    local_node_ids: np.ndarray
    local_node_global_ids: np.ndarray
    local_node_positions_um: np.ndarray
    local_node_radius_um: np.ndarray
    local_edges: np.ndarray
    local_edge_ids: np.ndarray
    local_edge_global_ids: np.ndarray
    local_edge_points_um: np.ndarray
    local_edge_radius_um: np.ndarray
    true_terminal_local_ids: tuple[int, ...]
    true_terminal_global_ids: tuple[int, ...]
    cut_ports: tuple[CutPort, ...]
    raw_component_count: int
    raw_total_vessel_length_um: float
    retained_component_length_um: float
    radius_features: dict[str, float] = field(default_factory=dict)
    structural_features: dict[str, float] = field(default_factory=dict)
    feature_vector: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    cluster_id: int = -1
    distance_to_cluster_center: float = float("nan")
    is_representative: bool = False
    selection_rank: int = -1
    boundary_conditions: Any | None = None
    flow_direction: Any | None = None
    pressure: Any | None = None
    flow_rate: Any | None = None
    velocity_field: Any | None = None
    mb_trajectories: Any | None = None

    @property
    def node_count(self) -> int:
        return int(len(self.local_node_ids))

    @property
    def edge_count(self) -> int:
        return int(len(self.local_edges))

    @property
    def branch_count(self) -> int:
        return int(round(self.structural_features.get("branch_count", 0.0)))

    @property
    def bifurcation_count(self) -> int:
        return int(round(self.structural_features.get("bifurcation_count", 0.0)))

    @property
    def true_terminal_count(self) -> int:
        return len(self.true_terminal_local_ids)

    @property
    def cut_port_count(self) -> int:
        return len(self.cut_ports)

    def manifest_row(self) -> dict[str, Any]:
        return {
            "roi_id": self.roi_id,
            "source_model_id": self.source_model_id,
            "source_mouse_id": self.source_mouse_id,
            "anchor_id": self.anchor_id,
            "anchor_x_um": self.anchor_position_um[0],
            "anchor_y_um": self.anchor_position_um[1],
            "anchor_z_um": self.anchor_position_um[2],
            "bbox_min_x_um": self.bbox_min_um[0],
            "bbox_min_y_um": self.bbox_min_um[1],
            "bbox_min_z_um": self.bbox_min_um[2],
            "bbox_max_x_um": self.bbox_max_um[0],
            "bbox_max_y_um": self.bbox_max_um[1],
            "bbox_max_z_um": self.bbox_max_um[2],
            "bbox_center_x_um": self.bbox_center_um[0],
            "bbox_center_y_um": self.bbox_center_um[1],
            "bbox_center_z_um": self.bbox_center_um[2],
            "bbox_size_x_um": self.bbox_size_um[0],
            "bbox_size_y_um": self.bbox_size_um[1],
            "bbox_size_z_um": self.bbox_size_um[2],
            "global_node_ids": ";".join(map(str, self.global_node_ids)),
            "global_edge_ids": ";".join(map(str, self.global_edge_ids)),
            "local_node_ids": ";".join(map(str, self.local_node_ids.tolist())),
            "local_edge_ids": ";".join(map(str, self.local_edge_ids.tolist())),
            "node_count": self.node_count,
            "edge_count": self.edge_count,
            "branch_count": self.branch_count,
            "bifurcation_count": self.bifurcation_count,
            "true_terminal_ids": ";".join(map(str, self.true_terminal_global_ids)),
            "cut_port_ids": ";".join(port.cut_port_id for port in self.cut_ports),
            "true_terminal_count": self.true_terminal_count,
            "cut_port_count": self.cut_port_count,
            "raw_component_count": self.raw_component_count,
            "raw_total_vessel_length_um": self.raw_total_vessel_length_um,
            "retained_component_length_um": self.retained_component_length_um,
            "cluster_id": self.cluster_id,
            "distance_to_center": self.distance_to_cluster_center,
            "selected": self.is_representative,
            "selection_rank": self.selection_rank,
        }

    def feature_row(self) -> dict[str, Any]:
        row = self.manifest_row()
        row.update(self.radius_features)
        row.update(self.structural_features)
        row.update(
            {
                "anchor_x": self.anchor_position_um[0],
                "anchor_y": self.anchor_position_um[1],
                "anchor_z": self.anchor_position_um[2],
                "total_vessel_length": self.structural_features.get(
                    "total_vessel_length_um", float("nan")
                ),
                "vessel_length_density": self.structural_features.get(
                    "vessel_length_density_um_per_um3", float("nan")
                ),
            }
        )
        return row


@dataclass(frozen=True, slots=True)
class ScalerState:
    feature_names: tuple[str, ...]
    median: np.ndarray
    iqr: np.ndarray

    def report(self) -> dict[str, Any]:
        return {
            "method": "robust",
            "feature_order": list(self.feature_names),
            "median": self.median.tolist(),
            "iqr": self.iqr.tolist(),
            "zero_iqr_replacement": 1.0,
        }


@dataclass(frozen=True, slots=True)
class ClusteringResult:
    method: str
    n_clusters: int
    feature_names: tuple[str, ...]
    scaled_features: np.ndarray
    assignments: np.ndarray
    centers: np.ndarray
    distances_to_center: np.ndarray
    inertia: float
    silhouette_score: float | None
    cluster_sizes: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class SelectionResult:
    selection_mode: str
    scientific_label: str
    selected_indices: tuple[int, ...]
    overlap_rejection_count: int
    requested_count: int


@dataclass(frozen=True, slots=True)
class SamplingExperiment:
    feature_mode: str
    feature_names: tuple[str, ...]
    scaler: ScalerState
    clustering: ClusteringResult
    selection: SelectionResult
    validation: dict[str, Any]


@dataclass(slots=True)
class SamplingRunResult:
    run_root: Path
    status: str
    candidates: list[ROIRecord]
    primary_experiment: SamplingExperiment
    comparison_experiments: dict[str, SamplingExperiment]
    summary_path: Path
    log_path: Path
    figure_paths: list[Path]
