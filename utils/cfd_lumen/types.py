"""Typed records shared by the CFD lumen reconstruction stages."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


class GeometryValidationError(RuntimeError):
    """Raised when source geometry or a quantitative CFD acceptance check fails."""


@dataclass(slots=True)
class BranchGeometry:
    branch_id: int
    local_node_ids: tuple[int, ...]
    source_global_nodes: tuple[int, ...]
    source_global_edges: tuple[int, ...]
    raw_points_um: np.ndarray
    raw_radius_um: np.ndarray
    points_um: np.ndarray = field(default_factory=lambda: np.empty((0, 3), dtype=float))
    radius_um: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    arc_length_um: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    construction_point_type: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.uint8)
    )
    construction_source_cut_port_id: tuple[str, ...] = ()
    construction_source_global_edge_id: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )
    construction_distance_from_core_boundary_um: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    construction_original_core_cut_position_um: np.ndarray = field(
        default_factory=lambda: np.empty((0, 3), dtype=float)
    )

    @property
    def length_um(self) -> float:
        points = self.points_um if len(self.points_um) else self.raw_points_um
        return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


@dataclass(frozen=True, slots=True)
class PortGeometry:
    port_id: int
    roi_id: str
    cut_port_id: str
    local_node_id: int
    global_edge_id: int
    original_position_um: np.ndarray
    cap_center_um: np.ndarray
    radius_um: float
    outward_tangent: np.ndarray
    extension_length_um: float
    overlap_length_um: float
    cylinder_start_um: np.ndarray
    cylinder_end_um: np.ndarray
    boundary_face: str
    boundary_role: str = "CUT_PORT_UNASSIGNED"
    source_core_cut_port_id: str | None = None
    original_core_cut_position_um: tuple[float, float, float] | None = None
    source_core_to_cfd_cut_length_um: float = 0.0

    def metadata(self, *, patch_area_um2: float | None = None) -> dict[str, Any]:
        return {
            "port_id": self.port_id,
            "roi_id": self.roi_id,
            "cut_port_id": self.cut_port_id,
            "global_edge_id": self.global_edge_id,
            "original_x_um": float(self.original_position_um[0]),
            "original_y_um": float(self.original_position_um[1]),
            "original_z_um": float(self.original_position_um[2]),
            "cap_x_um": float(self.cap_center_um[0]),
            "cap_y_um": float(self.cap_center_um[1]),
            "cap_z_um": float(self.cap_center_um[2]),
            "radius_um": self.radius_um,
            "diameter_um": 2.0 * self.radius_um,
            "outward_tx": float(self.outward_tangent[0]),
            "outward_ty": float(self.outward_tangent[1]),
            "outward_tz": float(self.outward_tangent[2]),
            "extension_length_um": self.extension_length_um,
            "patch_area_um2": patch_area_um2,
            "boundary_role": self.boundary_role,
            "source_core_cut_port_id": self.source_core_cut_port_id,
              "original_core_cut_position_um": self.original_core_cut_position_um,
              "source_core_to_cfd_cut_length_um": self.source_core_to_cfd_cut_length_um,
            "P_1D": None,
            "Q_1D": None,
        }


@dataclass(frozen=True, slots=True)
class CollisionEvent:
    branch_id_a: int
    branch_id_b: int
    segment_index_a: int
    segment_index_b: int
    distance_um: float
    radius_a_um: float
    radius_b_um: float
    clearance_um: float
    classification: str
    closest_a_um: tuple[float, float, float]
    closest_b_um: tuple[float, float, float]

    def report(self) -> dict[str, Any]:
        return {
            "branch_id_a": self.branch_id_a,
            "branch_id_b": self.branch_id_b,
            "segment_index_a": self.segment_index_a,
            "segment_index_b": self.segment_index_b,
            "distance_um": self.distance_um,
            "radius_a_um": self.radius_a_um,
            "radius_b_um": self.radius_b_um,
            "clearance_um": self.clearance_um,
            "classification": self.classification,
            "closest_a_x_um": self.closest_a_um[0],
            "closest_a_y_um": self.closest_a_um[1],
            "closest_a_z_um": self.closest_a_um[2],
            "closest_b_x_um": self.closest_b_um[0],
            "closest_b_y_um": self.closest_b_um[1],
            "closest_b_z_um": self.closest_b_um[2],
        }


@dataclass(frozen=True, slots=True)
class JunctionCollar:
    junction_node_id: int
    branch_id: int
    endpoint_index: int
    collar_distance_um: float
    overlap_length_um: float
    implicit_extent_um: float
    explicit_cap_distance_um: float
    collar_position_um: np.ndarray
    collar_radius_um: float
    tangent_away: np.ndarray


@dataclass(slots=True)
class LocalJunctionPatch:
    junction_node_id: int
    raw_mesh: Any
    clean_mesh: Any
    collars: list[JunctionCollar]
    centerline_points_um: np.ndarray
    centerline_radius_um: np.ndarray
    metadata: dict[str, Any]


@dataclass(slots=True)
class HybridBuildDetails:
    patches: dict[int, LocalJunctionPatch]
    merged_junction_meshes: dict[int, Any]
    trimmed_branches: list[BranchGeometry]
    merge_steps: list[dict[str, Any]]
    explicit_cap_planes: list[dict[str, Any]]
    runtime_s: dict[str, float]
    transition_backend: str = "manifold_boolean"
    face_region: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.uint8))
    transition_rows: list[dict[str, Any]] = field(default_factory=list)
    constructed_branches: list[BranchGeometry] = field(default_factory=list)
    port_extension_rows: list[dict[str, Any]] = field(default_factory=list)
    transition_fallback_reason: str | None = None


@dataclass(slots=True)
class ContextDomainResult:
    core_roi: Any
    cfd_roi: Any
    source_manifest: Path | None
    port_mappings: list[dict[str, Any]]
    added_global_edge_ids: tuple[int, ...]
    added_global_branch_ids: tuple[int, ...]
    added_centerline_length_um: float
    report: dict[str, Any]


@dataclass(slots=True)
class SurfaceBuildResult:
    mesh: Any
    backend_requested: str
    backend_used: str
    fallback_reason: str | None
    tube_primitive_count: int
    junction_primitive_count: int
    extension_primitive_count: int
    boolean_runtime_s: float
    implicit_grid: dict[str, Any] | None = None
    hybrid_details: HybridBuildDetails | None = None
    surface_continuity_version: str = "v4"


@dataclass(slots=True)
class LumenPrimitives:
    """Boolean-input solids retained for v2 root-cause diagnostics."""

    branch_tubes: dict[int, Any]
    junction_solids: dict[int, Any]
    port_extensions: dict[int, Any]

    @property
    def all_meshes(self) -> list[Any]:
        return [
            *self.branch_tubes.values(),
            *self.junction_solids.values(),
            *self.port_extensions.values(),
        ]


@dataclass(slots=True)
class PatchResult:
    patch_id: np.ndarray
    patch_type: np.ndarray
    port_id: np.ndarray
    port_rows: list[dict[str, Any]]
    detected_port_count: int
    all_ports_pass: bool


@dataclass(frozen=True, slots=True)
class RadiusFidelitySample:
    branch_id: int
    sample_index: int
    arc_length_um: float
    center_um: tuple[float, float, float]
    tangent: tuple[float, float, float]
    source_radius_um: float
    reconstructed_radius_um: float
    relative_error: float
    section_xy_um: tuple[tuple[float, float], ...] = ()

    def report(self) -> dict[str, Any]:
        return {
            "branch_id": self.branch_id,
            "sample_index": self.sample_index,
            "arc_length_um": self.arc_length_um,
            "center_x_um": self.center_um[0],
            "center_y_um": self.center_um[1],
            "center_z_um": self.center_um[2],
            "tangent_x": self.tangent[0],
            "tangent_y": self.tangent[1],
            "tangent_z": self.tangent[2],
            "source_radius_um": self.source_radius_um,
            "reconstructed_radius_um": self.reconstructed_radius_um,
            "relative_error": self.relative_error,
            "absolute_relative_error": abs(self.relative_error),
        }


@dataclass(slots=True)
class ROIProcessResult:
    roi_id: str
    status: str
    output_dir: Path
    summary: dict[str, Any]
    failure_reason: str | None = None


@dataclass(frozen=True, slots=True)
class CFDRunLayout:
    run_root: Path
    logs: Path
    config: Path
    manifests: Path
    figures: Path
    report: Path
    rois: Path
