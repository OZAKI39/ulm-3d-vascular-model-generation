"""Shared surface and boundary result records."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import trimesh


@dataclass(slots=True)
class TaggedSurface:
    vertices: np.ndarray
    faces: np.ndarray
    boundary_type: np.ndarray
    boundary_index: np.ndarray
    boundary_origin: np.ndarray
    face_kind: np.ndarray
    extension_index: np.ndarray
    extension_band: np.ndarray
    source_vertex_index: np.ndarray

    @classmethod
    def from_mesh(cls, mesh: trimesh.Trimesh) -> "TaggedSurface":
        count = len(mesh.faces)
        return cls(
            vertices=np.asarray(mesh.vertices, dtype=float).copy(),
            faces=np.asarray(mesh.faces, dtype=np.int64).copy(),
            boundary_type=np.zeros(count, dtype=np.uint8),
            boundary_index=np.full(count, -1, dtype=np.int32),
            boundary_origin=np.zeros(count, dtype=np.uint8),
            face_kind=np.zeros(count, dtype=np.uint8),
            extension_index=np.full(count, -1, dtype=np.int32),
            extension_band=np.full(count, -1, dtype=np.int32),
            source_vertex_index=np.arange(len(mesh.vertices), dtype=np.int64),
        )

    def mesh(self, *, compact: bool = False) -> trimesh.Trimesh:
        mesh = trimesh.Trimesh(self.vertices.copy(), self.faces.copy(), process=False)
        if compact:
            mesh.remove_unreferenced_vertices()
        return mesh

    def compact(self) -> None:
        used = np.unique(self.faces.ravel())
        mapping = np.full(len(self.vertices), -1, dtype=np.int64)
        mapping[used] = np.arange(len(used), dtype=np.int64)
        self.vertices = self.vertices[used]
        self.faces = mapping[self.faces]
        self.source_vertex_index = self.source_vertex_index[used]


@dataclass(frozen=True, slots=True)
class CutLoop:
    boundary_index: int
    vertex_ids: np.ndarray
    center_um: np.ndarray
    outward_normal: np.ndarray
    area_um2: float
    equivalent_radius_um: float
    plane_residual_um: float


@dataclass(frozen=True, slots=True)
class BoundarySurfaceResult:
    boundary_index: int
    port_id: str
    boundary_origin: str
    role: str
    source_radius_um: float
    extension_length_um: float
    actual_cap_area_um2: float
    equivalent_radius_um: float
    cap_planarity_error_um: float
    minimum_cap_normal_dot: float
    extension_length_error_um: float
    extension_axis_dot: float
    intermediate_ring_centerline_max_deviation_um: float
    intermediate_ring_centerline_p95_deviation_um: float
    intermediate_ring_centerline_mean_deviation_um: float
    intermediate_ring_centerline_worst_ring_index: int
    intermediate_ring_centerline_worst_ring_axial_station_um: float
    intermediate_ring_axial_station_max_error_um: float
    intermediate_ring_area_relative_error_max: float
    intermediate_ring_all_areas_finite_positive: bool
    intermediate_ring_all_polygons_simple_valid: bool
    intermediate_ring_all_orientations_consistent: bool
    local_original_median_edge_length_um: float
    target_edge_length_um: float
    ring_count: int
    transition_length_um: float
    proximal_ring_max_motion_um: float
    distal_ring_max_motion_um: float
    cut_loop_vertex_count: int
    cap_face_indices: np.ndarray
    side_face_indices: np.ndarray

    def report(self) -> dict[str, object]:
        return {
            "boundary_index": self.boundary_index,
            "port_id": self.port_id,
            "boundary_origin": self.boundary_origin,
            "role": self.role,
            "source_radius_um": self.source_radius_um,
            "extension_length_um": self.extension_length_um,
            "actual_cap_area_um2": self.actual_cap_area_um2,
            "equivalent_radius_um": self.equivalent_radius_um,
            "cap_planarity_error_um": self.cap_planarity_error_um,
            "minimum_cap_normal_dot": self.minimum_cap_normal_dot,
            "extension_length_error_um": self.extension_length_error_um,
            "extension_axis_dot": self.extension_axis_dot,
            "intermediate_ring_centerline_max_deviation_um": (
                self.intermediate_ring_centerline_max_deviation_um
            ),
            "maximum_intermediate_ring_center_drift_um": (
                self.intermediate_ring_centerline_max_deviation_um
            ),
            "P95_intermediate_ring_center_drift_um": (
                self.intermediate_ring_centerline_p95_deviation_um
            ),
            "mean_intermediate_ring_center_drift_um": (
                self.intermediate_ring_centerline_mean_deviation_um
            ),
            "worst_ring_index": (
                self.intermediate_ring_centerline_worst_ring_index
            ),
            "worst_ring_axial_station_um": (
                self.intermediate_ring_centerline_worst_ring_axial_station_um
            ),
            "intermediate_ring_axial_station_max_error_um": (
                self.intermediate_ring_axial_station_max_error_um
            ),
            "intermediate_ring_area_relative_error_max": (
                self.intermediate_ring_area_relative_error_max
            ),
            "intermediate_ring_all_areas_finite_positive": (
                self.intermediate_ring_all_areas_finite_positive
            ),
            "intermediate_ring_all_polygons_simple_valid": (
                self.intermediate_ring_all_polygons_simple_valid
            ),
            "intermediate_ring_all_orientations_consistent": (
                self.intermediate_ring_all_orientations_consistent
            ),
            "local_original_median_edge_length_um": self.local_original_median_edge_length_um,
            "target_edge_length_um": self.target_edge_length_um,
            "ring_count": self.ring_count,
            "transition_length_um": self.transition_length_um,
            "proximal_ring_max_motion_um": self.proximal_ring_max_motion_um,
            "distal_ring_max_motion_um": self.distal_ring_max_motion_um,
            "cut_loop_vertex_count": self.cut_loop_vertex_count,
            "cap_triangle_count": len(self.cap_face_indices),
            "side_triangle_count": len(self.side_face_indices),
        }
