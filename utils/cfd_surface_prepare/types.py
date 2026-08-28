"""Surface records shared by the formal local-cut path."""

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
