"""Low-level mesh conversion and quality measurements."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy


@dataclass(slots=True)
class MeshQuality:
    point_count: int
    triangle_count: int
    bounds_um: tuple[float, float, float, float, float, float]
    extent_um: tuple[float, float, float]
    diagonal_um: float
    surface_area_um2: float
    enclosed_volume_um3: float
    connected_component_count: int
    degenerate_triangle_count: int
    duplicate_triangle_count: int
    boundary_edge_count: int
    non_manifold_edge_count: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def polydata_arrays(mesh: vtk.vtkPolyData) -> tuple[np.ndarray, np.ndarray]:
    points_object = mesh.GetPoints()
    if points_object is None:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64)
    points = vtk_to_numpy(points_object.GetData()).astype(np.float64, copy=False)
    polys = mesh.GetPolys()
    if polys is None or polys.GetNumberOfCells() == 0:
        return points, np.empty((0, 3), dtype=np.int64)
    offsets = vtk_to_numpy(polys.GetOffsetsArray())
    if offsets.size and np.any(np.diff(offsets) != 3):
        raise ValueError("Mesh contains non-triangular polygon cells")
    faces = vtk_to_numpy(polys.GetConnectivityArray()).reshape((-1, 3)).astype(np.int64, copy=False)
    return points, faces


def build_polydata(points: np.ndarray, faces: np.ndarray) -> vtk.vtkPolyData:
    points = np.ascontiguousarray(points, dtype=np.float64)
    faces = np.ascontiguousarray(faces, dtype=np.int64).reshape((-1, 3))

    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(numpy_to_vtk(points, deep=True))
    vtk_cells = vtk.vtkCellArray()
    offsets = np.arange(0, 3 * (len(faces) + 1), 3, dtype=np.int64)
    vtk_cells.SetData(
        numpy_to_vtkIdTypeArray(offsets, deep=True),
        numpy_to_vtkIdTypeArray(faces.ravel(), deep=True),
    )
    output = vtk.vtkPolyData()
    output.SetPoints(vtk_points)
    output.SetPolys(vtk_cells)
    return output


def compact_polydata(points: np.ndarray, faces: np.ndarray) -> vtk.vtkPolyData:
    faces = np.asarray(faces, dtype=np.int64).reshape((-1, 3))
    if len(faces) == 0:
        return build_polydata(np.empty((0, 3)), faces)
    used, inverse = np.unique(faces.ravel(), return_inverse=True)
    return build_polydata(np.asarray(points)[used], inverse.reshape((-1, 3)))


def face_areas(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    a = points[faces[:, 0]]
    b = points[faces[:, 1]]
    c = points[faces[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def remove_degenerate_triangles(
    mesh: vtk.vtkPolyData, area_epsilon: float
) -> tuple[vtk.vtkPolyData, int]:
    points, faces = polydata_arrays(mesh)
    if len(faces) == 0:
        return mesh, 0
    repeated_vertex = (
        (faces[:, 0] == faces[:, 1])
        | (faces[:, 1] == faces[:, 2])
        | (faces[:, 2] == faces[:, 0])
    )
    areas = face_areas(points, faces)
    keep = (~repeated_vertex) & np.isfinite(areas) & (areas > area_epsilon)
    removed = int(np.count_nonzero(~keep))
    return compact_polydata(points, faces[keep]), removed


def edge_topology(faces: np.ndarray) -> tuple[int, int]:
    if len(faces) == 0:
        return 0, 0
    edges = np.concatenate(
        (faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]), axis=0
    )
    edges.sort(axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    return int(np.count_nonzero(counts == 1)), int(np.count_nonzero(counts > 2))


def duplicate_triangle_count(faces: np.ndarray) -> int:
    if len(faces) == 0:
        return 0
    canonical = np.sort(faces, axis=1)
    return int(len(canonical) - len(np.unique(canonical, axis=0)))


def _feature_edge_count(mesh: vtk.vtkPolyData, boundary: bool) -> int:
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(mesh)
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.BoundaryEdgesOff()
    edges.NonManifoldEdgesOff()
    if boundary:
        edges.BoundaryEdgesOn()
    else:
        edges.NonManifoldEdgesOn()
    edges.Update()
    return int(edges.GetOutput().GetNumberOfCells())


def _component_count(mesh: vtk.vtkPolyData) -> int:
    connectivity = vtk.vtkConnectivityFilter()
    connectivity.SetInputData(mesh)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.Update()
    return int(connectivity.GetNumberOfExtractedRegions())


def measure_mesh_quality(mesh: vtk.vtkPolyData, area_epsilon: float = 1.0e-10) -> MeshQuality:
    points, faces = polydata_arrays(mesh)
    if len(points) == 0 or len(faces) == 0:
        raise ValueError("Cannot measure an empty surface mesh")

    areas = face_areas(points, faces)
    bounds = tuple(float(value) for value in mesh.GetBounds())
    extent = (
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
    )
    a = points[faces[:, 0]]
    b = points[faces[:, 1]]
    c = points[faces[:, 2]]
    signed_volume = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0

    return MeshQuality(
        point_count=int(len(points)),
        triangle_count=int(len(faces)),
        bounds_um=bounds,
        extent_um=tuple(float(value) for value in extent),
        diagonal_um=float(np.linalg.norm(extent)),
        surface_area_um2=float(np.sum(areas)),
        enclosed_volume_um3=float(abs(signed_volume)),
        connected_component_count=_component_count(mesh),
        degenerate_triangle_count=int(np.count_nonzero(areas <= area_epsilon)),
        duplicate_triangle_count=duplicate_triangle_count(faces),
        boundary_edge_count=_feature_edge_count(mesh, boundary=True),
        non_manifold_edge_count=_feature_edge_count(mesh, boundary=False),
    )

