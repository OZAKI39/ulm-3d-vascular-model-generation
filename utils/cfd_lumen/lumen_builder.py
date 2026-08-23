"""Variable-radius primitives, balanced manifold union, and implicit fallback."""

from __future__ import annotations

import time
from typing import Iterable

import numpy as np
import pyvista as pv
import trimesh
import vtk
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .types import (
    BranchGeometry,
    GeometryValidationError,
    LumenPrimitives,
    PortGeometry,
    SurfaceBuildResult,
)


def _vtk_polydata_to_trimesh(polydata: vtk.vtkPolyData) -> trimesh.Trimesh:
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputData(polydata)
    triangles.PassLinesOff()
    triangles.PassVertsOff()
    triangles.Update()
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(triangles.GetOutput())
    clean.PointMergingOn()
    clean.Update()
    output = clean.GetOutput()
    wrapped = pv.wrap(output)
    vertices = np.asarray(wrapped.points, dtype=float)
    faces = np.asarray(wrapped.regular_faces, dtype=np.int64)
    if faces.size == 0:
        raise GeometryValidationError("A VTK solid primitive produced no triangle faces")
    return _clean_mesh(trimesh.Trimesh(vertices=vertices, faces=faces, process=False))


def _pyvista_to_trimesh(polydata: pv.PolyData) -> trimesh.Trimesh:
    return _vtk_polydata_to_trimesh(polydata.triangulate())


def _clean_mesh(mesh: trimesh.Trimesh) -> trimesh.Trimesh:
    mesh = mesh.copy()
    # VTK cleaning is applied before conversion and manifold3d already returns
    # indexed manifold output.  A second tolerance-based vertex merge here can
    # collapse intentionally distinct Boolean vertices and introduce pinholes.
    mesh.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(mesh, multibody=True)
    if mesh.is_watertight and mesh.volume < 0:
        mesh.invert()
    return mesh


def build_variable_radius_tube(branch: BranchGeometry, tube_sides: int) -> trimesh.Trimesh:
    points = vtk.vtkPoints()
    for point in branch.points_um:
        points.InsertNextPoint(*map(float, point))
    line = vtk.vtkPolyLine()
    line.GetPointIds().SetNumberOfIds(len(branch.points_um))
    for index in range(len(branch.points_um)):
        line.GetPointIds().SetId(index, index)
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(line)
    radii = vtk.vtkDoubleArray()
    radii.SetName("radius_um")
    radii.SetNumberOfComponents(1)
    for radius in branch.radius_um:
        radii.InsertNextValue(float(radius))
    centerline = vtk.vtkPolyData()
    centerline.SetPoints(points)
    centerline.SetLines(cells)
    centerline.GetPointData().SetScalars(radii)
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(centerline)
    tube.SetNumberOfSides(int(tube_sides))
    tube.SetRadius(1.0)
    tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube.SetCapping(True)
    tube.SidesShareVerticesOn()
    tube.Update()
    mesh = _vtk_polydata_to_trimesh(tube.GetOutput())
    if not mesh.is_watertight:
        raise GeometryValidationError(f"Branch {branch.branch_id} tube primitive is not watertight")
    return mesh


def _junction_mesh(position_um: np.ndarray, radius_um: float, tube_sides: int) -> trimesh.Trimesh:
    sphere = pv.Sphere(
        radius=float(radius_um),
        center=tuple(map(float, position_um)),
        theta_resolution=int(tube_sides),
        phi_resolution=max(8, int(tube_sides) // 2),
    )
    return _pyvista_to_trimesh(sphere)


def _port_extension_mesh(port: PortGeometry, tube_sides: int) -> trimesh.Trimesh:
    vector = port.cylinder_end_um - port.cylinder_start_um
    height = float(np.linalg.norm(vector))
    cylinder = pv.Cylinder(
        center=tuple(map(float, (port.cylinder_start_um + port.cylinder_end_um) * 0.5)),
        direction=tuple(map(float, vector / height)),
        radius=port.radius_um,
        height=height,
        resolution=int(tube_sides),
        capping=True,
    )
    mesh = _pyvista_to_trimesh(cylinder)
    if not mesh.is_watertight:
        raise GeometryValidationError(f"Port {port.port_id} extension primitive is not watertight")
    return mesh


def build_port_extension_mesh(port: PortGeometry, tube_sides: int) -> trimesh.Trimesh:
    """Public explicit port primitive; v3 keeps this geometry unchanged."""

    return _port_extension_mesh(port, tube_sides)


def build_lumen_primitives(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    tube_sides: int | None = None,
) -> LumenPrimitives:
    """Build, but do not union, the exact solids used by the explicit backend."""

    sides = int(tube_sides or config.geometry.tube_sides)
    degree = np.bincount(
        np.asarray(roi.local_edges, dtype=np.int64).ravel(), minlength=len(roi.local_node_ids)
    )
    junction_ids = (
        np.flatnonzero(degree >= 3) if config.junction.enabled else np.empty(0, dtype=np.int64)
    )
    return LumenPrimitives(
        branch_tubes={
            branch.branch_id: build_variable_radius_tube(branch, sides) for branch in branches
        },
        junction_solids={
            int(node_id): _junction_mesh(
                np.asarray(roi.local_node_positions_um[node_id], dtype=float),
                float(roi.local_node_radius_um[node_id]) * config.junction.radius_scale,
                sides,
            )
            for node_id in junction_ids
        },
        port_extensions={
            port.port_id: _port_extension_mesh(port, sides) for port in ports
        },
    )


def balanced_manifold_union(meshes: Iterable[trimesh.Trimesh]) -> trimesh.Trimesh:
    """Union solids pairwise in a balanced tree to limit accumulated complexity."""

    level = [_clean_mesh(mesh) for mesh in meshes]
    if not level:
        raise GeometryValidationError("No lumen primitives were constructed")
    for index, mesh in enumerate(level):
        if not mesh.is_watertight or abs(float(mesh.volume)) <= 1.0e-15:
            raise GeometryValidationError(f"Solid primitive {index} is not a positive watertight volume")
    while len(level) > 1:
        following: list[trimesh.Trimesh] = []
        for index in range(0, len(level), 2):
            if index + 1 >= len(level):
                following.append(level[index])
                continue
            united = trimesh.boolean.union(
                [level[index], level[index + 1]], engine="manifold", check_volume=True
            )
            if not isinstance(united, trimesh.Trimesh) or len(united.faces) == 0:
                raise GeometryValidationError("manifold3d returned an empty/non-mesh Boolean result")
            following.append(_clean_mesh(united))
        level = following
    return _clean_mesh(level[0])


def _capped_cylinder_sdf(points: np.ndarray, port: PortGeometry) -> np.ndarray:
    start = port.cylinder_start_um
    end = port.cylinder_end_um
    center = (start + end) * 0.5
    axis_vector = end - start
    half_length = float(np.linalg.norm(axis_vector)) * 0.5
    axis = axis_vector / (2.0 * half_length)
    relative = points - center
    axial = relative @ axis
    radial_vector = relative - axial[:, None] * axis[None, :]
    radial = np.linalg.norm(radial_vector, axis=1)
    q_radial = radial - port.radius_um
    q_axial = np.abs(axial) - half_length
    outside = np.linalg.norm(np.maximum(np.column_stack((q_radial, q_axial)), 0.0), axis=1)
    inside = np.minimum(np.maximum(q_radial, q_axial), 0.0)
    return outside + inside


def _implicit_lumen(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> tuple[trimesh.Trimesh, dict[str, object]]:
    sample_points = np.concatenate([branch.points_um for branch in branches], axis=0)
    sample_radii = np.concatenate([branch.radius_um for branch in branches], axis=0)
    graph_degree = np.zeros(len(roi.local_node_ids), dtype=np.int64)
    np.add.at(graph_degree, np.asarray(roi.local_edges)[:, 0], 1)
    np.add.at(graph_degree, np.asarray(roi.local_edges)[:, 1], 1)
    junction_ids = np.flatnonzero(graph_degree >= 3)
    if len(junction_ids):
        sample_points = np.vstack((sample_points, roi.local_node_positions_um[junction_ids]))
        sample_radii = np.concatenate((sample_radii, roi.local_node_radius_um[junction_ids]))
    minimum_radius = float(sample_radii.min())
    spacing = float(
        np.clip(
            2.0 * minimum_radius / config.implicit_fallback.cells_across_min_diameter,
            config.implicit_fallback.min_spacing_um,
            config.implicit_fallback.max_spacing_um,
        )
    )
    minimum = np.min(sample_points - sample_radii[:, None], axis=0)
    maximum = np.max(sample_points + sample_radii[:, None], axis=0)
    for port in ports:
        port_min = np.minimum(port.cylinder_start_um, port.cylinder_end_um) - port.radius_um
        port_max = np.maximum(port.cylinder_start_um, port.cylinder_end_um) + port.radius_um
        minimum = np.minimum(minimum, port_min)
        maximum = np.maximum(maximum, port_max)
    minimum -= 2.0 * spacing
    maximum += 2.0 * spacing
    axes = [np.arange(minimum[axis], maximum[axis] + 0.5 * spacing, spacing) for axis in range(3)]
    nx, ny, nz = (len(axis) for axis in axes)
    cell_count = int(nx * ny * nz)
    if cell_count > config.implicit_fallback.max_grid_cells:
        raise GeometryValidationError(
            f"Implicit grid would contain {cell_count:,} cells, exceeding configured limit "
            f"{config.implicit_fallback.max_grid_cells:,}"
        )
    field = np.empty((nz, ny, nx), dtype=np.float32)
    tree = cKDTree(sample_points)
    k_nearest = min(config.implicit_fallback.k_nearest, len(sample_points))
    plane_size = nx * ny
    slab_depth = max(1, config.implicit_fallback.chunk_size // max(1, plane_size))
    x_values, y_values, z_values = axes
    for z_start in range(0, nz, slab_depth):
        z_end = min(nz, z_start + slab_depth)
        zz, yy, xx = np.meshgrid(z_values[z_start:z_end], y_values, x_values, indexing="ij")
        query = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        distances, indices = tree.query(query, k=k_nearest, workers=-1)
        if k_nearest == 1:
            phi = distances - sample_radii[indices]
        else:
            phi = np.min(distances - sample_radii[indices], axis=1)
        for port in ports:
            phi = np.minimum(phi, _capped_cylinder_sdf(query, port))
        field[z_start:z_end] = phi.reshape((z_end - z_start, ny, nx)).astype(np.float32)
    if not (float(field.min()) < 0.0 < float(field.max())):
        raise GeometryValidationError("Implicit field does not bracket the lumen zero level")
    vertices_zyx, faces, _, _ = marching_cubes(field, level=0.0, spacing=(spacing, spacing, spacing))
    vertices_xyz = vertices_zyx[:, ::-1] + minimum
    mesh = _clean_mesh(trimesh.Trimesh(vertices=vertices_xyz, faces=faces, process=False))
    return mesh, {
        "dtype": "float32",
        "spacing_um": spacing,
        "shape_zyx": [nz, ny, nx],
        "cell_count": cell_count,
        "k_nearest": k_nearest,
    }


def build_lumen_surface(
    branches: list[BranchGeometry],
    roi: ROIRecord,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    tube_sides: int | None = None,
    backend: str | None = None,
    controlled_local_implicit: bool = False,
    transition_backend: str | None = None,
    continuous_port_extensions: bool | None = None,
) -> SurfaceBuildResult:
    sides = int(tube_sides or config.geometry.tube_sides)
    requested = backend or config.boolean.backend
    graph_degree = np.zeros(len(roi.local_node_ids), dtype=np.int64)
    np.add.at(graph_degree, np.asarray(roi.local_edges)[:, 0], 1)
    np.add.at(graph_degree, np.asarray(roi.local_edges)[:, 1], 1)
    junction_ids = np.flatnonzero(graph_degree >= 3) if config.junction.enabled else np.empty(0, dtype=int)
    tube_count = len(branches)
    junction_count = int(len(junction_ids))
    continuous_ports = (
        config.port_transition.backend == "continuous_centerline"
        if continuous_port_extensions is None
        else continuous_port_extensions
    )
    extension_count = 0 if continuous_ports and requested != "implicit" else len(ports)
    if requested == "implicit":
        started = time.perf_counter()
        mesh, grid = _implicit_lumen(branches, roi, ports, config)
        return SurfaceBuildResult(
            mesh=mesh,
            backend_requested=requested,
            backend_used="implicit",
            fallback_reason=None,
            tube_primitive_count=tube_count,
            junction_primitive_count=junction_count,
            extension_primitive_count=extension_count,
            boolean_runtime_s=time.perf_counter() - started,
            implicit_grid=grid,
        )
    if config.junction.enabled and config.junction.backend == "local_implicit":
        started = time.perf_counter()
        try:
            # Lazy import avoids a module cycle: hybrid_merge reuses the exact
            # explicit branch/port primitive builders defined above.
            from .hybrid_merge import build_hybrid_lumen

            mesh, details = build_hybrid_lumen(
                branches,
                roi,
                ports,
                config,
                tube_sides=sides,
                controlled_local_implicit=controlled_local_implicit,
                transition_backend=transition_backend,
                continuous_port_extensions=continuous_ports,
            )
            if details.transition_backend == "loop_stitch":
                backend_used = "hybrid_loop_stitch"
            elif details.transition_fallback_reason:
                backend_used = "hybrid_manifold_boolean_fallback"
            else:
                backend_used = (
                    "hybrid_controlled_local_implicit"
                    if controlled_local_implicit
                    else "hybrid_local_implicit"
                )
            return SurfaceBuildResult(
                mesh=mesh,
                backend_requested=requested,
                backend_used=backend_used,
                fallback_reason=None,
                tube_primitive_count=tube_count,
                junction_primitive_count=junction_count,
                extension_primitive_count=extension_count,
                boolean_runtime_s=details.runtime_s["reconstruction_total"],
                hybrid_details=details,
                surface_continuity_version=(
                    "v5"
                    if continuous_ports or details.transition_backend == "loop_stitch"
                    else "v4"
                ),
            )
        except Exception as exc:
            if not config.boolean.allow_implicit_fallback:
                raise
            fallback_reason = f"HybridFailure: {type(exc).__name__}: {exc}"
            mesh, grid = _implicit_lumen(branches, roi, ports, config)
            return SurfaceBuildResult(
                mesh=mesh,
                backend_requested=requested,
                backend_used="implicit",
                fallback_reason=fallback_reason,
                tube_primitive_count=tube_count,
                junction_primitive_count=junction_count,
                extension_primitive_count=extension_count,
                boolean_runtime_s=time.perf_counter() - started,
                implicit_grid=grid,
            )
    started = time.perf_counter()
    try:
        primitives = build_lumen_primitives(branches, roi, ports, config, tube_sides=sides)
        mesh = balanced_manifold_union(primitives.all_meshes)
        return SurfaceBuildResult(
            mesh=mesh,
            backend_requested=requested,
            backend_used="manifold",
            fallback_reason=None,
            tube_primitive_count=tube_count,
            junction_primitive_count=junction_count,
            extension_primitive_count=extension_count,
            boolean_runtime_s=time.perf_counter() - started,
        )
    except Exception as exc:
        if not config.boolean.allow_implicit_fallback:
            raise
        fallback_reason = f"{type(exc).__name__}: {exc}"
        mesh, grid = _implicit_lumen(branches, roi, ports, config)
        return SurfaceBuildResult(
            mesh=mesh,
            backend_requested=requested,
            backend_used="implicit",
            fallback_reason=fallback_reason,
            tube_primitive_count=tube_count,
            junction_primitive_count=junction_count,
            extension_primitive_count=extension_count,
            boolean_runtime_s=time.perf_counter() - started,
            implicit_grid=grid,
        )
