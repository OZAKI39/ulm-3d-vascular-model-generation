"""Same-camera v6/v7 visual evidence required by the v7 protocol."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyvista as pv
import trimesh


def _polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def _camera(meshes: tuple[trimesh.Trimesh, ...]):
    points = np.vstack([np.asarray(mesh.vertices, dtype=float) for mesh in meshes])
    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    center = 0.5 * (minimum + maximum)
    extent = max(float(np.max(maximum - minimum)), 1.0)
    direction = np.asarray((1.35, -1.55, 1.05), dtype=float)
    direction /= np.linalg.norm(direction)
    return [tuple(center + 2.8 * extent * direction), tuple(center), (0.0, 0.0, 1.0)]


def paired_surface_figure(
    v6_mesh: trimesh.Trimesh,
    v7_mesh: trimesh.Trimesh,
    path: Path,
    *,
    mode: str,
) -> Path:
    plotter = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1900, 820))
    camera = _camera((v6_mesh, v7_mesh))
    for column, (title, mesh, color) in enumerate(
        (("v6 topology-valid hybrid", v6_mesh, "#38BDF8"), ("v7 unified PolyBall", v7_mesh, "#22C55E"))
    ):
        plotter.subplot(0, column)
        plotter.add_text(title, font_size=13)
        if mode == "silhouette":
            plotter.add_mesh(
                _polydata(mesh),
                color="#F8FAFC",
                smooth_shading=True,
                silhouette={"color": color, "line_width": 4.0},
            )
        else:
            plotter.add_mesh(
                _polydata(mesh),
                color=color,
                smooth_shading=mode == "smooth",
                show_edges=mode == "wireframe",
                edge_color="#111827",
                line_width=0.6,
                specular=0.18,
            )
        plotter.camera_position = camera
        plotter.add_axes()
    plotter.link_views()
    path.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(path))
    plotter.close()
    return path


def generate_v7_comparison_figures(
    v6_mesh: trimesh.Trimesh,
    v7_mesh: trimesh.Trimesh,
    root: Path,
) -> list[Path]:
    return [
        paired_surface_figure(v6_mesh, v7_mesh, root / "v6_v7_flat_same_camera.png", mode="flat"),
        paired_surface_figure(v6_mesh, v7_mesh, root / "v6_v7_smooth_same_camera.png", mode="smooth"),
        paired_surface_figure(v6_mesh, v7_mesh, root / "v6_v7_wireframe_same_camera.png", mode="wireframe"),
        paired_surface_figure(v6_mesh, v7_mesh, root / "v6_v7_silhouette_same_camera.png", mode="silhouette"),
    ]
