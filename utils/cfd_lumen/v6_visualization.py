"""Deterministic v5/v6 continuity visualizations with shared cameras."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import trimesh
from matplotlib.collections import LineCollection

from .surface_qc import _orthogonal_basis


def _trimesh_to_polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def _camera(meshes: list[trimesh.Trimesh]) -> list[tuple[float, float, float]]:
    points = np.vstack([np.asarray(mesh.vertices, dtype=float) for mesh in meshes])
    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    center = 0.5 * (minimum + maximum)
    extent = max(float(np.max(maximum - minimum)), 1.0)
    direction = np.asarray((1.35, -1.55, 1.05), dtype=float)
    direction /= np.linalg.norm(direction)
    position = center + 2.8 * extent * direction
    return [tuple(position), tuple(center), (0.0, 0.0, 1.0)]


def _paired_surface(
    v5: trimesh.Trimesh,
    v6: trimesh.Trimesh,
    path: Path,
    *,
    smooth: bool,
    wireframe: bool,
) -> Path:
    plotter = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1900, 820))
    camera = _camera([v5, v6])
    for column, (name, mesh, color) in enumerate(
        (("v5 loop-stitch", v5, "#38BDF8"), ("v6 continuous field", v6, "#22C55E"))
    ):
        plotter.subplot(0, column)
        plotter.add_text(name, font_size=13)
        plotter.add_mesh(
            _trimesh_to_polydata(mesh),
            color=color,
            smooth_shading=smooth,
            show_edges=wireframe,
            edge_color="#111827",
            line_width=0.7,
            specular=0.18,
        )
        plotter.camera_position = camera
        plotter.add_axes()
    plotter.link_views()
    path.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(path))
    plotter.close()
    return path


def v5_v6_flat_shading(
    v5: trimesh.Trimesh, v6: trimesh.Trimesh, path: Path
) -> Path:
    return _paired_surface(v5, v6, path, smooth=False, wireframe=False)


def v5_v6_smooth_shading(
    v5: trimesh.Trimesh, v6: trimesh.Trimesh, path: Path
) -> Path:
    return _paired_surface(v5, v6, path, smooth=True, wireframe=False)


def v5_v6_wireframe(
    v5: trimesh.Trimesh, v6: trimesh.Trimesh, path: Path
) -> Path:
    return _paired_surface(v5, v6, path, smooth=False, wireframe=True)


def v5_v6_wireframe_overlay(
    v5: trimesh.Trimesh, v6: trimesh.Trimesh, path: Path
) -> Path:
    """Overlay both triangulations in one unchanged camera for seam inspection."""

    plotter = pv.Plotter(off_screen=True, window_size=(1500, 1050))
    plotter.add_mesh(
        _trimesh_to_polydata(v5),
        style="wireframe",
        color="#0284C7",
        opacity=0.55,
        line_width=0.65,
        label="v5 loop-stitch",
    )
    plotter.add_mesh(
        _trimesh_to_polydata(v6),
        style="wireframe",
        color="#16A34A",
        opacity=0.55,
        line_width=0.65,
        label="v6 continuous field",
    )
    plotter.add_text("v5 (blue) / v6 (green) wireframe overlay", font_size=13)
    plotter.camera_position = _camera([v5, v6])
    plotter.add_legend()
    plotter.add_axes()
    path.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(path))
    plotter.close()
    return path


def _edge_polydata(mesh: trimesh.Trimesh, edges: np.ndarray) -> pv.PolyData:
    if not len(edges):
        return pv.PolyData()
    points = np.asarray(mesh.vertices, dtype=float)[np.asarray(edges, dtype=np.int64)].reshape(
        (-1, 3)
    )
    lines = np.column_stack(
        (
            np.full(len(edges), 2, dtype=np.int64),
            np.arange(2 * len(edges), dtype=np.int64).reshape((-1, 2)),
        )
    ).ravel()
    return pv.PolyData(points, lines=lines)


def v5_interface_edges_figure(
    mesh: trimesh.Trimesh,
    interface_edges: np.ndarray,
    path: Path,
) -> Path:
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 1050))
    plotter.add_mesh(
        _trimesh_to_polydata(mesh),
        color="#CBD5E1",
        opacity=0.72,
        smooth_shading=False,
        show_edges=True,
        edge_color="#94A3B8",
        line_width=0.25,
    )
    if len(interface_edges):
        plotter.add_mesh(
            _edge_polydata(mesh, interface_edges),
            color="#DC2626",
            line_width=5.0,
        )
    plotter.add_text("v5 loop-stitch interface edges (red)", font_size=13)
    plotter.camera_position = _camera([mesh])
    plotter.add_axes()
    path.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(path))
    plotter.close()
    return path


def junction_field_transition_figure(
    mesh: trimesh.Trimesh,
    face_region: np.ndarray,
    path: Path,
) -> Path:
    polydata = _trimesh_to_polydata(mesh)
    polydata.cell_data["v6_region"] = np.asarray(face_region, dtype=np.int32)
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 1050))
    plotter.add_mesh(
        polydata,
        scalars="v6_region",
        categories=True,
        cmap=["#0284C7", "#F59E0B", "#64748B", "#22C55E"],
        clim=(1, 4),
        smooth_shading=False,
        show_edges=True,
        edge_color="#334155",
        line_width=0.25,
    )
    plotter.add_text(
        "v6: CORE / quintic TRANSITION / EXPLICIT / PURE_BRANCH", font_size=12
    )
    plotter.camera_position = _camera([mesh])
    plotter.add_axes()
    path.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(path))
    plotter.close()
    return path


def port_tangent_profile_figure(
    diagnostics: list[dict[str, Any]], path: Path
) -> Path:
    labels = [f"P{row['port_id']}" for row in diagnostics]
    x = np.arange(len(labels))
    figure, axis = plt.subplots(figsize=(10, 5.5))
    axis.bar(
        x - 0.19,
        [row["v5_interface_tangent_jump_deg"] for row in diagnostics],
        0.38,
        label="v5 actual interface",
        color="#38BDF8",
    )
    axis.bar(
        x + 0.19,
        [row["v6_interface_tangent_jump_deg"] for row in diagnostics],
        0.38,
        label="v6 actual interface",
        color="#22C55E",
    )
    axis.scatter(
        x,
        [row["v5_tangent_jump_deg"] for row in diagnostics],
        marker="x",
        color="#475569",
        label="weighted-fit diagnostic",
        zorder=3,
    )
    axis.set_xticks(x, labels)
    axis.set(
        ylabel="source/extension tangent jump (deg)",
        title="Port C1 tangent continuity (bars: actual CUT_PORT seam)",
    )
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220)
    plt.close(figure)
    return path


def port_radius_profile_figure(
    diagnostics: list[dict[str, Any]], path: Path
) -> Path:
    figure, axes = plt.subplots(
        len(diagnostics), 1, figsize=(10, max(3.2 * len(diagnostics), 5.5)), squeeze=False
    )
    for axis, row in zip(axes[:, 0], diagnostics):
        distance = row["distance_profile_um"]
        axis.plot(distance, row["radius_profile_before"], "--", color="#38BDF8", label="v5 constant")
        axis.plot(distance, row["radius_profile_after"], color="#22C55E", label="v6 Hermite")
        axis.set(
            ylabel="radius (um)",
            title=(
                f"P{row['port_id']}: source dr/ds={row['source_radius_slope_um_per_um']:.5f}"
            ),
        )
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[-1, 0].set_xlabel("distance into CFD extension (um)")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220)
    plt.close(figure)
    return path


def interface_dihedral_figure(
    rows: list[dict[str, Any]], path: Path, *, title: str
) -> Path:
    selected = [
        row
        for row in rows
        if row.get("junction_node_id", "ALL") == "ALL"
        and row.get("dihedral_max_deg") is not None
    ]
    labels = [f"{row['version']}\n{row['interface']}" for row in selected]
    x = np.arange(len(selected))
    figure, axis = plt.subplots(figsize=(max(10, 1.4 * len(labels)), 6))
    axis.bar(x - 0.22, [row["dihedral_p95_deg"] for row in selected], 0.22, label="P95")
    axis.bar(x, [row["dihedral_p99_deg"] for row in selected], 0.22, label="P99")
    axis.bar(x + 0.22, [row["dihedral_max_deg"] for row in selected], 0.22, label="max")
    axis.set_xticks(x, labels, rotation=25, ha="right", fontsize=8)
    axis.set(ylabel="dihedral angle (deg)", title=title)
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220)
    plt.close(figure)
    return path


def silhouette_comparison_figure(
    v5: trimesh.Trimesh,
    v6: trimesh.Trimesh,
    path: Path,
) -> Path:
    directions = (
        ("x", np.asarray((1.0, 0.0, 0.0))),
        ("y", np.asarray((0.0, 1.0, 0.0))),
        ("z", np.asarray((0.0, 0.0, 1.0))),
        ("xyz", np.asarray((1.0, 1.0, 1.0))),
        ("x-y+z", np.asarray((1.0, -1.0, 1.0))),
        ("-x+y+z", np.asarray((-1.0, 1.0, 1.0))),
    )
    figure, axes = plt.subplots(2, 6, figsize=(22, 7.5))
    for row_index, (version, mesh, color) in enumerate(
        (("v5", v5, "#0284C7"), ("v6", v6, "#16A34A"))
    ):
        adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
        edges = np.asarray(mesh.face_adjacency_edges, dtype=np.int64)
        normals = np.asarray(mesh.face_normals, dtype=float)
        vertices = np.asarray(mesh.vertices, dtype=float)
        for column, (name, direction) in enumerate(directions):
            direction = direction / np.linalg.norm(direction)
            facing = normals @ direction
            first_facing = facing[adjacency[:, 0]]
            second_facing = facing[adjacency[:, 1]]
            epsilon = 1.0e-7
            selected = ((first_facing > epsilon) & (second_facing < -epsilon)) | (
                (first_facing < -epsilon) & (second_facing > epsilon)
            )
            first, second = _orthogonal_basis(direction)
            projected = np.column_stack((vertices @ first, vertices @ second))
            segments = projected[edges[selected]]
            axes[row_index, column].add_collection(
                LineCollection(segments, colors=color, linewidths=0.45)
            )
            axes[row_index, column].autoscale()
            axes[row_index, column].set_aspect("equal")
            axes[row_index, column].axis("off")
            axes[row_index, column].set_title(f"{version} / {name}", fontsize=9)
    figure.suptitle("Six-view silhouette comparison", fontsize=15)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=220)
    plt.close(figure)
    return path
