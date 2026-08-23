"""Data-derived figures required by the v2 root-cause protocol."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import trimesh

from .types import BranchGeometry, LumenPrimitives, PortGeometry
from .visualization import _equal_3d, _save, _surface_collection


def _surface_axis(mesh: trimesh.Trimesh, title: str) -> tuple[plt.Figure, Any]:
    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    _surface_collection(axis, mesh, color="#80CBC4", alpha=0.28)
    _equal_3d(axis, np.asarray(mesh.vertices))
    axis.set_title(title)
    return figure, axis


def _edge_segments(mesh: trimesh.Trimesh, edges: np.ndarray) -> np.ndarray:
    return np.asarray(mesh.vertices)[np.asarray(edges, dtype=np.int64)] if len(edges) else np.empty((0, 2, 3))


def generate_diagnostic_figures(
    mesh: trimesh.Trimesh,
    implicit_mesh: trimesh.Trimesh,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    primitives: LumenPrimitives,
    port_rows: list[dict[str, Any]],
    junction_rows: list[dict[str, Any]],
    defects: dict[str, Any],
    artifacts: dict[str, Any],
    backend_rows: list[dict[str, Any]],
    figures: Path,
) -> list[Path]:
    figures.mkdir(parents=True, exist_ok=True)
    output: list[Path] = []

    figure, axis = _surface_axis(mesh, "Port step overview: post-Boolean surface")
    for port in ports:
        cut = port.original_position_um
        axis.scatter(*cut, s=70, marker="^", color="#F59E0B")
        axis.quiver(*cut, *(port.outward_tangent * max(port.radius_um * 3.0, 1.0)), color="#EF4444")
        axis.text(*cut, f"P{port.port_id}")
    output.append(_save(figure, figures / "01_port_step_overview.png"))

    figure, axes = plt.subplots(1, 2, figsize=(13, 5))
    ids = [row["port_id"] for row in port_rows]
    width = 0.25
    axes[0].bar(np.asarray(ids) - width, [row["source_radius_um"] for row in port_rows], width, label="source")
    axes[0].bar(ids, [row["branch_mesh_radius_um"] for row in port_rows], width, label="branch mesh")
    axes[0].bar(np.asarray(ids) + width, [row["extension_mesh_radius_um"] for row in port_rows], width, label="extension mesh")
    axes[0].set(xlabel="port_id", ylabel="equivalent radius (um)", title="Measured cross-sections")
    axes[0].legend()
    axes[1].bar(ids, [100.0 * row["step_area_relative_error"] for row in port_rows], color="#F97316")
    axes[1].set(xlabel="port_id", ylabel="symmetric area error (%)", title="Geometric step metric")
    output.append(_save(figure, figures / "02_port_cross_sections.png"))

    figure, axis = plt.subplots(figsize=(10, 5))
    labels = ["source cut", "branch mesh", "extension input", "extension mesh"]
    for row in port_rows:
        values = [row["source_radius_um"], row["branch_mesh_radius_um"], row["extension_radius_input_um"], row["extension_mesh_radius_um"]]
        axis.plot(labels, values, marker="o", label=f"port {row['port_id']}")
    axis.set_ylabel("radius (um)")
    axis.set_title("Radius profile across each CUT_PORT interface")
    axis.legend()
    output.append(_save(figure, figures / "03_port_radius_profile.png"))

    figure = plt.figure(figsize=(11, 8))
    axis = figure.add_subplot(111, projection="3d")
    colors = plt.get_cmap("tab10")
    for row, port in zip(port_rows, ports):
        branch_id = int(row["branch_id"])
        _surface_collection(axis, primitives.branch_tubes[branch_id], color="#38BDF8", alpha=0.18, maximum_faces=8000)
        _surface_collection(axis, primitives.port_extensions[port.port_id], color="#FB923C", alpha=0.28, maximum_faces=8000)
        cut = port.original_position_um
        axis.scatter(*cut, color=colors(port.port_id), s=60, marker="^")
        axis.scatter(*port.cylinder_start_um, color="#111827", s=35, marker="o")
        axis.scatter(*port.cap_center_um, color="#EF4444", s=35, marker="s")
        axis.quiver(*cut, *(port.outward_tangent * max(2.0 * port.radius_um, 1.0)), color="#DC2626")
    all_points = np.vstack([mesh.vertices, *[port.cap_center_um for port in ports]])
    _equal_3d(axis, all_points)
    axis.set_title("Pre-Boolean branch (blue) / extension (orange) overlap")
    output.append(_save(figure, figures / "04_extension_overlap.png"))

    figure = plt.figure(figsize=(11, 8))
    axis = figure.add_subplot(111, projection="3d")
    palette = plt.get_cmap("tab20")
    for branch in branches:
        _surface_collection(axis, primitives.branch_tubes[branch.branch_id], color=palette(branch.branch_id), alpha=0.18, maximum_faces=5000)
        axis.plot(*branch.points_um.T, linewidth=1.6, color=palette(branch.branch_id))
        midpoint = branch.points_um[len(branch.points_um) // 2]
        axis.text(*midpoint, f"B{branch.branch_id}", fontsize=8)
    for node_id, solid in primitives.junction_solids.items():
        _surface_collection(axis, solid, color="#E11D48", alpha=0.38, maximum_faces=5000)
        center = solid.centroid
        axis.scatter(*center, color="#7F1D1D", s=50)
        axis.text(*center, f"J{node_id}")
    _equal_3d(axis, np.asarray(mesh.vertices))
    axis.set_title("Junction primitives before Boolean")
    output.append(_save(figure, figures / "05_junction_primitives.png"))

    figure, axis = _surface_axis(mesh, "Junction neighborhoods after Boolean")
    seen: set[int] = set()
    for row in junction_rows:
        node_id = int(row["junction_node_id"])
        if node_id in seen:
            continue
        seen.add(node_id)
        point = np.asarray((row["junction_x_um"], row["junction_y_um"], row["junction_z_um"]))
        axis.scatter(*point, color="#E11D48", s=75)
        axis.text(*point, f"J{node_id}")
    output.append(_save(figure, figures / "06_junction_post_boolean.png"))

    for name, title, filename, color in (
        ("boundary_edges", "Boundary edges only", "07_boundary_edges.png", "#EF4444"),
        ("nonmanifold_edges", "Non-manifold edges only", "08_nonmanifold_edges.png", "#A855F7"),
    ):
        figure, axis = _surface_axis(mesh, title)
        for segment in _edge_segments(mesh, artifacts[name]):
            axis.plot(*segment.T, linewidth=3.0, color=color)
        axis.text2D(0.02, 0.95, f"count = {len(artifacts[name])}", transform=axis.transAxes)
        output.append(_save(figure, figures / filename))

    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    for index, component in enumerate(mesh.split(only_watertight=False)):
        _surface_collection(axis, component, color=palette(index), alpha=0.38)
    _equal_3d(axis, np.asarray(mesh.vertices))
    axis.set_title(f"Surface connected components: {defects['surface_connected_component_count']}")
    output.append(_save(figure, figures / "09_surface_components.png"))

    faces = np.column_stack((np.full(len(mesh.faces), 3, dtype=np.int64), mesh.faces)).ravel()
    polydata = pv.PolyData(np.asarray(mesh.vertices), faces)
    plotter = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(1800, 600))
    for column, (smooth, wireframe, title) in enumerate(
        ((False, False, "flat shading"), (True, False, "smooth shading"), (False, True, "wireframe overlay"))
    ):
        plotter.subplot(0, column)
        plotter.add_mesh(polydata, color="#7DD3FC", smooth_shading=smooth, opacity=1.0)
        if wireframe:
            plotter.add_mesh(polydata, style="wireframe", color="#111827", line_width=1.0)
        plotter.add_text(title, font_size=11)
        plotter.view_isometric()
    normal_path = figures / "10_normals_overlay.png"
    plotter.show(screenshot=normal_path, auto_close=True)
    output.append(normal_path)

    repaired = mesh.copy()
    trimesh.repair.fix_normals(repaired, multibody=True)
    figure, axis = _surface_axis(repaired, "Recomputed normals preview (formal mesh unchanged)")
    sample = np.linspace(0, len(repaired.faces) - 1, min(250, len(repaired.faces))).astype(int)
    centers = repaired.triangles_center[sample]
    normals = repaired.face_normals[sample]
    scale = max(float(np.ptp(repaired.vertices, axis=0).max()) * 0.015, 0.1)
    axis.quiver(centers[:, 0], centers[:, 1], centers[:, 2], normals[:, 0], normals[:, 1], normals[:, 2], length=scale, color="#DC2626")
    output.append(_save(figure, figures / "10_recomputed_normals_preview.png"))

    aspect = np.asarray(artifacts["aspect_ratio"])
    finite = aspect[np.isfinite(aspect)]
    threshold = float(np.percentile(finite, 95)) if len(finite) else float("inf")
    bad = np.asarray(artifacts.get("junction_bad_faces", np.flatnonzero(aspect >= threshold)))
    figure, axis = _surface_axis(mesh, "Junction triangle quality: aspect >= 20 or min angle <= 1 deg")
    _surface_collection(axis, mesh, face_indices=bad, color="#EF4444", alpha=0.9, maximum_faces=10000)
    output.append(_save(figure, figures / "11_triangle_quality.png"))

    figure, axis = plt.subplots(figsize=(12, 7))
    axis.axis("off")
    lines = [
        "V2 GEOMETRY DEFECT SUMMARY",
        "",
        f"watertight: {defects['watertight']}",
        f"boundary edges: {defects['boundary_edge_count']}",
        f"non-manifold edges: {defects['non_manifold_edge_count']}",
        f"components: {defects['surface_connected_component_count']}",
        f"self-intersection pairs: {defects['self_intersection_count']}",
        f"suspected internal faces: {defects['suspected_internal_face_count']}",
        f"duplicate faces: {defects['duplicate_face_count']}",
        f"flipped faces: {defects['number_of_flipped_faces']}",
        "",
        f"max port endpoint position error (um): {max((r['endpoint_position_error_um'] for r in port_rows), default=float('nan')):.6g}",
        f"max port area-step error: {max((r['step_area_relative_error'] for r in port_rows), default=float('nan')):.6g}",
    ]
    axis.text(0.03, 0.96, "\n".join(lines), va="top", family="monospace", fontsize=14)
    output.append(_save(figure, figures / "12_defect_summary.png"))

    figure = plt.figure(figsize=(15, 7))
    for index, (candidate, title) in enumerate(((mesh, "explicit / manifold"), (implicit_mesh, "implicit fallback")), start=1):
        axis = figure.add_subplot(1, 2, index, projection="3d")
        _surface_collection(axis, candidate, color="#67E8F9" if index == 1 else "#C4B5FD", alpha=0.35)
        _equal_3d(axis, np.asarray(candidate.vertices))
        row = backend_rows[index - 1]
        axis.set_title(f"{title}\ntri={row['triangle_count']}, radius P95={row['radius_p95_absolute_relative_error']:.4g}")
    output.append(_save(figure, figures / "12_explicit_vs_implicit.png"))
    return output
