"""Same-camera visual evidence for v8 stage and k/r comparisons."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh

from .unified_polyball import JunctionBlendSpec


def mesh_polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def local_box_crop(
    mesh: trimesh.Trimesh, center: np.ndarray, half_extent: float
) -> trimesh.Trimesh:
    cropped = mesh.copy()
    for axis in range(3):
        normal = np.zeros(3, dtype=float)
        normal[axis] = 1.0
        cropped = trimesh.intersections.slice_mesh_plane(
            cropped,
            plane_normal=normal,
            plane_origin=center - half_extent * normal,
            cap=False,
        )
        if cropped is None or not len(cropped.faces):
            return trimesh.Trimesh()
        cropped = trimesh.intersections.slice_mesh_plane(
            cropped,
            plane_normal=-normal,
            plane_origin=center + half_extent * normal,
            cap=False,
        )
        if cropped is None or not len(cropped.faces):
            return trimesh.Trimesh()
    cropped.merge_vertices(digits_vertex=8)
    cropped.remove_unreferenced_vertices()
    return cropped


def shared_problem_camera(
    meshes: tuple[trimesh.Trimesh, ...],
) -> list[tuple[float, float, float]]:
    points = np.vstack([np.asarray(mesh.vertices, dtype=float) for mesh in meshes])
    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    center = 0.5 * (minimum + maximum)
    extent = max(float(np.max(maximum - minimum)), 1.0)
    direction = np.asarray((1.35, -1.55, 1.05), dtype=float)
    direction /= np.linalg.norm(direction)
    return [
        tuple(center + 2.8 * extent * direction),
        tuple(center),
        (0.0, 0.0, 1.0),
    ]


def _add_surface(
    plotter: pv.Plotter,
    mesh: trimesh.Trimesh,
    *,
    mode: str,
    color: str,
) -> None:
    data = mesh_polydata(mesh)
    if mode == "silhouette":
        plotter.add_mesh(
            data,
            color="#F8FAFC",
            smooth_shading=False,
            silhouette={"color": "#EF4444", "line_width": 4.0},
        )
    else:
        plotter.add_mesh(
            data,
            color=color,
            smooth_shading=False,
            show_edges=mode == "wireframe",
            edge_color="#111827",
            line_width=0.65,
            specular=0.12,
        )


def generate_stage_localization_figures(
    stage_meshes: dict[str, trimesh.Trimesh],
    spec: JunctionBlendSpec,
    root: Path,
) -> tuple[list[Path], list[tuple[float, float, float]], dict[str, trimesh.Trimesh]]:
    root.mkdir(parents=True, exist_ok=True)
    stage_order = (
        "S0_raw_flying_edges",
        "S1_newton_projected",
        "S2_pyacvd_before_second_projection",
        "S3_final_projected_before_port_clip",
    )
    half_extent = 1.15 * spec.blend_length_um
    local = {
        stage: local_box_crop(
            stage_meshes[stage], np.asarray(spec.center_world_um), half_extent
        )
        for stage in stage_order
    }
    camera = shared_problem_camera(tuple(local[stage] for stage in stage_order))
    paths: list[Path] = []
    colors = ("#38BDF8", "#A78BFA", "#F59E0B", "#22C55E")
    for stage_index, (stage, color) in enumerate(zip(stage_order, colors)):
        plotter = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 820))
        for column, mode in enumerate(("flat", "wireframe")):
            plotter.subplot(0, column)
            plotter.add_text(f"{stage} | {mode}", font_size=12)
            _add_surface(plotter, local[stage], mode=mode, color=color)
            plotter.camera_position = camera
            plotter.add_axes()
        plotter.link_views()
        # Keep filenames short enough for Windows MAX_PATH with long ROI IDs.
        path = root / f"S{stage_index}_flat_wire.png"
        plotter.screenshot(str(path))
        plotter.close()
        paths.append(path)
    for mode in ("flat", "wireframe"):
        plotter = pv.Plotter(shape=(1, 4), off_screen=True, window_size=(2600, 650))
        for column, (stage, color) in enumerate(zip(stage_order, colors)):
            plotter.subplot(0, column)
            plotter.add_text(stage.split("_", 1)[0], font_size=12)
            _add_surface(plotter, local[stage], mode=mode, color=color)
            plotter.camera_position = camera
        plotter.link_views()
        path = root / f"stages_{mode}.png"
        plotter.screenshot(str(path))
        plotter.close()
        paths.append(path)
    return paths, camera, local


def generate_k_sensitivity_figures(
    meshes: dict[str, trimesh.Trimesh],
    spec: JunctionBlendSpec,
    root: Path,
    camera: list[tuple[float, float, float]] | None = None,
) -> tuple[list[Path], list[tuple[float, float, float]], dict[str, trimesh.Trimesh]]:
    root.mkdir(parents=True, exist_ok=True)
    labels = ("v7_hard_min", "v8_k0p10r", "v8_k0p20r", "v8_k0p30r")
    half_extent = 1.15 * spec.blend_length_um
    local = {
        label: local_box_crop(
            meshes[label], np.asarray(spec.center_world_um), half_extent
        )
        for label in labels
    }
    shared_camera = camera or shared_problem_camera(tuple(local[label] for label in labels))
    colors = ("#38BDF8", "#22C55E", "#F59E0B", "#EF4444")
    paths: list[Path] = []
    for mode in ("flat", "wireframe", "silhouette"):
        plotter = pv.Plotter(shape=(1, 4), off_screen=True, window_size=(2600, 650))
        for column, (label, color) in enumerate(zip(labels, colors)):
            plotter.subplot(0, column)
            plotter.add_text(label, font_size=11)
            _add_surface(plotter, local[label], mode=mode, color=color)
            plotter.camera_position = shared_camera
        plotter.link_views()
        path = root / f"k_sensitivity_{mode}.png"
        plotter.screenshot(str(path))
        plotter.close()
        paths.append(path)
    return paths, shared_camera, local


def ownership_and_overlay_figures(
    mesh: trimesh.Trimesh,
    winner: np.ndarray,
    switch_edges: np.ndarray,
    defect_edges: np.ndarray,
    spec: JunctionBlendSpec,
    root: Path,
    camera: list[tuple[float, float, float]],
) -> list[Path]:
    root.mkdir(parents=True, exist_ok=True)
    surface = mesh_polydata(mesh)
    surface.point_data["winner_branch"] = np.asarray(winner, dtype=np.int64)
    paths: list[Path] = []

    plotter = pv.Plotter(off_screen=True, window_size=(1100, 900))
    plotter.add_mesh(
        surface,
        scalars="winner_branch",
        categories=True,
        cmap="tab10",
        smooth_shading=False,
        show_edges=True,
        edge_color="#1F2937",
        line_width=0.4,
    )
    plotter.camera_position = camera
    path = root / "polyball_ownership_winner_branch.png"
    plotter.screenshot(str(path))
    plotter.close()
    paths.append(path)

    def lines(edges: np.ndarray) -> pv.PolyData:
        if not len(edges):
            return pv.PolyData()
        points = np.asarray(mesh.vertices, dtype=float)[edges].reshape((-1, 3))
        connectivity = np.column_stack(
            (
                np.full(len(edges), 2, dtype=np.int64),
                2 * np.arange(len(edges), dtype=np.int64),
                2 * np.arange(len(edges), dtype=np.int64) + 1,
            )
        ).ravel()
        return pv.PolyData(points, lines=connectivity)

    plotter = pv.Plotter(off_screen=True, window_size=(1100, 900))
    plotter.add_mesh(surface, color="#CBD5E1", opacity=0.45, smooth_shading=False)
    if len(defect_edges):
        plotter.add_mesh(lines(defect_edges), color="#F59E0B", line_width=6.0)
    if len(switch_edges):
        plotter.add_mesh(lines(switch_edges), color="#EF4444", line_width=3.0)
    plotter.add_text(
        "orange: saw-tooth/high-dihedral defect | red: ownership switching",
        font_size=11,
    )
    plotter.camera_position = camera
    path = root / "ownership_switch_sawtooth_overlay.png"
    plotter.screenshot(str(path))
    plotter.close()
    paths.append(path)
    return paths


def camera_report(camera: list[tuple[float, float, float]]) -> dict[str, Any]:
    return {
        "position_um": list(camera[0]),
        "focal_point_um": list(camera[1]),
        "view_up": list(camera[2]),
        "projection": "perspective",
        "flat_shading": True,
    }
