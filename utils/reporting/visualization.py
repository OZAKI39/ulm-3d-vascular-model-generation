"""Headless PNG visualizations used for functional acceptance."""

from __future__ import annotations

from pathlib import Path
import shutil
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import vtk
from matplotlib.colors import ListedColormap
from scipy import ndimage

from ..config import VisualizationConfig
from ..mesh.cleanup import ComponentRecord
from .acceptance import AcceptanceResult


def _camera_position(bounds: Iterable[float]) -> list[tuple[float, float, float]]:
    values = tuple(float(value) for value in bounds)
    center = np.asarray(
        [(values[0] + values[1]) / 2, (values[2] + values[3]) / 2, (values[4] + values[5]) / 2]
    )
    diagonal = float(
        np.linalg.norm([values[1] - values[0], values[3] - values[2], values[5] - values[4]])
    )
    position = center + diagonal * np.asarray([1.35, 1.35, 1.0])
    return [tuple(position), tuple(center), (0.0, 0.0, 1.0)]


def _display_mesh(mesh: vtk.vtkPolyData, config: VisualizationConfig) -> pv.PolyData:
    wrapped = pv.wrap(mesh)
    face_count = wrapped.n_cells
    if face_count <= config.max_display_faces:
        return wrapped
    reduction = 1.0 - config.max_display_faces / face_count
    try:
        return wrapped.decimate_pro(reduction, preserve_topology=True)
    except Exception:
        return wrapped


def _render_single_mesh(
    mesh: vtk.vtkPolyData,
    path: Path,
    color: str,
    camera: list[tuple[float, float, float]],
    config: VisualizationConfig,
) -> None:
    plotter = pv.Plotter(off_screen=True, window_size=(config.window_width, config.window_height))
    plotter.set_background("white")
    plotter.add_mesh(_display_mesh(mesh, config), color=color, smooth_shading=True)
    plotter.camera_position = camera
    plotter.show(screenshot=str(path), auto_close=True)


def render_mesh_visualizations(
    before: vtk.vtkPolyData,
    after: vtk.vtkPolyData,
    removed: vtk.vtkPolyData | None,
    removed_small_fragments: vtk.vtkPolyData | None,
    removed_island_networks: vtk.vtkPolyData | None,
    components: list[ComponentRecord],
    output_dir: Path,
    config: VisualizationConfig,
) -> list[Path]:
    if not config.enabled:
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    camera = _camera_position(before.GetBounds())
    paths: list[Path] = []

    before_path = output_dir / "mesh_before.png"
    after_path = output_dir / "mesh_after.png"
    _render_single_mesh(before, before_path, "#397a5d", camera, config)
    _render_single_mesh(after, after_path, "#376da3", camera, config)
    paths.extend((before_path, after_path))
    main_only_path = output_dir / "main_network_only.png"
    shutil.copyfile(after_path, main_only_path)
    paths.append(main_only_path)

    compare_path = output_dir / "mesh_before_after.png"
    plotter = pv.Plotter(
        shape=(1, 2), off_screen=True, window_size=(config.window_width * 2, config.window_height)
    )
    plotter.set_background("white")
    plotter.subplot(0, 0)
    plotter.add_text("Before cleanup", font_size=14, color="black")
    plotter.add_mesh(_display_mesh(before, config), color="#397a5d", smooth_shading=True)
    plotter.camera_position = camera
    plotter.subplot(0, 1)
    plotter.add_text("After cleanup", font_size=14, color="black")
    plotter.add_mesh(_display_mesh(after, config), color="#376da3", smooth_shading=True)
    plotter.camera_position = camera
    plotter.show(screenshot=str(compare_path), auto_close=True)
    paths.append(compare_path)
    connectivity_compare_path = output_dir / "mesh_connectivity_before_after.png"
    shutil.copyfile(compare_path, connectivity_compare_path)
    paths.append(connectivity_compare_path)

    decisions_path = output_dir / "component_decisions.png"
    plotter = pv.Plotter(off_screen=True, window_size=(config.window_width, config.window_height))
    plotter.set_background("white")
    plotter.add_mesh(
        _display_mesh(after, config),
        color="#2d936c",
        smooth_shading=True,
        label="Main network kept",
    )
    if removed_island_networks is not None and removed_island_networks.GetNumberOfCells():
        island_display = _display_mesh(removed_island_networks, config)
        plotter.add_mesh(
            island_display,
            color="#e69f24",
            smooth_shading=True,
            label="Disconnected island networks",
        )
    if removed_small_fragments is not None and removed_small_fragments.GetNumberOfCells():
        fragment_display = _display_mesh(removed_small_fragments, config)
        plotter.add_mesh(fragment_display, color="#d64545", smooth_shading=False)
        removed_markers = pv.PolyData(fragment_display.cell_centers().points)
        plotter.add_mesh(
            removed_markers,
            color="#e22d2d",
            point_size=6.0,
            render_points_as_spheres=True,
            label="Small fragments (enlarged markers)",
        )
    if removed is not None and removed.GetNumberOfCells():
        plotter.add_legend(bcolor="white", face=None)
        plotter.add_text(
            "Orange surfaces are disconnected islands. Red markers enlarge tiny fragments.",
            position="lower_left",
            font_size=10,
            color="black",
        )
    plotter.camera_position = camera
    plotter.show(screenshot=str(decisions_path), auto_close=True)
    paths.append(decisions_path)

    distribution_path = output_dir / "component_distribution.png"
    areas = np.asarray([item.surface_area_um2 for item in components])
    diagonals = np.asarray([item.bbox_diagonal_um for item in components])
    main_mask = np.asarray([item.component_type == "main_network" for item in components])
    fragment_mask = np.asarray([item.component_type == "small_fragment" for item in components])
    island_mask = np.asarray([item.component_type == "island_network" for item in components])
    secondary_mask = ~(main_mask | fragment_mask | island_mask)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    positive = areas[areas > 0]
    axes[0].hist(positive, bins=np.logspace(np.log10(positive.min()), np.log10(positive.max()), 50))
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Component surface area (um^2)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Connected-component size distribution")
    axes[1].scatter(
        diagonals[main_mask], areas[main_mask], s=55, c="#2d936c", label="Main network"
    )
    axes[1].scatter(
        diagonals[island_mask], areas[island_mask], s=13, c="#e69f24", label="Island network"
    )
    axes[1].scatter(
        diagonals[fragment_mask], areas[fragment_mask], s=13, c="#d64545", label="Small fragment"
    )
    if np.any(secondary_mask):
        axes[1].scatter(
            diagonals[secondary_mask],
            areas[secondary_mask],
            s=13,
            c="#4f81bd",
            label="Retained secondary",
        )
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Bounding-box diagonal (um)")
    axes[1].set_ylabel("Surface area (um^2)")
    axes[1].legend()
    axes[1].set_title("Main-network selection decisions")
    fig.savefig(distribution_path, dpi=170)
    plt.close(fig)
    paths.append(distribution_path)
    return paths


def _best_slices(data: np.ndarray) -> list[np.ndarray]:
    indices = [
        int(np.argmax(data.sum(axis=(1, 2)))),
        int(np.argmax(data.sum(axis=(0, 2)))),
        int(np.argmax(data.sum(axis=(0, 1)))),
    ]
    return [data[indices[0], :, :].T, data[:, indices[1], :].T, data[:, :, indices[2]].T]


def _mips(data: np.ndarray) -> list[np.ndarray]:
    return [data.max(axis=0).T, data.max(axis=1).T, data.max(axis=2).T]


def render_voxel_visualizations(mask: np.ndarray, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    labels = ("X view", "Y view", "Z view")
    output: list[Path] = []
    for filename, views, title in (
        ("voxel_mask_views.png", _best_slices(mask), "Densest mask slices"),
        ("voxel_mask_mip.png", _mips(mask), "Mask maximum-intensity projections"),
    ):
        path = output_dir / filename
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
        for axis, image, label in zip(axes, views, labels):
            axis.imshow(image, cmap="gray", origin="lower", interpolation="nearest")
            axis.set_title(label)
            axis.axis("off")
        fig.suptitle(title)
        fig.savefig(path, dpi=170)
        plt.close(fig)
        output.append(path)
    return output


def render_voxel_connectivity_visualization(
    final_mask: np.ndarray,
    removed_islands: np.ndarray,
    output_dir: Path,
) -> Path:
    before_mask = final_mask | removed_islands
    path = output_dir / "voxel_connectivity_before_after.png"
    removed_color_map = ListedColormap(["black", "#f0a020"])
    rows = (
        ("Before voxel-island filtering", _mips(before_mask), "gray"),
        ("Removed voxel islands", _mips(removed_islands), removed_color_map),
        ("Final main voxel network", _mips(final_mask), "gray"),
    )
    fig, axes = plt.subplots(3, 3, figsize=(15, 13), constrained_layout=True)
    for row_index, (row_title, views, color_map) in enumerate(rows):
        for column_index, (image, label) in enumerate(
            zip(views, ("X view", "Y view", "Z view"))
        ):
            axes[row_index, column_index].imshow(
                image, cmap=color_map, origin="lower", interpolation="nearest"
            )
            axes[row_index, column_index].set_title(f"{row_title} - {label}")
            axes[row_index, column_index].axis("off")
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def render_skeleton_visualizations(
    mask: np.ndarray,
    skeleton: np.ndarray,
    origin_lps_um: tuple[float, float, float],
    spacing_um: tuple[float, float, float],
    output_dir: Path,
    config: VisualizationConfig,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    overlay_path = output_dir / "skeleton_overlay.png"
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for axis, mask_image, skeleton_image, label in zip(
        axes, _mips(mask), _mips(skeleton), ("X view", "Y view", "Z view")
    ):
        axis.imshow(mask_image, cmap="gray", origin="lower", interpolation="nearest")
        overlay = np.ma.masked_where(~skeleton_image.astype(bool), skeleton_image)
        axis.imshow(overlay, cmap="autumn", origin="lower", interpolation="nearest", alpha=0.9)
        axis.set_title(label)
        axis.axis("off")
    fig.suptitle("Coarse skeleton (red) over voxel mask")
    fig.savefig(overlay_path, dpi=170)
    component_count = ndimage.label(
        skeleton, structure=ndimage.generate_binary_structure(3, 3)
    )[1]
    connectivity_path = output_dir / "skeleton_connectivity.png"
    fig.suptitle(f"Final coarse skeleton: {component_count} connected component(s)")
    fig.savefig(connectivity_path, dpi=170)
    plt.close(fig)

    points_path = output_dir / "skeleton_3d.png"
    indices = np.argwhere(skeleton)
    if len(indices) > config.max_skeleton_points:
        sample = np.linspace(0, len(indices) - 1, config.max_skeleton_points, dtype=np.int64)
        indices = indices[sample]
    coordinates = np.asarray(origin_lps_um) + indices * np.asarray(spacing_um)
    cloud = pv.PolyData(coordinates)
    plotter = pv.Plotter(off_screen=True, window_size=(config.window_width, config.window_height))
    plotter.set_background("white")
    plotter.add_mesh(cloud, color="#c93434", point_size=2.5, render_points_as_spheres=False)
    mask_bounds = (
        origin_lps_um[0],
        origin_lps_um[0] + spacing_um[0] * (mask.shape[0] - 1),
        origin_lps_um[1],
        origin_lps_um[1] + spacing_um[1] * (mask.shape[1] - 1),
        origin_lps_um[2],
        origin_lps_um[2] + spacing_um[2] * (mask.shape[2] - 1),
    )
    plotter.camera_position = _camera_position(mask_bounds)
    plotter.show(screenshot=str(points_path), auto_close=True)
    return [overlay_path, connectivity_path, points_path]


def render_acceptance_dashboard(result: AcceptanceResult, path: Path) -> Path:
    colors = {"PASS": "#2d936c", "WARNING": "#d89b26", "FAIL": "#d64545"}
    labels = [item.name for item in result.checks]
    values = [1] * len(labels)
    bar_colors = [colors[item.status] for item in result.checks]
    height = max(6.0, 0.46 * len(labels) + 1.8)
    fig, axis = plt.subplots(figsize=(12, height), constrained_layout=True)
    y = np.arange(len(labels))
    axis.barh(y, values, color=bar_colors, height=0.62)
    axis.set_yticks(y, labels=labels)
    axis.set_xlim(0, 1.05)
    axis.set_xticks([])
    axis.invert_yaxis()
    for index, item in enumerate(result.checks):
        axis.text(0.03, index, item.status, va="center", ha="left", color="white", weight="bold")
    axis.set_title(f"Automatic acceptance status: {result.overall_status}", fontsize=16)
    for side in ("top", "right", "bottom", "left"):
        axis.spines[side].set_visible(False)
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path
