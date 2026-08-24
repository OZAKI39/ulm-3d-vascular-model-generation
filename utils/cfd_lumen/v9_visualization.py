"""Same-camera v7/v8/v9 visual comparisons."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import trimesh

from .types import BranchGeometry
from .unified_polyball import JunctionBlendSpec
from .v8_visualization import local_box_crop, mesh_polydata, shared_problem_camera


def _save(plotter: pv.Plotter, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(path)
    plotter.close()
    return path


def generate_method_comparison(
    meshes: dict[str, trimesh.Trimesh],
    spec: JunctionBlendSpec,
    root: Path,
    *,
    camera: list[tuple[float, float, float]] | None = None,
) -> tuple[list[Path], list[tuple[float, float, float]], dict[str, trimesh.Trimesh]]:
    root.mkdir(parents=True, exist_ok=True)
    local = {
        label: local_box_crop(
            mesh,
            np.asarray(spec.center_world_um),
            1.25 * spec.blend_length_um,
        )
        for label, mesh in meshes.items()
    }
    selected_camera = camera or shared_problem_camera(tuple(local.values()))
    labels = list(local)
    paths: list[Path] = []
    for mode in ("flat", "wireframe", "silhouette"):
        plotter = pv.Plotter(
            shape=(1, len(labels)),
            off_screen=True,
            window_size=(520 * len(labels), 520),
            border=False,
        )
        for column, label in enumerate(labels):
            plotter.subplot(0, column)
            data = mesh_polydata(local[label])
            if mode == "silhouette":
                plotter.set_background("white")
                plotter.add_mesh(
                    data,
                    color="black",
                    smooth_shading=False,
                    show_edges=False,
                    lighting=False,
                )
            else:
                plotter.set_background("#101216")
                plotter.add_mesh(
                    data,
                    color="#87c9ff",
                    smooth_shading=False,
                    show_edges=mode == "wireframe",
                    edge_color="#18222b",
                    line_width=0.6,
                )
            plotter.add_text(label, font_size=10, color="white" if mode != "silhouette" else "black")
            plotter.camera_position = selected_camera
        plotter.link_views()
        paths.append(_save(plotter, root / f"methods_{mode}.png"))
    return paths, selected_camera, local


def centerline_fidelity_figure(
    source: list[BranchGeometry],
    smooth: list[BranchGeometry],
    path: Path,
) -> Path:
    source_by_id = {branch.branch_id: branch for branch in source}
    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    for branch in smooth:
        raw = source_by_id[branch.branch_id].raw_points_um
        axis.plot(
            raw[:, 0], raw[:, 1], raw[:, 2], "o--", color="#ef7d32", ms=3, lw=0.8
        )
        axis.plot(
            branch.points_um[:, 0],
            branch.points_um[:, 1],
            branch.points_um[:, 2],
            color="#2a9d8f",
            lw=1.4,
        )
    axis.set_xlabel("X (um)")
    axis.set_ylabel("Y (um)")
    axis.set_zlabel("Z (um)")
    axis.set_title("Source SWC constraints (orange) vs CFD-derived C1 spline (green)")
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def tangent_defect_correlation_figure(
    rows: list[dict[str, Any]], path: Path
) -> Path:
    matched = [
        row
        for row in rows
        if row.get("source_tangent_kink_deg") is not None
    ]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    if matched:
        kink = np.asarray([row["source_tangent_kink_deg"] for row in matched])
        gradient = np.asarray([row["gradient_jump_angle_deg"] for row in matched])
        dihedral = np.asarray([row["surface_dihedral_angle_deg"] for row in matched])
        defect = np.asarray([bool(row["is_sawtooth_defect"]) for row in matched])
        colors = np.where(defect, "#d62828", "#457b9d")
        axes[0].scatter(kink, gradient, c=colors, s=14, alpha=0.75)
        axes[1].scatter(kink, dihedral, c=colors, s=14, alpha=0.75)
    axes[0].set(xlabel="source tangent kink (deg)", ylabel="field gradient jump (deg)")
    axes[1].set(xlabel="source tangent kink (deg)", ylabel="surface dihedral (deg)")
    for axis in axes:
        axis.grid(alpha=0.2)
    figure.suptitle("Segment switch correlation (red = saw-tooth defect)")
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path
