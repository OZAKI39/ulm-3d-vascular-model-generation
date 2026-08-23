"""Headless, data-derived acceptance visualizations for CFD lumen geometry."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import trimesh
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from utils.sampling.sampling_types import ROIRecord

from .types import (
    BranchGeometry,
    CollisionEvent,
    PatchResult,
    PortGeometry,
    RadiusFidelitySample,
)


def _equal_3d(axis: object, points: np.ndarray) -> None:
    if len(points) == 0:
        return
    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    center = (minimum + maximum) * 0.5
    radius = max(float((maximum - minimum).max()) * 0.55, 1.0)
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_xlabel("X (um)")
    axis.set_ylabel("Y (um)")
    axis.set_zlabel("Z (um)")


def _save(figure: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return path


def _surface_collection(
    axis: object,
    mesh: trimesh.Trimesh,
    *,
    face_indices: np.ndarray | None = None,
    color: str = "#4FC3F7",
    alpha: float = 0.25,
    edgecolor: str = "none",
    maximum_faces: int = 25_000,
) -> None:
    indices = np.arange(len(mesh.faces)) if face_indices is None else np.asarray(face_indices, dtype=int)
    if len(indices) > maximum_faces:
        indices = indices[np.linspace(0, len(indices) - 1, maximum_faces).astype(int)]
    triangles = mesh.vertices[mesh.faces[indices]]
    collection = Poly3DCollection(
        triangles,
        facecolor=color,
        edgecolor=edgecolor,
        alpha=alpha,
        linewidth=0.15,
    )
    axis.add_collection3d(collection)


def source_geometry_figure(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    path: Path,
) -> Path:
    figure = plt.figure(figsize=(9, 7))
    axis = figure.add_subplot(111, projection="3d")
    all_radii = np.concatenate([branch.raw_radius_um for branch in branches])
    norm = plt.Normalize(float(all_radii.min()), float(all_radii.max()))
    cmap = plt.get_cmap("viridis")
    for branch in branches:
        segments = np.stack((branch.raw_points_um[:-1], branch.raw_points_um[1:]), axis=1)
        radii = (branch.raw_radius_um[:-1] + branch.raw_radius_um[1:]) * 0.5
        axis.add_collection3d(Line3DCollection(segments, colors=cmap(norm(radii)), linewidths=2.0))
    if roi.cut_ports:
        cut = np.asarray([port.intersection_position_um for port in roi.cut_ports], dtype=float)
        axis.scatter(cut[:, 0], cut[:, 1], cut[:, 2], marker="^", s=55, color="#F59E0B", label="CUT_PORT")
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    figure.colorbar(mappable, ax=axis, shrink=0.65, label="source radius (um)")
    axis.set_title(f"Source SWC ROI geometry: {roi.roi_id}")
    axis.legend(loc="upper right")
    _equal_3d(axis, roi.local_node_positions_um)
    return _save(figure, path)


def resampling_figure(branches: list[BranchGeometry], path: Path) -> Path:
    figure = plt.figure(figsize=(9, 7))
    axis = figure.add_subplot(111, projection="3d")
    for index, branch in enumerate(branches):
        axis.plot(*branch.raw_points_um.T, color="#9CA3AF", linewidth=1.2, alpha=0.9, label="raw" if index == 0 else None)
        axis.plot(*branch.points_um.T, color="#EF4444", linewidth=0.8, alpha=0.8, label="resampled" if index == 0 else None)
        axis.scatter(*branch.points_um.T, color="#EF4444", s=3, alpha=0.45)
    all_points = np.concatenate([branch.raw_points_um for branch in branches])
    axis.set_title("Raw vs arc-length-resampled centerline")
    axis.legend()
    _equal_3d(axis, all_points)
    return _save(figure, path)


def lumen_overlay_figure(
    mesh: trimesh.Trimesh,
    roi: ROIRecord,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    path: Path,
) -> Path:
    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    _surface_collection(axis, mesh, alpha=0.22)
    for branch in branches:
        axis.plot(*branch.points_um.T, color="#0F172A", linewidth=1.0)
    degree = np.bincount(np.asarray(roi.local_edges).ravel(), minlength=roi.node_count)
    junctions = roi.local_node_positions_um[degree >= 3]
    if len(junctions):
        axis.scatter(*junctions.T, color="#DC2626", s=25, label="bifurcation")
    for index, port in enumerate(ports):
        extension = np.vstack((port.cylinder_start_um, port.cylinder_end_um))
        axis.plot(*extension.T, color="#F59E0B", linewidth=2.2, label="CUT_PORT extension" if index == 0 else None)
    axis.set_title("Source centerline + reconstructed CFD lumen")
    if len(junctions) or ports:
        axis.legend()
    _equal_3d(axis, mesh.vertices)
    return _save(figure, path)


def ports_figure(
    mesh: trimesh.Trimesh,
    patch: PatchResult,
    ports: list[PortGeometry],
    path: Path,
) -> Path:
    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    _surface_collection(axis, mesh, alpha=0.08, color="#94A3B8")
    palette = plt.get_cmap("tab20")
    for port in ports:
        faces = np.flatnonzero(patch.port_id == port.port_id)
        _surface_collection(axis, mesh, face_indices=faces, color=palette(port.port_id % 20), alpha=0.95)
        arrow_length = max(3.0 * port.radius_um, 1.0)
        axis.quiver(
            *port.cap_center_um,
            *(port.outward_tangent * arrow_length),
            color="#DC2626",
            arrow_length_ratio=0.25,
            linewidth=1.5,
        )
        axis.scatter(*port.original_position_um, color="#F59E0B", marker="^", s=35)
        axis.text(*port.cap_center_um, f" P{port.port_id}", fontsize=8)
    axis.set_title("Flat CFD port patches and geometric outward normals")
    _equal_3d(axis, mesh.vertices)
    return _save(figure, path)


def wireframe_figure(mesh: trimesh.Trimesh, path: Path) -> Path:
    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    _surface_collection(
        axis,
        mesh,
        alpha=0.10,
        color="#E2E8F0",
        edgecolor="#334155",
        maximum_faces=18_000,
    )
    axis.set_title("CFD lumen triangulated surface wireframe")
    _equal_3d(axis, mesh.vertices)
    return _save(figure, path)


def radius_fidelity_figure(samples: list[RadiusFidelitySample], path: Path) -> Path:
    figure, axis = plt.subplots(figsize=(7.5, 7))
    if samples:
        source = np.asarray([sample.source_radius_um for sample in samples])
        reconstructed = np.asarray([sample.reconstructed_radius_um for sample in samples])
        errors = np.abs((reconstructed - source) / source)
        scatter = axis.scatter(source, reconstructed, c=errors, cmap="magma", s=28)
        bounds = (min(source.min(), reconstructed.min()), max(source.max(), reconstructed.max()))
        axis.plot(bounds, bounds, "--", color="#475569", label="y = x")
        figure.colorbar(scatter, ax=axis, label="absolute relative error")
        axis.text(
            0.04,
            0.96,
            f"median={np.median(errors):.3f}\nP95={np.percentile(errors, 95):.3f}\nmax={errors.max():.3f}",
            transform=axis.transAxes,
            va="top",
            bbox={"facecolor": "white", "alpha": 0.85},
        )
        axis.legend()
    else:
        axis.text(0.5, 0.5, "No valid cross-sections", transform=axis.transAxes, ha="center")
    axis.set_xlabel("source SWC radius (um)")
    axis.set_ylabel("reconstructed equivalent radius (um)")
    axis.set_title("Radius fidelity")
    axis.grid(alpha=0.25)
    return _save(figure, path)


def cross_section_figure(samples: list[RadiusFidelitySample], path: Path) -> Path:
    examples = [sample for sample in samples if sample.section_xy_um][:6]
    count = max(1, len(examples))
    columns = min(3, count)
    rows = int(np.ceil(count / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(4.5 * columns, 4.2 * rows), squeeze=False)
    for axis in axes.ravel():
        axis.set_visible(False)
    for axis, sample in zip(axes.ravel(), examples):
        axis.set_visible(True)
        section = np.asarray(sample.section_xy_um, dtype=float)
        axis.plot(section[:, 0], section[:, 1], color="#0284C7", linewidth=1.8, label="CFD section")
        angles = np.linspace(0, 2 * np.pi, 200)
        axis.plot(
            sample.source_radius_um * np.cos(angles),
            sample.source_radius_um * np.sin(angles),
            "--",
            color="#F97316",
            label="target circle",
        )
        axis.scatter([0], [0], color="#111827", s=18, label="centerline")
        axis.set_aspect("equal")
        axis.set_title(f"B{sample.branch_id}, s={sample.arc_length_um:.1f} um")
        axis.set_xlabel("u (um)")
        axis.set_ylabel("v (um)")
        axis.grid(alpha=0.2)
    if examples:
        axes.ravel()[0].legend(fontsize=8)
    else:
        axes.ravel()[0].set_visible(True)
        axes.ravel()[0].text(0.5, 0.5, "No valid cross-sections", ha="center")
    figure.suptitle("Representative lumen cross-sections")
    return _save(figure, path)


def _profile_series(
    rows: list[dict[str, object]],
    field: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    offsets = np.asarray(sorted({float(row["offset_diameters"]) for row in rows}))
    means: list[float] = []
    maxima: list[float] = []
    for offset in offsets:
        values = np.asarray(
            [
                abs(float(row[field]))
                for row in rows
                if float(row["offset_diameters"]) == offset
                and row.get(field) is not None
            ],
            dtype=float,
        )
        means.append(float(np.mean(values)) if len(values) else np.nan)
        maxima.append(float(np.max(values)) if len(values) else np.nan)
    return offsets, np.asarray(means), np.asarray(maxima)


def port_transition_profile_figure(
    v4_rows: list[dict[str, object]],
    v5_rows: list[dict[str, object]],
    path: Path,
    *,
    metric: str,
) -> Path:
    figure, axis = plt.subplots(figsize=(8.2, 5.4))
    if metric == "area":
        field = "area_relative_error"
        ylabel = "absolute area relative error"
        title = "Port transition area profile"
    else:
        field = "normal_jump_p99_deg"
        ylabel = "local normal jump P99 (deg)"
        title = "Port surface-normal profile"
    for label, rows, color in (
        ("v4 separate cylinder", v4_rows, "#DC2626"),
        ("v5 continuous centerline", v5_rows, "#0284C7"),
    ):
        offsets, means, maxima = _profile_series(rows, field)
        axis.plot(offsets, means, "o-", color=color, linewidth=2, label=f"{label}: mean")
        axis.plot(offsets, maxima, "--", color=color, alpha=0.7, label=f"{label}: max")
    axis.axvline(0.0, color="#475569", linewidth=1, linestyle=":")
    axis.set_xticks((-1.0, 0.0, 1.0, 2.0), ("cut - 1D", "cut", "cut + 1D", "cut + 2D"))
    axis.set_xlabel("position relative to active source cut")
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    return _save(figure, path)


def _local_face_ids(mesh: trimesh.Trimesh, center: np.ndarray, radius: float) -> np.ndarray:
    centers = np.asarray(mesh.triangles_center, dtype=float)
    selected = np.flatnonzero(np.linalg.norm(centers - center[None, :], axis=1) <= radius)
    return selected if len(selected) else np.arange(len(mesh.faces), dtype=np.int64)


def port_v4_v5_figure(
    v4_mesh: trimesh.Trimesh,
    v5_mesh: trimesh.Trimesh,
    port: PortGeometry,
    path: Path,
) -> Path:
    center = np.asarray(port.original_position_um, dtype=float)
    extent = 5.0 * port.radius_um
    figure = plt.figure(figsize=(12, 5.5))
    for index, (mesh, title) in enumerate(
        ((v4_mesh, "v4: branch + cylinder Boolean"), (v5_mesh, "v5: one continuous tube")),
        start=1,
    ):
        axis = figure.add_subplot(1, 2, index, projection="3d")
        faces = _local_face_ids(mesh, center, extent)
        _surface_collection(
            axis,
            mesh,
            face_indices=faces,
            alpha=0.92,
            color="#38BDF8",
            edgecolor="#334155",
            maximum_faces=8_000,
        )
        axis.scatter(*center, color="#F59E0B", marker="^", s=45)
        axis.set_xlim(center[0] - extent, center[0] + extent)
        axis.set_ylim(center[1] - extent, center[1] + extent)
        axis.set_zlim(center[2] - extent, center[2] + extent)
        axis.view_init(elev=22, azim=-58)
        axis.set_title(title)
    figure.suptitle(f"Port {port.cut_port_id}: identical camera / tessellation")
    return _save(figure, path)


def junction_transition_wireframe_figure(
    mesh: trimesh.Trimesh,
    face_region: np.ndarray,
    path: Path,
) -> Path:
    transition = np.flatnonzero(np.asarray(face_region) == 2)
    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    _surface_collection(
        axis,
        mesh,
        face_indices=transition,
        alpha=0.85,
        color="#F59E0B",
        edgecolor="#1E293B",
        maximum_faces=18_000,
    )
    points = mesh.vertices[np.unique(mesh.faces[transition])] if len(transition) else mesh.vertices
    axis.set_title("v5 TRANSITION_COLLAR boundary-loop strip wireframe")
    _equal_3d(axis, points)
    return _save(figure, path)


def junction_transition_normals_figure(
    mesh: trimesh.Trimesh,
    face_region: np.ndarray,
    path: Path,
) -> Path:
    transition = np.flatnonzero(np.asarray(face_region) == 2)
    if len(transition) > 800:
        transition = transition[np.linspace(0, len(transition) - 1, 800).astype(int)]
    figure = plt.figure(figsize=(10, 8))
    axis = figure.add_subplot(111, projection="3d")
    _surface_collection(axis, mesh, face_indices=transition, alpha=0.35, color="#A7F3D0")
    if len(transition):
        centers = np.asarray(mesh.triangles_center)[transition]
        normals = np.asarray(mesh.face_normals)[transition]
        scale = max(float(np.ptp(centers, axis=0).max()) * 0.035, 0.1)
        axis.quiver(
            centers[:, 0], centers[:, 1], centers[:, 2],
            normals[:, 0], normals[:, 1], normals[:, 2],
            length=scale, normalize=True, color="#DC2626", linewidth=0.45,
        )
        _equal_3d(axis, centers)
    axis.set_title("True transition face normals (geometry unchanged)")
    return _save(figure, path)


def junction_v4_v5_figure(
    v4_mesh: trimesh.Trimesh,
    v5_mesh: trimesh.Trimesh,
    center_um: np.ndarray,
    radius_um: float,
    junction_label: str,
    path: Path,
) -> Path:
    center = np.asarray(center_um, dtype=float)
    figure = plt.figure(figsize=(12, 5.5))
    for index, (mesh, title) in enumerate(
        ((v4_mesh, "v4 manifold overlap"), (v5_mesh, "v5 loop stitching")), start=1
    ):
        axis = figure.add_subplot(1, 2, index, projection="3d")
        faces = _local_face_ids(mesh, center, radius_um)
        _surface_collection(
            axis,
            mesh,
            face_indices=faces,
            alpha=0.9,
            color="#67E8F9",
            edgecolor="#334155",
            maximum_faces=12_000,
        )
        axis.set_xlim(center[0] - radius_um, center[0] + radius_um)
        axis.set_ylim(center[1] - radius_um, center[1] + radius_um)
        axis.set_zlim(center[2] - radius_um, center[2] + radius_um)
        axis.view_init(elev=24, azim=-60)
        axis.set_title(title)
    figure.suptitle(f"{junction_label}: identical camera / shading / tube sides")
    return _save(figure, path)


def collision_figure(
    branches: list[BranchGeometry],
    events: Iterable[CollisionEvent],
    path: Path,
) -> Path:
    events = list(events)
    figure = plt.figure(figsize=(9, 7))
    axis = figure.add_subplot(111, projection="3d")
    for branch in branches:
        axis.plot(*branch.points_um.T, color="#64748B", linewidth=1.0, alpha=0.7)
    for event in events:
        first = np.asarray(event.closest_a_um)
        second = np.asarray(event.closest_b_um)
        color = "#DC2626" if event.classification == "HARD_COLLISION" else "#F59E0B"
        axis.scatter(*first, color=color, s=45)
        axis.scatter(*second, color=color, s=45)
        axis.plot(*np.vstack((first, second)).T, color=color, linewidth=2)
    if not events:
        axis.text2D(0.03, 0.96, "No hard collision / near contact", transform=axis.transAxes)
    axis.set_title("Non-adjacent branch collision QC")
    all_points = np.concatenate([branch.points_um for branch in branches])
    _equal_3d(axis, all_points)
    return _save(figure, path)


def run_summary_figure(rows: list[dict[str, object]], path: Path) -> Path:
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    labels = [str(row["roi_id"])[-18:] for row in rows]
    triangles = [float(row.get("triangle_count", 0) or 0) for row in rows]
    p95 = [float(row.get("radius_p95_error", np.nan) or np.nan) for row in rows]
    axes[0].bar(range(len(rows)), triangles, color="#0284C7")
    axes[0].set_title("Surface triangles")
    axes[0].set_xticks(range(len(rows)), labels, rotation=55, ha="right", fontsize=7)
    axes[1].bar(range(len(rows)), p95, color="#F97316")
    axes[1].set_title("Radius fidelity P95 absolute error")
    axes[1].set_xticks(range(len(rows)), labels, rotation=55, ha="right", fontsize=7)
    figure.suptitle("CFD lumen run summary")
    return _save(figure, path)
