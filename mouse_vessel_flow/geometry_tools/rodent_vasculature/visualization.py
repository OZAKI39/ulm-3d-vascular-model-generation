"""Acceptance visualizations with explicit parent-to-current arrows."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .model import DirectedVascularGraph


def sample_direction_arrows(
    result: DirectedVascularGraph, maximum: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    candidates: list[tuple[np.ndarray, np.ndarray, int]] = []
    for branch in result.branches:
        points, directions = branch.derived_points_um, branch.direction_vectors_xyz
        if len(points) < 2:
            continue
        local_indices = np.linspace(0, len(points) - 1, min(4, len(points)), dtype=int)
        for index in np.unique(local_indices):
            if np.linalg.norm(directions[index]) > 0:
                candidates.append((points[index], directions[index], branch.depth))
    if len(candidates) > maximum:
        keep = np.linspace(0, len(candidates) - 1, maximum, dtype=int)
        candidates = [candidates[index] for index in keep]
    if not candidates:
        return np.empty((0, 3)), np.empty((0, 3)), np.empty(0, dtype=int)
    return (
        np.asarray([item[0] for item in candidates]),
        np.asarray([item[1] for item in candidates]),
        np.asarray([item[2] for item in candidates]),
    )


def _draw_branches_3d(axis: plt.Axes, result: DirectedVascularGraph) -> None:
    maximum_depth = max((branch.depth for branch in result.branches), default=1)
    colors = plt.get_cmap("viridis")
    for branch in result.branches:
        points = branch.derived_points_um
        if len(points):
            axis.plot(
                points[:, 0], points[:, 1], points[:, 2],
                color=colors(branch.depth / max(maximum_depth, 1)), linewidth=0.7, alpha=0.85,
            )


def save_direction_3d(
    result: DirectedVascularGraph, path: Path, *, max_arrows: int
) -> Path:
    points, directions, depths = sample_direction_arrows(result, max_arrows)
    figure = plt.figure(figsize=(10, 8), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    _draw_branches_3d(axis, result)
    if len(points):
        typical = np.median(
            [branch.derived_length_um for branch in result.branches if branch.derived_length_um > 0]
            or [10.0]
        )
        axis.quiver(
            points[:, 0], points[:, 1], points[:, 2],
            directions[:, 0], directions[:, 1], directions[:, 2],
            length=max(2.0, min(12.0, float(typical) * 0.12)), normalize=True,
            colors=plt.get_cmap("plasma")(depths / max(int(depths.max()), 1)),
            linewidth=0.8, arrow_length_ratio=0.45,
        )
    axis.set_xlabel("X (µm)")
    axis.set_ylabel("Y (µm)")
    axis.set_zlabel("Z (µm)")
    axis.set_title("SWC structural direction: parent → current node\n(arrows are not measured flow velocity)")
    try:
        axis.set_box_aspect(np.ptp(result.swc.points_um, axis=0).clip(min=1.0))
    except AttributeError:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _backgrounds(
    image_volume: np.ndarray | None,
    mask_volume: np.ndarray | None,
    spacing_xyz: tuple[float, float, float],
) -> dict[str, tuple[np.ndarray, list[float]]] | None:
    volume = mask_volume if mask_volume is not None else image_volume
    if volume is None:
        return None
    z, y, x = volume.shape
    sx, sy, sz = spacing_xyz
    return {
        "XY": (np.max(volume, axis=0), [0, x * sx, y * sy, 0]),
        "XZ": (np.max(volume, axis=1), [0, x * sx, z * sz, 0]),
        "YZ": (np.max(volume, axis=2), [0, y * sy, z * sz, 0]),
    }


def save_direction_orthogonal(
    result: DirectedVascularGraph,
    path: Path,
    *,
    max_arrows: int,
    spacing_xyz_um: tuple[float, float, float],
    image_volume: np.ndarray | None = None,
    mask_volume: np.ndarray | None = None,
) -> Path:
    arrow_points, arrow_directions, depths = sample_direction_arrows(result, max_arrows)
    backgrounds = _backgrounds(image_volume, mask_volume, spacing_xyz_um)
    projections = (("XY", 0, 1), ("XZ", 0, 2), ("YZ", 1, 2))
    figure, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    maximum_depth = max((branch.depth for branch in result.branches), default=1)
    for axis, (name, horizontal, vertical) in zip(axes, projections):
        if backgrounds:
            plane, extent = backgrounds[name]
            axis.imshow(plane, cmap="gray", extent=extent, alpha=0.45, interpolation="nearest")
        for branch in result.branches:
            points = branch.derived_points_um
            axis.plot(
                points[:, horizontal], points[:, vertical], linewidth=0.65,
                color=plt.get_cmap("viridis")(branch.depth / max(maximum_depth, 1)), alpha=0.8,
            )
        if len(arrow_points):
            axis.quiver(
                arrow_points[:, horizontal], arrow_points[:, vertical],
                arrow_directions[:, horizontal], arrow_directions[:, vertical],
                depths, cmap="plasma", angles="xy", scale_units="xy", scale=0.12,
                width=0.003, headwidth=4.5, headlength=5.5, pivot="tail",
            )
        axis.set_title(f"{name} projection: parent → current")
        axis.set_xlabel(f"{'XYZ'[horizontal]} (µm)")
        axis.set_ylabel(f"{'XYZ'[vertical]} (µm)")
        axis.set_aspect("equal", adjustable="box")
    figure.suptitle("Arrow direction is inferred from SWC parent topology, not pressure/velocity", fontsize=13)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def save_branch_hierarchy(result: DirectedVascularGraph, path: Path) -> Path:
    figure, axis = plt.subplots(figsize=(11, 8), constrained_layout=True)
    positions = {
        int(node): (float(data["x_um"]), float(data["y_um"]))
        for node, data in result.junction_graph.nodes(data=True)
    }
    for branch in result.branches:
        points = branch.derived_points_um
        axis.plot(points[:, 0], points[:, 1], color="0.7", linewidth=0.6, zorder=1)
        if len(points) >= 2:
            middle = max(0, len(points) // 2 - 1)
            start, end = points[middle], points[middle + 1]
            axis.annotate(
                "", xy=(end[0], end[1]), xytext=(start[0], start[1]),
                arrowprops={"arrowstyle": "-|>", "color": "tab:blue", "lw": 0.8}, zorder=2,
            )
    role_style = {
        "inferred_inlet": ("tab:green", "^"),
        "inferred_outlet": ("tab:red", "v"),
        "divergence_junction": ("tab:orange", "o"),
        "convergence_junction": ("tab:purple", "s"),
    }
    for role, (color, marker) in role_style.items():
        nodes = [node for node, data in result.junction_graph.nodes(data=True) if data["role"] == role]
        if nodes:
            axis.scatter(
                [positions[node][0] for node in nodes], [positions[node][1] for node in nodes],
                s=22, c=color, marker=marker, label=role.replace("_", " "), zorder=3,
            )
    axis.set_title("Directed branch topology (XY): arrows follow SWC parent → current")
    axis.set_xlabel("X (µm)")
    axis.set_ylabel("Y (µm)")
    axis.set_aspect("equal", adjustable="box")
    axis.legend(loc="best", fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def create_visualizations(
    result: DirectedVascularGraph,
    output_dir: Path,
    *,
    max_arrows: int,
    spacing_xyz_um: tuple[float, float, float],
    image_volume: np.ndarray | None = None,
    mask_volume: np.ndarray | None = None,
) -> list[Path]:
    return [
        save_direction_3d(result, output_dir / "direction_parent_to_current_3d.png", max_arrows=max_arrows),
        save_direction_orthogonal(
            result,
            output_dir / "direction_parent_to_current_orthogonal.png",
            max_arrows=max_arrows,
            spacing_xyz_um=spacing_xyz_um,
            image_volume=image_volume,
            mask_volume=mask_volume,
        ),
        save_branch_hierarchy(result, output_dir / "directed_branch_topology_xy.png"),
    ]
