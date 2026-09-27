"""Human-review figures for NNE2 segmentation, graphs and direction inference."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from matplotlib.colors import Normalize, to_rgba
from matplotlib.collections import LineCollection

from ..graph.model import HierarchicalGraphResult
from .model import NNE2HierarchyResult
from .catalog import NNE2Record
from .segmentation import SegmentationResult


def _percentile_image(image: np.ndarray) -> np.ndarray:
    low, high = np.percentile(image, (1.0, 99.5))
    if high <= low:
        return np.zeros_like(image, dtype=float)
    return np.clip((image - low) / (high - low), 0.0, 1.0)


def _line_collection(
    axis: plt.Axes,
    lines: list[np.ndarray],
    dimensions: tuple[int, int],
    *,
    colors: str | list[object],
    linewidths: float | list[float],
    alpha: float = 1.0,
    linestyles: str | list[str] = "solid",
) -> None:
    projected = [line[:, dimensions] for line in lines if len(line)]
    if not projected:
        return
    axis.add_collection(
        LineCollection(
            projected,
            colors=colors,
            linewidths=linewidths,
            alpha=alpha,
            linestyles=linestyles,
            rasterized=True,
        )
    )
    all_points = np.concatenate(projected, axis=0)
    axis.set_xlim(float(all_points[:, 0].min()), float(all_points[:, 0].max()))
    axis.set_ylim(float(all_points[:, 1].min()), float(all_points[:, 1].max()))


def render_segmentation_overview(
    segmentation: SegmentationResult,
    skeleton_zyx: np.ndarray,
    output: Path,
) -> Path:
    raw_mip = np.max(segmentation.normalized_zyx, axis=0)
    score_mip = np.max(segmentation.vessel_score_zyx, axis=0)
    mask_occupancy = np.mean(segmentation.mask_zyx, axis=0)
    mask_mip = mask_occupancy > 0
    skeleton_mip = np.max(skeleton_zyx, axis=0)
    middle = segmentation.normalized_zyx.shape[0] // 2
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)
    axes[0, 0].imshow(_percentile_image(raw_mip), cmap="gray")
    axes[0, 0].set_title("Normalized volume - XY maximum projection")
    axes[0, 1].imshow(_percentile_image(score_mip), cmap="magma")
    axes[0, 1].set_title("Vessel score - XY maximum projection")
    occupancy_view = axes[0, 2].imshow(mask_occupancy, cmap="inferno")
    axes[0, 2].set_title("Cleaned mask: fraction of Z planes")
    fig.colorbar(occupancy_view, ax=axes[0, 2], fraction=0.046, pad=0.04)
    axes[1, 0].imshow(_percentile_image(segmentation.normalized_zyx[middle]), cmap="gray")
    axes[1, 0].contour(segmentation.mask_zyx[middle], levels=[0.5], colors="lime", linewidths=0.5)
    axes[1, 0].set_title(f"Middle plane z={middle + 1}: mask outline")
    axes[1, 1].imshow(_percentile_image(raw_mip), cmap="gray")
    axes[1, 1].imshow(np.ma.masked_where(~mask_mip, mask_mip), cmap="spring", alpha=0.35)
    axes[1, 1].set_title("Raw projection + mask")
    axes[1, 2].imshow(_percentile_image(raw_mip), cmap="gray")
    yy, xx = np.nonzero(skeleton_mip)
    if len(xx):
        step = max(1, len(xx) // 50_000)
        axes[1, 2].scatter(xx[::step], yy[::step], s=0.2, c="cyan", alpha=0.8)
    axes[1, 2].set_title("Raw projection + coarse centerline")
    for axis in axes.ravel():
        axis.set_axis_off()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def render_step1_component_cleanup(
    segmentation: SegmentationResult,
    output: Path,
) -> Path:
    """Show exactly what Step 1 considered, retained, and removed."""
    candidate = np.max(segmentation.candidate_mask_zyx, axis=0)
    cleaned = np.max(segmentation.mask_zyx, axis=0)
    removed = np.max(segmentation.removed_mask_zyx, axis=0)
    background = _percentile_image(np.max(segmentation.normalized_zyx, axis=0))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    axes[0].imshow(background, cmap="gray")
    axes[0].imshow(np.ma.masked_where(~candidate, candidate), cmap="autumn", alpha=0.45)
    axes[0].set_title("Step 1 candidate mask (before cleanup)")
    axes[1].imshow(background, cmap="gray")
    axes[1].imshow(np.ma.masked_where(~cleaned, cleaned), cmap="winter", alpha=0.5)
    axes[1].set_title(f"Retained components: {segmentation.component_count_after}")
    axes[2].imshow(background, cmap="gray")
    axes[2].imshow(np.ma.masked_where(~removed, removed), cmap="Reds", alpha=0.8)
    axes[2].set_title(f"Removed voxels: {segmentation.removed_voxel_count:,}")
    for axis in axes:
        axis.set_axis_off()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)
    return output


def render_step2_radius_and_centerline(
    segmentation: SegmentationResult,
    skeleton_zyx: np.ndarray,
    radius_zyx_um: np.ndarray,
    output: Path,
) -> Path:
    """Show Step 2 radius estimates and centerline in three simple views."""
    radius_mip = np.max(radius_zyx_um, axis=0)
    skeleton_xy = np.max(skeleton_zyx, axis=0)
    skeleton_xz = np.max(skeleton_zyx, axis=1)
    mask_xz = np.max(segmentation.mask_zyx, axis=1)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    view = axes[0].imshow(radius_mip, cmap="viridis")
    axes[0].set_title("Coarse vessel radius (XY maximum)")
    fig.colorbar(view, ax=axes[0], label="radius (um)", fraction=0.046, pad=0.04)
    axes[1].imshow(np.max(segmentation.mask_zyx, axis=0), cmap="gray")
    yy, xx = np.nonzero(skeleton_xy)
    axes[1].scatter(xx, yy, s=0.15, c="cyan")
    axes[1].set_title("Cleaned mask + centerline (XY)")
    axes[2].imshow(mask_xz, cmap="gray", aspect="auto")
    zz, xx = np.nonzero(skeleton_xz)
    axes[2].scatter(xx, zz, s=0.15, c="orange")
    axes[2].set_title("Cleaned mask + centerline (XZ)")
    for axis in axes:
        axis.set_axis_off()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)
    return output


def render_anchor_registration(
    hierarchy: NNE2HierarchyResult,
    records: list[NNE2Record],
    normalized_stack_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    output: Path,
    *,
    max_records: int = 6,
) -> Path:
    record_by_id = {item.record_id: item for item in records}
    anchors = hierarchy.anchors[:max_records]
    fig, axes = plt.subplots(
        len(anchors),
        2,
        figsize=(10, max(4, 4 * len(anchors))),
        squeeze=False,
        constrained_layout=True,
    )
    sx, sy, _ = spacing_xyz_um
    for row, anchor in enumerate(anchors):
        record = record_by_id[anchor.record_id]
        assert record.reference_file is not None
        with Image.open(record.reference_file) as image:
            reference = np.asarray(image.convert("RGB"))
        axes[row, 0].imshow(reference)
        axes[row, 0].set_title(
            f"Record {anchor.record_id}, BO={anchor.branching_order}: reference image"
        )
        z_index = max(0, min(normalized_stack_zyx.shape[0] - 1, anchor.matched_stack_index - 1))
        axes[row, 1].imshow(_percentile_image(normalized_stack_zyx[z_index]), cmap="gray")
        x_px = anchor.seed_xyz_um[0] / sx
        y_px = anchor.seed_xyz_um[1] / sy
        axes[row, 1].scatter(x_px, y_px, marker="+", s=150, c="red", linewidth=2)
        axes[row, 1].set_title(
            f"Matched stack frame {anchor.matched_stack_index}; score={anchor.registration_score:.3f}; "
            f"centerline distance={anchor.branch_distance_um:.1f} um"
        )
        for axis in axes[row]:
            axis.set_axis_off()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def render_tree_map(record: NNE2Record, output: Path) -> Path:
    if record.map_file is None:
        raise ValueError("A complete NNE2 record must have a subject map")
    with Image.open(record.map_file) as image:
        subject_map = np.asarray(image.convert("RGB"))
    fig, axis = plt.subplots(figsize=(8, 8), constrained_layout=True)
    axis.imshow(subject_map)
    axis.set_title(
        f"Subject {record.subject_id} surface map - inspect label Tree ID {record.tree_id}"
    )
    axis.set_axis_off()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)
    return output


def render_component_decisions(
    hierarchy: NNE2HierarchyResult,
    graph: HierarchicalGraphResult,
    output: Path,
) -> Path:
    excluded = set(hierarchy.excluded_branch_ids)
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
    retained_lines = [
        branch.points_raw_lps_um for branch in graph.branches
        if branch.branch_id not in excluded
    ]
    excluded_lines = [
        branch.points_raw_lps_um for branch in graph.branches
        if branch.branch_id in excluded
    ]
    for axis, dimensions in zip(axes, ((0, 1), (0, 2))):
        _line_collection(
            axis, retained_lines, dimensions, colors="#2676c8", linewidths=0.35, alpha=0.28
        )
        _line_collection(
            axis, excluded_lines, dimensions, colors="#d62728", linewidths=0.9, alpha=0.85
        )
    root = graph.nodes[hierarchy.root_node_id].representative_lps_um
    axes[0].scatter(root[0], root[1], marker="X", s=100, c="lime", edgecolor="black")
    axes[1].scatter(root[0], root[2], marker="X", s=100, c="lime", edgecolor="black")
    axes[0].set_title(
        f"XY: retained root component (blue), excluded islands (red={len(excluded)})"
    )
    axes[1].set_title("XZ: green X is the diving-trunk root")
    axes[0].set_xlabel("stack X (um)")
    axes[0].set_ylabel("stack Y (um)")
    axes[1].set_xlabel("stack X (um)")
    axes[1].set_ylabel("stack Z (um)")
    for axis in axes:
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)
    return output


def render_undirected_graph(
    graph: HierarchicalGraphResult,
    output: Path,
    *,
    max_branches: int = 12_000,
) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
    branches = graph.branches
    stride = max(1, int(np.ceil(len(branches) / max_branches)))
    lines = [branch.points_raw_lps_um for branch in branches[::stride]]
    _line_collection(
        axes[0], lines, (0, 1), colors="#1565c0", linewidths=0.35, alpha=0.65
    )
    _line_collection(
        axes[1], lines, (0, 2), colors="#00897b", linewidths=0.35, alpha=0.65
    )
    axes[0].set_title(f"Undirected branch graph - XY ({len(branches)} branches)")
    axes[1].set_title("Undirected branch graph - XZ")
    axes[0].set_xlabel("stack X (um)")
    axes[0].set_ylabel("stack Y (um)")
    axes[1].set_xlabel("stack X (um)")
    axes[1].set_ylabel("stack Z (um)")
    for axis in axes:
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def render_hierarchy(
    hierarchy: NNE2HierarchyResult,
    graph: HierarchicalGraphResult,
    output: Path,
    *,
    max_branches: int = 8_000,
) -> Path:
    source = {item.branch_id: item for item in graph.branches}
    shown = hierarchy.branches
    stride = max(1, int(np.ceil(len(shown) / max_branches)))
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    max_order = max(
        (item.branching_order for item in shown if item.branching_order is not None),
        default=1,
    )
    norm = Normalize(vmin=0, vmax=max(1, max_order))
    lines: list[np.ndarray] = []
    colors: list[object] = []
    styles: list[str] = []
    widths: list[float] = []
    arrow_points: list[np.ndarray] = []
    arrow_deltas: list[np.ndarray] = []
    for item in shown[::stride]:
        branch = source[item.branch_id]
        points = branch.points_raw_lps_um
        if item.upstream_node == branch.node_v:
            points = points[::-1]
        color = "#999999" if item.branching_order is None else cmap(norm(item.branching_order))
        linestyle = "--" if item.is_cross_link else "-"
        alpha = 0.35 if item.confidence in {"low", "unresolved"} else 0.8
        lines.append(points)
        colors.append(to_rgba(color, alpha))
        styles.append(linestyle)
        widths.append(0.6)
        if item.upstream_node is not None and len(points) >= 2 and len(arrow_points) < 600:
            middle = len(points) // 2
            left = max(0, middle - 1)
            right = min(len(points) - 1, middle + 1)
            delta = points[right] - points[left]
            arrow_points.append(points[middle])
            arrow_deltas.append(delta)
    for axis, dimensions in zip(axes, ((0, 1), (0, 2))):
        _line_collection(
            axis,
            lines,
            dimensions,
            colors=colors,
            linewidths=widths,
            linestyles=styles,
        )
    if arrow_points:
        arrow_points_array = np.asarray(arrow_points)
        arrow_deltas_array = np.asarray(arrow_deltas)
        axes[0].quiver(
            arrow_points_array[:, 0], arrow_points_array[:, 1],
            arrow_deltas_array[:, 0], arrow_deltas_array[:, 1],
            color="#202020", angles="xy", scale_units="xy", scale=1, width=0.0014, alpha=0.55,
        )
        axes[1].quiver(
            arrow_points_array[:, 0], arrow_points_array[:, 2],
            arrow_deltas_array[:, 0], arrow_deltas_array[:, 2],
            color="#202020", angles="xy", scale_units="xy", scale=1, width=0.0014, alpha=0.55,
        )
    for anchor in hierarchy.anchors:
        x, y, z = anchor.seed_xyz_um
        marker = "*" if anchor.branching_order == 0 else "o"
        axes[0].scatter(x, y, marker=marker, s=70, edgecolor="black", linewidth=0.5, label=None, c=[cmap(norm(anchor.branching_order))])
        axes[1].scatter(x, z, marker=marker, s=70, edgecolor="black", linewidth=0.5, c=[cmap(norm(anchor.branching_order))])
    root_node = graph.nodes[hierarchy.root_node_id]
    rx, ry, rz = root_node.representative_lps_um
    axes[0].scatter(rx, ry, marker="X", s=110, c="red", edgecolor="white", linewidth=0.8)
    axes[1].scatter(rx, rz, marker="X", s=110, c="red", edgecolor="white", linewidth=0.8)
    axes[0].set_title(f"{hierarchy.tree_key}: parent -> child, XY")
    axes[1].set_title("Parent -> child, XZ (red X = root)")
    axes[0].set_xlabel("stack X (um)")
    axes[0].set_ylabel("stack Y (um)")
    axes[1].set_xlabel("stack X (um)")
    axes[1].set_ylabel("stack Z (um)")
    for axis in axes:
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.15)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(scalar, ax=axes, label="Branching Order")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)
    return output


def render_step3_diagnostics(
    graph: HierarchicalGraphResult,
    output: Path,
    *,
    short_branch_warning_um: float,
    large_junction_warning_voxels: int,
    high_degree_warning: int,
    max_branches: int = 6_000,
) -> Path:
    """Fast, compact Step 3 review board for large NNE2 branch graphs."""
    stride = max(1, int(np.ceil(len(graph.branches) / max_branches)))
    selected = graph.branches[::stride]
    lines = [branch.points_raw_lps_um for branch in selected]
    colors = ["#d95f02" if branch.cycle_ids else "#2b6ca3" for branch in selected]
    lengths = np.asarray([branch.arc_length_raw_um[-1] for branch in graph.branches])
    radii = np.asarray([np.mean(branch.coarse_radius_raw_um) for branch in graph.branches])
    degrees = np.asarray([node.graph_degree for node in graph.nodes])
    short_count = int(np.count_nonzero(lengths < short_branch_warning_um))
    large_count = sum(
        len(node.voxel_indices_xyz) >= large_junction_warning_voxels for node in graph.nodes
    )
    high_degree_count = int(np.count_nonzero(degrees >= high_degree_warning))
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
    for axis, dimensions, title in zip(
        axes[0], ((0, 1), (0, 2), (1, 2)), ("XY view", "XZ view", "YZ view")
    ):
        _line_collection(axis, lines, dimensions, colors=colors, linewidths=0.42, alpha=0.65)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(f"{title}: orange branches belong to stored cycles")
        axis.set_axis_off()
    axes[1, 0].hist(lengths, bins=45, color="#4c78a8")
    axes[1, 0].axvline(short_branch_warning_um, color="red", linestyle="--")
    axes[1, 0].set_title(f"Branch length; short-warning count={short_count:,}")
    axes[1, 0].set_xlabel("length (um)")
    axes[1, 1].hist(radii, bins=45, color="#72a66a")
    axes[1, 1].set_title("Coarse mean radius (navigation estimate)")
    axes[1, 1].set_xlabel("radius (um)")
    axes[1, 2].axis("off")
    axes[1, 2].text(
        0.02,
        0.95,
        "\n".join(
            (
                f"Nodes: {len(graph.nodes):,}",
                f"Branches: {len(graph.branches):,}",
                f"Skeleton components: {graph.skeleton_component_count:,}",
                f"Independent cycles: {graph.cycle_rank:,}",
                f"Short branches: {short_count:,}",
                f"Large junction regions: {large_count:,}",
                f"High-degree nodes: {high_degree_count:,}",
                f"Round-trip missing/extra: {graph.missing_voxel_count}/{graph.extra_voxel_count}",
            )
        ),
        va="top",
        fontsize=12,
    )
    axes[1, 2].set_title("Step 3 acceptance summary")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)
    return output
