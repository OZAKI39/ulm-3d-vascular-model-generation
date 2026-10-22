"""Visual acceptance artifacts for the hierarchical vascular graph."""

from __future__ import annotations

from itertools import cycle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from ..config import HierarchicalGraphConfig
from .model import HierarchicalGraphResult


NODE_COLORS = {
    "terminal": "#2f72b7",
    "junction": "#d94841",
    "complex_junction": "#8f3faf",
    "connector": "#666666",
    "cycle_anchor": "#e69f24",
    "isolated": "#111111",
}


def _set_equal_3d(axis: plt.Axes, points: np.ndarray) -> None:
    if not len(points):
        return
    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    center = (minimum + maximum) / 2
    radius = max(float(np.max(maximum - minimum)) / 2, 1.0)
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1))


def _all_raw_points(result: HierarchicalGraphResult) -> np.ndarray:
    return np.concatenate([item.points_raw_lps_um for item in result.branches], axis=0)


def _project(points: np.ndarray, view: int) -> tuple[np.ndarray, np.ndarray]:
    axes = ((1, 2), (0, 2), (0, 1))[view]
    return points[:, axes[0]], points[:, axes[1]]


def _render_graph_3d(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "hierarchical_graph_3d.png"
    figure = plt.figure(figsize=(13, 10), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    regular = [
        item.points_smoothed_lps_um for item in result.branches if not item.cycle_ids
    ]
    cyclic = [
        item.points_smoothed_lps_um for item in result.branches if item.cycle_ids
    ]
    if regular:
        axis.add_collection3d(
            Line3DCollection(regular, colors="#557f98", linewidths=0.55, alpha=0.7),
            autolim=False,
        )
    if cyclic:
        axis.add_collection3d(
            Line3DCollection(cyclic, colors="#e69f24", linewidths=1.2, alpha=0.95),
            autolim=False,
        )
    for node_type, color in NODE_COLORS.items():
        coordinates = np.asarray(
            [item.representative_lps_um for item in result.nodes if item.node_type == node_type]
        )
        if len(coordinates):
            axis.scatter(
                coordinates[:, 0],
                coordinates[:, 1],
                coordinates[:, 2],
                s=10 if node_type == "terminal" else 18,
                c=color,
                label=node_type.replace("_", " "),
                depthshade=False,
            )
    _set_equal_3d(axis, _all_raw_points(result))
    axis.set_xlabel("LPS X (um)")
    axis.set_ylabel("LPS Y (um)")
    axis.set_zlabel("LPS Z (um)")
    axis.set_title("Whole-network hierarchical vascular graph\nOrange branches belong to an independent cycle basis")
    axis.legend(loc="upper right", fontsize=8)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_three_views(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "hierarchical_graph_three_views.png"
    figure, axes = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    for view, axis in enumerate(axes):
        for branch in result.branches:
            x, y = _project(branch.points_smoothed_lps_um, view)
            axis.plot(x, y, color="#557f98", linewidth=0.45, alpha=0.7)
        for node_type, color in NODE_COLORS.items():
            coordinates = np.asarray(
                [item.representative_lps_um for item in result.nodes if item.node_type == node_type]
            )
            if len(coordinates):
                x, y = _project(coordinates, view)
                axis.scatter(x, y, s=5 if node_type == "terminal" else 10, c=color)
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(("X view", "Y view", "Z view")[view])
    figure.suptitle("Junction/terminal nodes and full branch geometry")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_skeleton_overlay(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "graph_skeleton_overlay.png"
    node_mask = np.zeros_like(result.source_skeleton, dtype=bool)
    for node in result.nodes:
        node_mask[tuple(node.voxel_indices_xyz.T)] = True
    figure, axes = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    for view, axis in enumerate(axes):
        skeleton_projection = result.source_skeleton.max(axis=view).T
        reconstructed_projection = result.reconstructed_skeleton.max(axis=view).T
        node_projection = node_mask.max(axis=view).T
        axis.imshow(skeleton_projection, cmap="gray", origin="lower", interpolation="nearest")
        axis.imshow(
            np.ma.masked_where(~reconstructed_projection, reconstructed_projection),
            cmap="Greens",
            origin="lower",
            interpolation="nearest",
            alpha=0.45,
        )
        axis.imshow(
            np.ma.masked_where(~node_projection, node_projection),
            cmap="autumn",
            origin="lower",
            interpolation="nearest",
            alpha=0.9,
        )
        axis.axis("off")
        axis.set_title(("X view", "Y view", "Z view")[view])
    figure.suptitle("Raw skeleton (white), graph reconstruction (green), node regions (orange)")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_difference(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "graph_reconstruction_difference.png"
    missing = result.source_skeleton & ~result.reconstructed_skeleton
    extra = result.reconstructed_skeleton & ~result.source_skeleton
    figure, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    for view in range(3):
        axes[0, view].imshow(missing.max(axis=view).T, cmap="Reds", origin="lower")
        axes[0, view].set_title(f"Missing - {( 'X', 'Y', 'Z')[view]} view")
        axes[1, view].imshow(extra.max(axis=view).T, cmap="Purples", origin="lower")
        axes[1, view].set_title(f"Extra - {( 'X', 'Y', 'Z')[view]} view")
        axes[0, view].axis("off")
        axes[1, view].axis("off")
    figure.suptitle(
        f"Round-trip difference: missing={result.missing_voxel_count}, extra={result.extra_voxel_count}"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_branch_as_node(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "branch_as_node_graph.png"
    midpoints = {
        item.branch_id: item.points_smoothed_lps_um[len(item.points_smoothed_lps_um) // 2]
        for item in result.branches
    }
    figure, axes = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    for view, axis in enumerate(axes):
        relation_segments = []
        for left, right in result.branch_as_node_graph.edges():
            relation_segments.append(np.vstack((midpoints[int(left)], midpoints[int(right)])))
        for segment in relation_segments:
            x, y = _project(segment, view)
            axis.plot(x, y, color="#b6bdc4", linewidth=0.35, alpha=0.45)
        points = np.asarray([midpoints[index] for index in sorted(midpoints)])
        x, y = _project(points, view)
        axis.scatter(x, y, s=4, c="#386f9e", alpha=0.8)
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(("X view", "Y view", "Z view")[view])
    figure.suptitle("Branch-as-node graph: one blue point represents one complete branch")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_cycles(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "cycle_overview.png"
    palette = cycle(plt.cm.tab20(np.linspace(0, 1, 20)))
    colors = {record.cycle_id: next(palette) for record in result.cycles}
    figure, axes = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    for view, axis in enumerate(axes):
        for branch in result.branches:
            x, y = _project(branch.points_smoothed_lps_um, view)
            if branch.cycle_ids:
                color = colors[branch.cycle_ids[0]]
                axis.plot(x, y, color=color, linewidth=1.0, alpha=0.95)
            else:
                axis.plot(x, y, color="#c8cdd1", linewidth=0.3, alpha=0.35)
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(("X view", "Y view", "Z view")[view])
    figure.suptitle(
        f"Cycle-preserving aggregation: independent cycle rank = {result.cycle_rank}"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_distributions(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "branch_morphology_distributions.png"
    lengths = np.asarray([item.arc_length_raw_um[-1] for item in result.branches])
    tortuosities = np.asarray([item.tortuosity_smoothed for item in result.branches])
    radii = np.asarray([np.mean(item.coarse_radius_raw_um) for item in result.branches])
    degrees = np.asarray([item.graph_degree for item in result.nodes])
    figure, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    axes[0, 0].hist(lengths, bins=40, color="#4c7fa5")
    axes[0, 0].set_title("Raw branch length")
    axes[0, 0].set_xlabel("um")
    finite_tortuosity = tortuosities[np.isfinite(tortuosities)]
    axes[0, 1].hist(finite_tortuosity, bins=40, color="#6d9f71")
    axes[0, 1].set_title("Smoothed tortuosity")
    axes[0, 1].set_xlabel("length / straight distance")
    axes[1, 0].hist(radii, bins=40, color="#b080c0")
    axes[1, 0].set_title("Coarse mean radius (navigation only)")
    axes[1, 0].set_xlabel("um")
    unique, counts = np.unique(degrees, return_counts=True)
    axes[1, 1].bar(unique, counts, color="#d77855")
    axes[1, 1].set_title("Graph-node degree")
    axes[1, 1].set_xlabel("incident branch endpoints")
    axes[1, 1].set_ylabel("node count")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_radius_profiles(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "coarse_radius_profiles.png"
    selected = sorted(
        result.branches, key=lambda item: item.arc_length_raw_um[-1], reverse=True
    )[:12]
    figure, axis = plt.subplots(figsize=(12, 7), constrained_layout=True)
    for branch in selected:
        total = max(float(branch.arc_length_raw_um[-1]), np.finfo(float).eps)
        axis.plot(
            branch.arc_length_raw_um / total,
            branch.coarse_radius_raw_um,
            linewidth=1.0,
            label=f"branch {branch.branch_id}",
        )
    axis.set_xlabel("Normalized position along branch")
    axis.set_ylabel("Coarse radius (um)")
    axis.set_title("Coarse radius profiles of the 12 longest branches\nNavigation estimates, not final radius truth")
    axis.legend(ncol=3, fontsize=7)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_raw_vs_smoothed(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "raw_vs_smoothed_centerlines.png"
    selected = sorted(
        result.branches, key=lambda item: item.arc_length_raw_um[-1], reverse=True
    )[:6]
    figure = plt.figure(figsize=(14, 10))
    figure.subplots_adjust(top=0.86, bottom=0.03, wspace=0.04, hspace=0.18)
    for plot_index, branch in enumerate(selected, start=1):
        axis = figure.add_subplot(2, 3, plot_index, projection="3d")
        raw = branch.points_raw_lps_um
        smooth = branch.points_smoothed_lps_um
        axis.plot(raw[:, 0], raw[:, 1], raw[:, 2], color="#888888", linewidth=1, label="raw")
        axis.plot(
            smooth[:, 0], smooth[:, 1], smooth[:, 2], color="#d94b45", linewidth=1.4, label="smoothed"
        )
        _set_equal_3d(axis, raw)
        axis.set_title(
            f"Branch {branch.branch_id}; max shift {branch.smoothing_max_deviation_um:.2f} um",
            fontsize=9,
        )
        axis.set_axis_off()
    figure.legend(
        ["raw", "smoothed"],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=2,
    )
    figure.suptitle(
        "Raw geometry is retained; smoothing is a separate derived sequence",
        y=0.985,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_warning_locations(
    result: HierarchicalGraphResult,
    config: HierarchicalGraphConfig,
    output_dir: Path,
) -> Path:
    """Show exactly where threshold-based manual-review warnings occur."""
    path = output_dir / "graph_warning_locations.png"
    short_branches = [
        branch
        for branch in result.branches
        if branch.arc_length_raw_um[-1] < config.short_branch_warning_um
    ]
    high_degree_nodes = [
        node
        for node in result.nodes
        if node.graph_degree >= config.high_degree_warning
    ]
    large_nodes = [
        node
        for node in result.nodes
        if len(node.voxel_indices_xyz) >= config.large_junction_warning_voxels
    ]

    figure, axes = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    for view, axis in enumerate(axes):
        for branch in result.branches:
            x, y = _project(branch.points_smoothed_lps_um, view)
            axis.plot(x, y, color="#c8cdd1", linewidth=0.35, alpha=0.45)
        for branch in short_branches:
            x, y = _project(branch.points_smoothed_lps_um, view)
            axis.plot(x, y, color="#e69f24", linewidth=2.0, alpha=1.0)
        if large_nodes:
            coordinates = np.asarray(
                [node.representative_lps_um for node in large_nodes], dtype=float
            )
            x, y = _project(coordinates, view)
            axis.scatter(x, y, s=55, c="#d94841", marker="s", zorder=4)
        if high_degree_nodes:
            coordinates = np.asarray(
                [node.representative_lps_um for node in high_degree_nodes], dtype=float
            )
            x, y = _project(coordinates, view)
            axis.scatter(x, y, s=75, c="#8f3faf", marker="*", zorder=5)
            for node, x_value, y_value in zip(
                high_degree_nodes, x, y, strict=True
            ):
                axis.annotate(
                    f"node {node.node_id}\ndegree {node.graph_degree}",
                    (x_value, y_value),
                    xytext=(6, 6),
                    textcoords="offset points",
                    fontsize=7,
                    color="#5b216b",
                )
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(("X view", "Y view", "Z view")[view])

    handles = [
        Line2D([0], [0], color="#c8cdd1", lw=1, label="other branches"),
        Line2D(
            [0],
            [0],
            color="#e69f24",
            lw=3,
            label=f"short branches < {config.short_branch_warning_um:g} um",
        ),
        Line2D(
            [0],
            [0],
            marker="*",
            color="none",
            markerfacecolor="#8f3faf",
            markeredgecolor="#8f3faf",
            markersize=11,
            label=f"node degree >= {config.high_degree_warning}",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor="#d94841",
            markeredgecolor="#d94841",
            markersize=7,
            label=(
                "junction region >= "
                f"{config.large_junction_warning_voxels} voxels"
            ),
        ),
    ]
    figure.legend(handles=handles, loc="lower center", ncol=4, fontsize=8)
    figure.suptitle(
        "Manual-review locations: "
        f"short branches={len(short_branches)}, "
        f"high-degree nodes={len(high_degree_nodes)}, "
        f"large junction regions={len(large_nodes)}"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_junction_geometry(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "junction_geometry.png"
    selected = sorted(
        result.junctions,
        key=lambda item: len(item.incident_directions),
        reverse=True,
    )[:6]
    nodes = {item.node_id: item for item in result.nodes}
    figure = plt.figure(figsize=(14, 10), constrained_layout=True)
    for plot_index, junction in enumerate(selected, start=1):
        axis = figure.add_subplot(2, 3, plot_index, projection="3d")
        origin = np.asarray(nodes[junction.node_id].representative_lps_um)
        for incident in junction.incident_directions:
            direction = np.asarray(incident["outward_direction_lps"])
            axis.quiver(*origin, *direction, length=12.0, normalize=True)
        axis.scatter(*origin, c="#d94841", s=35)
        axis.set_title(
            f"Node {junction.node_id}: degree {nodes[junction.node_id].graph_degree}", fontsize=9
        )
        axis.set_xlim(origin[0] - 15, origin[0] + 15)
        axis.set_ylim(origin[1] - 15, origin[1] + 15)
        axis.set_zlim(origin[2] - 15, origin[2] + 15)
        axis.set_axis_off()
    figure.suptitle("Outward branch directions at selected junctions (flow direction unknown)")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_fidelity(result: HierarchicalGraphResult, output_dir: Path) -> Path:
    path = output_dir / "representation_fidelity_dashboard.png"
    labels = ["Source voxels", "Represented voxels", "Missing", "Extra", "Duplicated interior"]
    values = [
        result.skeleton_voxel_count,
        result.represented_voxel_count,
        result.missing_voxel_count,
        result.extra_voxel_count,
        result.duplicate_interior_voxel_count,
    ]
    colors = ["#4c7fa5", "#2d936c", "#d64545", "#8f3faf", "#e69f24"]
    figure, axis = plt.subplots(figsize=(11, 6), constrained_layout=True)
    bars = axis.bar(labels, values, color=colors)
    axis.bar_label(bars, fmt="%d")
    axis.set_ylabel("Voxel count")
    axis.set_title(
        f"Skeleton-to-graph round-trip fidelity\ncomponents: skeleton={result.skeleton_component_count}, graph={result.graph_component_count}; cycle rank={result.cycle_rank}"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def render_hierarchical_graph_visualizations(
    result: HierarchicalGraphResult,
    config: HierarchicalGraphConfig,
    output_dir: Path,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    return [
        _render_graph_3d(result, output_dir),
        _render_three_views(result, output_dir),
        _render_skeleton_overlay(result, output_dir),
        _render_difference(result, output_dir),
        _render_branch_as_node(result, output_dir),
        _render_cycles(result, output_dir),
        _render_distributions(result, output_dir),
        _render_radius_profiles(result, output_dir),
        _render_raw_vs_smoothed(result, output_dir),
        _render_junction_geometry(result, output_dir),
        _render_fidelity(result, output_dir),
        _render_warning_locations(result, config, output_dir),
    ]
