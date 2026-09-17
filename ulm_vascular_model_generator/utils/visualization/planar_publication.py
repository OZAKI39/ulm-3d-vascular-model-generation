"""Flat, quantitative publication figure for an X-Z planar SWC tree."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, PowerNorm
from matplotlib.lines import Line2D
import numpy as np


def is_xz_planar(tree, tolerance=1e-6):
    """Return True when all nodes lie in one constant-Y plane."""
    coordinates = tree.coordinates
    if len(coordinates) == 0:
        return False
    reference_span = max(
        float(np.ptp(coordinates[:, 0])),
        float(np.ptp(coordinates[:, 2])),
        1.0,
    )
    return float(np.ptp(coordinates[:, 1])) <= tolerance * reference_span


def _tree_topology(tree):
    children = {node_id: [] for node_id in tree.nodes}
    roots = []
    for node_id, node in tree.nodes.items():
        if node.parent_id < 0:
            roots.append(node_id)
        else:
            children[node.parent_id].append(node_id)

    root_id = roots[0]
    generation = {root_id: 0}
    path_distance = {root_id: 0.0}
    stack = [root_id]
    while stack:
        parent_id = stack.pop()
        parent = tree.nodes[parent_id]
        for child_id in children[parent_id]:
            child = tree.nodes[child_id]
            generation[child_id] = generation[parent_id] + 1
            path_distance[child_id] = (
                path_distance[parent_id]
                + float(np.linalg.norm(child.xyz - parent.xyz))
            )
            stack.append(child_id)

    bifurcations = [
        node_id for node_id, child_ids in children.items() if len(child_ids) > 1
    ]
    terminals = [
        node_id for node_id, child_ids in children.items() if not child_ids
    ]
    return root_id, children, bifurcations, terminals, generation, path_distance


def _style(base_font):
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": base_font,
            "axes.labelsize": base_font,
            "axes.titlesize": base_font + 1.0,
            "xtick.labelsize": base_font - 1.5,
            "ytick.labelsize": base_font - 1.5,
            "legend.fontsize": base_font - 2.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def _relative_geometry(tree):
    x_min = min(node.xyz[0] for node in tree.nodes.values())
    z_min = min(node.xyz[2] for node in tree.nodes.values())
    points = {
        node_id: np.array([node.xyz[0] - x_min, node.xyz[2] - z_min])
        for node_id, node in tree.nodes.items()
    }
    segments = np.asarray(
        [[points[parent_id], points[child_id]] for parent_id, child_id in tree.edges]
    )
    radii = np.asarray(
        [tree.nodes[child_id].radius_um for _, child_id in tree.edges],
        dtype=float,
    )
    return points, segments, radii, float(x_min), float(z_min)


def _add_flow_arrows(ax, tree, points, generation):
    """Add a few uncluttered parent-to-child arrows across the tree."""
    maximum_generation = max(generation.values())
    target_generations = sorted(
        {1, max(2, maximum_generation // 3), max(3, 2 * maximum_generation // 3)}
    )
    for target in target_generations:
        candidates = [
            edge for edge in tree.edges if generation[edge[1]] == target
        ]
        if not candidates:
            continue
        parent_id, child_id = max(
            candidates,
            key=lambda edge: tree.nodes[edge[1]].radius_um,
        )
        start = points[parent_id]
        end = points[child_id]
        arrow_start = start + 0.42 * (end - start)
        arrow_end = start + 0.70 * (end - start)
        ax.annotate(
            "",
            xy=arrow_end,
            xytext=arrow_start,
            arrowprops={
                "arrowstyle": "-|>",
                "color": "#20252a",
                "lw": 0.9,
                "mutation_scale": 8,
            },
            zorder=7,
        )


def _add_scale_bar(ax):
    x_lower, x_upper = ax.get_xlim()
    z_lower, z_upper = ax.get_ylim()
    x_span = x_upper - x_lower
    z_span = z_upper - z_lower
    span = max(x_span, z_span, 1.0)
    target = 0.22 * span
    magnitude = 10.0 ** np.floor(np.log10(target))
    length = max(
        value for value in (1.0, 2.0, 5.0, 10.0) if value * magnitude <= target
    ) * magnitude
    x0 = x_lower + 0.08 * x_span
    z0 = z_lower + 0.08 * z_span
    ax.plot([x0, x0 + length], [z0, z0], color="#20252a", lw=2.2)
    ax.text(
        x0 + 0.5 * length,
        z0 + 0.04 * z_span,
        f"{length:g} \u00b5m",
        ha="center",
        va="bottom",
    )


def _representative_bifurcation(bifurcations, path_distance):
    distances = np.asarray([path_distance[node_id] for node_id in bifurcations])
    middle = float(np.median(distances))
    return min(
        bifurcations,
        key=lambda node_id: abs(path_distance[node_id] - middle),
    )


def _bifurcation_opening_angles(tree, children):
    """Return the largest daughter-to-daughter opening angle at each branch."""
    result = {}
    for node_id, child_ids in children.items():
        if len(child_ids) < 2:
            continue
        vectors = [
            tree.nodes[child_id].xyz - tree.nodes[node_id].xyz
            for child_id in child_ids
        ]
        pair_angles = []
        for first_index, first in enumerate(vectors):
            for second in vectors[first_index + 1 :]:
                denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
                if denominator <= 1.0e-12:
                    continue
                cosine = float(np.dot(first, second)) / denominator
                pair_angles.append(
                    float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
                )
        if pair_angles:
            result[node_id] = max(pair_angles)
    return result


def _add_opening_angle_annotation(ax, center, daughter_points, angle_deg):
    """Draw a compact angle arc inside a representative bifurcation."""
    if len(daughter_points) < 2:
        return
    vectors = [np.asarray(point) - center for point in daughter_points[:2]]
    lengths = [float(np.linalg.norm(vector)) for vector in vectors]
    if min(lengths) <= 1.0e-12:
        return
    first_angle = float(np.arctan2(vectors[0][1], vectors[0][0]))
    second_angle = float(np.arctan2(vectors[1][1], vectors[1][0]))
    signed_sweep = (second_angle - first_angle + np.pi) % (2.0 * np.pi) - np.pi
    arc_angles = first_angle + np.linspace(0.0, signed_sweep, 80)
    radius = 0.24 * min(lengths)
    arc = center + radius * np.column_stack(
        (np.cos(arc_angles), np.sin(arc_angles))
    )
    ax.plot(arc[:, 0], arc[:, 1], color="#d1495b", linewidth=1.6, zorder=7)
    middle_angle = first_angle + 0.5 * signed_sweep
    label_position = center + 1.42 * radius * np.asarray(
        [np.cos(middle_angle), np.sin(middle_angle)]
    )
    ax.text(
        label_position[0],
        label_position[1],
        f"{angle_deg:.1f}°",
        ha="center",
        va="center",
        color="#a62f43",
        fontweight="bold",
        zorder=8,
    )


def _write_caption(
    path,
    tree,
    y_value,
    x_min,
    z_min,
    bifurcation_count,
    terminal_count,
    median_opening_angle,
):
    caption = (
        "Planar vascular-tree morphology and radius distribution. "
        f"(a) The network is a genuine two-dimensional X-Z structure "
        f"(all nodes have Y = {y_value:.1f} \u00b5m), not a projection of a "
        "three-dimensional network. The orthographic data-space view preserves "
        "lengths and angles. Coordinates are shown relative to "
        f"(Xmin, Zmin) = ({x_min:.1f}, {z_min:.1f}) \u00b5m. Segment line width "
        "is constant and only color represents local radius. Filled circles "
        "mark true bifurcations, open circles mark terminal nodes, the diamond "
        "marks the root, and arrows indicate the parent-to-child flow direction. "
        "(b) Representative bifurcation with its daughter opening angle. "
        "(c) Segment-radius distribution. (d) Daughter opening-angle "
        f"distribution (median {median_opening_angle:.1f} degrees). "
        f"The tree contains {tree.edge_count} segments, {bifurcation_count} "
        f"bifurcations, and {terminal_count} terminal nodes."
    )
    path.write_text(caption + "\n", encoding="utf-8")


def render_planar_publication_figure(
    tree,
    output_path,
    color_map="viridis",
    font_size=32,
):
    """
    Save a four-panel 2D figure with unambiguous radius and topology encoding.

    The PNG named by output_path is accompanied by vector PDF/SVG files and a
    ready-to-edit English caption.
    """
    if not is_xz_planar(tree):
        raise ValueError("The planar publication figure requires constant Y coordinates.")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    base_font = max(10.0, min(0.38 * font_size, 18.0))
    _style(base_font)

    (
        root_id,
        children,
        bifurcations,
        terminals,
        generation,
        path_distance,
    ) = _tree_topology(tree)
    points, segments, radii, x_min, z_min = _relative_geometry(tree)
    radius_min = float(np.min(radii))
    radius_max = float(np.max(radii))
    norm = (
        PowerNorm(gamma=0.65, vmin=radius_min, vmax=radius_max)
        if radius_max > radius_min
        else Normalize(vmin=radius_min - 0.5, vmax=radius_max + 0.5)
    )
    cmap = mpl.colormaps[color_map]
    opening_by_node = _bifurcation_opening_angles(tree, children)
    opening_angles = np.asarray(list(opening_by_node.values()), dtype=float)
    edge_lengths = np.asarray(
        [
            float(np.linalg.norm(tree.nodes[child_id].xyz - tree.nodes[parent_id].xyz))
            for parent_id, child_id in tree.edges
        ],
        dtype=float,
    )

    figure = plt.figure(figsize=(13.4, 8.2), constrained_layout=True)
    grid = figure.add_gridspec(
        3,
        4,
        width_ratios=(1.0, 1.0, 1.0, 1.06),
        height_ratios=(1.0, 1.0, 1.0),
    )
    main_ax = figure.add_subplot(grid[:, :3])
    zoom_ax = figure.add_subplot(grid[0, 3])
    histogram_ax = figure.add_subplot(grid[1, 3])
    angle_ax = figure.add_subplot(grid[2, 3])

    network_halo = LineCollection(
        segments,
        colors="#e3e7e9",
        linewidths=4.15,
        capstyle="round",
        joinstyle="round",
        zorder=1,
    )
    main_ax.add_collection(network_halo)
    network = LineCollection(
        segments,
        cmap=cmap,
        norm=norm,
        linewidths=2.35,
        capstyle="round",
        joinstyle="round",
        zorder=2,
    )
    network.set_array(radii)
    main_ax.add_collection(network)

    bifurcation_xy = np.asarray([points[node_id] for node_id in bifurcations])
    terminal_xy = np.asarray([points[node_id] for node_id in terminals])
    root_xy = points[root_id]
    main_ax.scatter(
        bifurcation_xy[:, 0],
        bifurcation_xy[:, 1],
        s=13,
        facecolor="#1f2529",
        edgecolor="white",
        linewidth=0.35,
        zorder=5,
    )
    main_ax.scatter(
        terminal_xy[:, 0],
        terminal_xy[:, 1],
        s=18,
        facecolor="white",
        edgecolor="#1f2529",
        linewidth=0.75,
        zorder=5,
    )
    main_ax.scatter(
        [root_xy[0]],
        [root_xy[1]],
        marker="D",
        s=50,
        facecolor="#d1495b",
        edgecolor="white",
        linewidth=0.8,
        zorder=8,
    )
    _add_flow_arrows(main_ax, tree, points, generation)

    x_span = max(point[0] for point in points.values())
    z_span = max(point[1] for point in points.values())
    tick_step = 500.0
    x_tick_max = np.ceil(x_span / tick_step) * tick_step
    z_tick_max = np.ceil(z_span / tick_step) * tick_step
    main_ax.set_xlim(-0.035 * x_tick_max, x_tick_max + 0.03 * x_tick_max)
    main_ax.set_ylim(-0.035 * z_tick_max, z_tick_max + 0.03 * z_tick_max)
    main_ax.set_xticks(np.arange(0.0, x_tick_max + 1.0, tick_step))
    main_ax.set_yticks(np.arange(0.0, z_tick_max + 1.0, tick_step))
    main_ax.set_aspect("equal", adjustable="box")
    main_ax.set_xlabel(r"$X-X_{\min}$ ($\mathrm{\mu m}$)")
    main_ax.set_ylabel(r"$Z-Z_{\min}$ ($\mathrm{\mu m}$)")
    main_ax.set_title(
        "(a) Planar X\u2013Z vascular tree",
        loc="left",
        fontweight="bold",
        pad=24,
    )
    main_ax.text(
        0.0,
        1.006,
        f"Y = {tree.nodes[root_id].xyz[1]:.0f} \u00b5m; "
        f"median opening = {np.median(opening_angles):.1f}°; "
        f"minimum segment = {np.min(edge_lengths):.1f} \u00b5m",
        transform=main_ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=base_font - 1.0,
        color="#3b4348",
    )
    main_ax.grid(color="#d9dde0", linewidth=0.45, alpha=0.55)

    legend_handles = [
        Line2D(
            [],
            [],
            marker="D",
            linestyle="none",
            markerfacecolor="#d1495b",
            markeredgecolor="white",
            markersize=6,
            label="Root / inlet",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markerfacecolor="#1f2529",
            markeredgecolor="white",
            markersize=5,
            label="Bifurcation",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="#1f2529",
            markersize=5,
            label="Terminal",
        ),
        Line2D(
            [],
            [],
            color="#20252a",
            marker=">",
            markevery=[1],
            markersize=5,
            label="Flow direction",
        ),
    ]
    main_ax.legend(
        handles=legend_handles,
        loc="lower left",
        ncol=2,
        frameon=True,
        facecolor="white",
        edgecolor="#c7cccf",
        framealpha=0.92,
    )

    color_bar = figure.colorbar(
        network,
        ax=main_ax,
        orientation="horizontal",
        fraction=0.035,
        pad=0.06,
        aspect=32,
    )
    color_bar.set_label(r"Local radius, $r$ ($\mathrm{\mu m}$)")
    color_bar.outline.set_linewidth(0.7)

    median_opening = float(np.median(opening_angles))
    selected_id = min(
        opening_by_node,
        key=lambda node_id: abs(opening_by_node[node_id] - median_opening),
    )
    local_edges = [
        edge
        for edge in tree.edges
        if edge[0] == selected_id or edge[1] == selected_id
    ]
    local_segments = np.asarray(
        [[points[parent_id], points[child_id]] for parent_id, child_id in local_edges]
    )
    local_radii = np.asarray(
        [tree.nodes[child_id].radius_um for _, child_id in local_edges]
    )
    local_network = LineCollection(
        local_segments,
        cmap=cmap,
        norm=norm,
        linewidths=4.0,
        capstyle="round",
        zorder=2,
    )
    local_network.set_array(local_radii)
    zoom_ax.add_collection(local_network)
    selected_xy = points[selected_id]
    zoom_ax.scatter(
        [selected_xy[0]],
        [selected_xy[1]],
        s=30,
        facecolor="#1f2529",
        edgecolor="white",
        linewidth=0.6,
        zorder=5,
    )
    local_points = local_segments.reshape(-1, 2)
    local_width = max(float(np.ptp(local_points[:, 0])), 1.0)
    local_height = max(float(np.ptp(local_points[:, 1])), 1.0)
    local_margin = 0.25 * max(local_width, local_height)
    zoom_ax.set_xlim(
        float(np.min(local_points[:, 0])) - local_margin,
        float(np.max(local_points[:, 0])) + local_margin,
    )
    zoom_ax.set_ylim(
        float(np.min(local_points[:, 1])) - local_margin,
        float(np.max(local_points[:, 1])) + local_margin,
    )
    zoom_ax.set_aspect("equal", adjustable="box")
    zoom_ax.set_xticks([])
    zoom_ax.set_yticks([])
    zoom_ax.spines[["top", "right", "bottom", "left"]].set_visible(True)
    zoom_ax.spines[["top", "right", "bottom", "left"]].set_color("#9ca4a8")
    zoom_ax.set_title("(b) Representative bifurcation", loc="left", fontweight="bold")
    _add_opening_angle_annotation(
        zoom_ax,
        selected_xy,
        [points[child_id] for child_id in children[selected_id]],
        opening_by_node[selected_id],
    )
    _add_scale_bar(zoom_ax)

    unique_radii, radius_counts = np.unique(np.round(radii, decimals=8), return_counts=True)
    if len(unique_radii) <= 15:
        bar_width = (
            0.62 * float(np.min(np.diff(unique_radii)))
            if len(unique_radii) > 1
            else 1.0
        )
        histogram_ax.bar(
            unique_radii,
            radius_counts,
            width=bar_width,
            color=cmap(norm(unique_radii)),
            edgecolor="white",
            linewidth=0.7,
        )
    else:
        bin_count = max(6, min(12, int(np.ceil(np.sqrt(len(radii))))))
        histogram_ax.hist(
            radii,
            bins=bin_count,
            color="#3e7f8f",
            edgecolor="white",
            linewidth=0.7,
        )
    median_radius = float(np.median(radii))
    histogram_ax.axvline(
        median_radius,
        color="#d1495b",
        linewidth=1.4,
        linestyle="--",
        label=f"Median = {median_radius:.1f} \u00b5m",
    )
    histogram_ax.set_xlabel(r"Local radius, $r$ ($\mathrm{\mu m}$)")
    histogram_ax.set_ylabel("Segment count")
    histogram_ax.set_title("(c) Radius distribution", loc="left", fontweight="bold")
    histogram_ax.legend(frameon=False, loc="upper right")
    histogram_ax.grid(axis="y", color="#d9dde0", linewidth=0.45)

    angle_bin_count = max(5, min(10, int(np.ceil(np.sqrt(len(opening_angles))))))
    angle_ax.hist(
        opening_angles,
        bins=angle_bin_count,
        color="#6b78b4",
        edgecolor="white",
        linewidth=0.7,
    )
    angle_ax.axvline(
        median_opening,
        color="#d1495b",
        linewidth=1.4,
        linestyle="--",
        label=f"Median = {median_opening:.1f}°",
    )
    angle_ax.set_xlabel("Daughter opening angle (°)")
    angle_ax.set_ylabel("Bifurcation count")
    angle_ax.set_title("(d) Opening-angle distribution", loc="left", fontweight="bold")
    angle_ax.legend(frameon=False, loc="upper left")
    angle_ax.grid(axis="y", color="#d9dde0", linewidth=0.45)

    for suffix, kwargs in (
        (output_path.suffix.lower().lstrip(".") or "png", {"dpi": 300}),
        ("pdf", {}),
        ("svg", {}),
    ):
        target = output_path.with_suffix(f".{suffix}")
        figure.savefig(target, bbox_inches="tight", **kwargs)
    plt.close(figure)

    caption_path = output_path.with_name(f"{output_path.stem}_caption.txt")
    _write_caption(
        caption_path,
        tree,
        float(tree.nodes[root_id].xyz[1]),
        x_min,
        z_min,
        len(bifurcations),
        len(terminals),
        median_opening,
    )
    return [
        output_path,
        output_path.with_suffix(".pdf"),
        output_path.with_suffix(".svg"),
        caption_path,
    ]
