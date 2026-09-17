"""Publication figure for a genuinely three-dimensional vascular SWC tree."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np


def _tree_topology(tree):
    children = {node_id: [] for node_id in tree.nodes}
    roots = []
    for node_id, node in tree.nodes.items():
        if node.parent_id < 0:
            roots.append(node_id)
        else:
            children[node.parent_id].append(node_id)
    if len(roots) != 1:
        raise ValueError(
            f"A publication tree requires one root; found {len(roots)}."
        )
    bifurcations = [
        node_id for node_id, child_ids in children.items() if len(child_ids) > 1
    ]
    terminals = [
        node_id for node_id, child_ids in children.items() if not child_ids
    ]
    return roots[0], bifurcations, terminals


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
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def _relative_geometry(tree):
    coordinates = tree.coordinates
    origin = np.min(coordinates, axis=0)
    points = {
        node_id: node.xyz.astype(float) - origin
        for node_id, node in tree.nodes.items()
    }
    segments = np.asarray(
        [[points[parent_id], points[child_id]] for parent_id, child_id in tree.edges],
        dtype=float,
    )
    radii = np.asarray(
        [tree.nodes[child_id].radius_um for _, child_id in tree.edges],
        dtype=float,
    )
    return points, segments, radii, origin


def _equal_3d_limits(ax, coordinates):
    minima = np.min(coordinates, axis=0)
    maxima = np.max(coordinates, axis=0)
    centers = 0.5 * (minima + maxima)
    half_span = 0.53 * max(float(np.max(maxima - minima)), 1.0)
    ax.set_xlim(centers[0] - half_span, centers[0] + half_span)
    ax.set_ylim(centers[1] - half_span, centers[1] + half_span)
    ax.set_zlim(centers[2] - half_span, centers[2] + half_span)
    ax.set_box_aspect((1.0, 1.0, 1.0))


def _projection_panel(ax, segments, radii, axes, labels, cmap, norm, title):
    projected = segments[:, :, axes]
    collection = LineCollection(
        projected,
        cmap=cmap,
        norm=norm,
        linewidths=1.7,
        capstyle="round",
        joinstyle="round",
    )
    collection.set_array(radii)
    ax.add_collection(collection)
    flat = projected.reshape(-1, 2)
    spans = np.ptp(flat, axis=0)
    margins = np.maximum(0.04 * spans, 1.0)
    ax.set_xlim(float(np.min(flat[:, 0]) - margins[0]), float(np.max(flat[:, 0]) + margins[0]))
    ax.set_ylim(float(np.min(flat[:, 1]) - margins[1]), float(np.max(flat[:, 1]) + margins[1]))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(color="#d9dde0", linewidth=0.45, alpha=0.55)
    return collection


def _write_caption(
    path,
    tree,
    origin,
    spans,
    bifurcation_count,
    terminal_count,
):
    caption = (
        "Volumetric vascular-tree morphology and orthogonal projections. "
        "(a) Three-dimensional isometric view of the genuine volumetric network; "
        "segment color represents local radius. The diamond marks the root, "
        "filled circles mark bifurcations, and open circles mark terminals. "
        "(b-d) Orthographic X-Y, X-Z, and Y-Z projections of the same unmodified "
        "coordinates. Coordinates are shown relative to "
        f"(Xmin, Ymin, Zmin) = ({origin[0]:.1f}, {origin[1]:.1f}, "
        f"{origin[2]:.1f}) µm. Coordinate spans are "
        f"({spans[0]:.1f}, {spans[1]:.1f}, {spans[2]:.1f}) µm. "
        f"The tree contains {tree.edge_count} segments, {bifurcation_count} "
        f"bifurcations, and {terminal_count} terminal nodes."
    )
    path.write_text(caption + "\n", encoding="utf-8")


def render_volumetric_publication_figure(
    tree,
    output_path,
    color_map="viridis",
    font_size=32,
):
    """
    Save a 3D isometric panel plus three quantitative orthogonal projections.

    The PNG named by ``output_path`` is accompanied by PDF, SVG, and an
    editable caption text file.
    """
    coordinates = tree.coordinates
    if coordinates.size == 0:
        raise ValueError("The volumetric publication figure requires a non-empty tree.")
    spans = np.ptp(coordinates, axis=0)
    reference_span = max(float(spans[0]), float(spans[2]), 1.0)
    if float(spans[1]) <= 1.0e-6 * reference_span:
        raise ValueError(
            "The volumetric publication figure requires non-constant Y coordinates."
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    base_font = max(10.0, min(0.36 * font_size, 17.0))
    _style(base_font)

    root_id, bifurcations, terminals = _tree_topology(tree)
    points, segments, radii, origin = _relative_geometry(tree)
    norm = Normalize(vmin=float(np.min(radii)), vmax=float(np.max(radii)))
    cmap = mpl.colormaps[color_map]

    figure = plt.figure(figsize=(13.8, 9.4), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, width_ratios=(1.35, 1.0, 1.0))
    ax3d = figure.add_subplot(grid[:, 0], projection="3d")
    xy_ax = figure.add_subplot(grid[0, 1])
    xz_ax = figure.add_subplot(grid[0, 2])
    yz_ax = figure.add_subplot(grid[1, 1])
    info_ax = figure.add_subplot(grid[1, 2])

    network = Line3DCollection(
        segments,
        cmap=cmap,
        norm=norm,
        linewidths=2.1,
    )
    network.set_array(radii)
    ax3d.add_collection3d(network)
    all_relative = np.vstack(list(points.values()))
    _equal_3d_limits(ax3d, all_relative)

    bif_xyz = np.vstack([points[node_id] for node_id in bifurcations])
    term_xyz = np.vstack([points[node_id] for node_id in terminals])
    root_xyz = points[root_id]
    ax3d.scatter(
        bif_xyz[:, 0], bif_xyz[:, 1], bif_xyz[:, 2],
        s=12, c="#1f2529", depthshade=False,
    )
    ax3d.scatter(
        term_xyz[:, 0], term_xyz[:, 1], term_xyz[:, 2],
        s=15, facecolors="white", edgecolors="#1f2529",
        linewidths=0.65, depthshade=False,
    )
    ax3d.scatter(
        [root_xyz[0]], [root_xyz[1]], [root_xyz[2]],
        marker="D", s=48, c="#d1495b", edgecolors="white",
        linewidths=0.8, depthshade=False,
    )
    ax3d.set_xlabel(r"$X-X_{\min}$ ($\mathrm{\mu m}$)")
    ax3d.set_ylabel(r"$Y-Y_{\min}$ ($\mathrm{\mu m}$)")
    ax3d.set_zlabel(r"$Z-Z_{\min}$ ($\mathrm{\mu m}$)")
    ax3d.set_title("(a) Volumetric vascular tree", loc="left", fontweight="bold")
    ax3d.view_init(elev=24.0, azim=-52.0)
    ax3d.grid(True, linewidth=0.4, alpha=0.45)

    _projection_panel(
        xy_ax, segments, radii, (0, 1),
        (r"$X-X_{\min}$ ($\mathrm{\mu m}$)", r"$Y-Y_{\min}$ ($\mathrm{\mu m}$)"),
        cmap, norm, "(b) X–Y projection",
    )
    _projection_panel(
        xz_ax, segments, radii, (0, 2),
        (r"$X-X_{\min}$ ($\mathrm{\mu m}$)", r"$Z-Z_{\min}$ ($\mathrm{\mu m}$)"),
        cmap, norm, "(c) X–Z projection",
    )
    _projection_panel(
        yz_ax, segments, radii, (1, 2),
        (r"$Y-Y_{\min}$ ($\mathrm{\mu m}$)", r"$Z-Z_{\min}$ ($\mathrm{\mu m}$)"),
        cmap, norm, "(d) Y–Z projection",
    )

    info_ax.axis("off")
    info_ax.set_title("(e) Geometry summary", loc="left", fontweight="bold")
    summary = (
        f"Segments: {tree.edge_count}\n"
        f"Bifurcations: {len(bifurcations)}\n"
        f"Terminals: {len(terminals)}\n\n"
        f"Span X: {spans[0]:.1f} µm\n"
        f"Span Y: {spans[1]:.1f} µm\n"
        f"Span Z: {spans[2]:.1f} µm\n\n"
        f"Radius range:\n{np.min(radii):.2f}–{np.max(radii):.2f} µm"
    )
    info_ax.text(0.04, 0.94, summary, va="top", ha="left", linespacing=1.5)

    color_bar = figure.colorbar(
        network,
        ax=[ax3d, xy_ax, xz_ax, yz_ax],
        orientation="horizontal",
        fraction=0.035,
        pad=0.04,
        aspect=38,
    )
    color_bar.set_label(r"Local radius, $r$ ($\mathrm{\mu m}$)")

    destinations = [
        output_path,
        output_path.with_suffix(".pdf"),
        output_path.with_suffix(".svg"),
    ]
    for destination in destinations:
        figure.savefig(destination, dpi=300, bbox_inches="tight")
    plt.close(figure)

    caption_path = output_path.with_name(output_path.stem + "_caption.txt")
    _write_caption(
        caption_path,
        tree,
        origin,
        spans,
        len(bifurcations),
        len(terminals),
    )
    return [*destinations, caption_path]
