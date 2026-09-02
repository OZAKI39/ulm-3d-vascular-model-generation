"""Draw publication figures for 2-D domain construction and the final CFD mesh.

Both figures are derived from the geometry and field archives produced by
``generate_microbubble_trajectories.py``.  No synthetic vessel geometry or
surrogate grid is introduced.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import ConnectionPatch, Patch, Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import ScaledTranslation

from draw_vascular_grid_construction import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_VESSEL_ARCHIVE,
    PACKAGE_DIR,
    _cell_edges,
    _configure_style,
    _latest_field_archive,
    _representative_bifurcation,
    _segments_in_box,
)


DOMAIN_BASENAME = "two_dimensional_computational_domain"
GRID_BASENAME = "final_cfd_grid"
COMBINED_BASENAME = "two_dimensional_domain_with_local_cfd_grid"
ACCENT = "#C23B4A"
ZOOM_ACCENT = "#0072B2"
INK = "#222222"
MUTED = "#707070"


def _load_inputs(
    field_archive: Path,
    vessel_archive: Path,
) -> tuple[dict[str, np.ndarray | float], dict[str, np.ndarray]]:
    field_keys = {
        "x_coordinates_um",
        "z_coordinates_um",
        "spacing_um",
        "fixed_y_um",
        "lumen_mask",
        "wall_mask",
        "distance_to_wall_um",
        "continuous_wall_start_xz_um",
        "continuous_wall_end_xz_um",
    }
    with np.load(field_archive, allow_pickle=False) as archive:
        missing = sorted(field_keys.difference(archive.files))
        if missing:
            raise KeyError(f"Field archive is missing required arrays: {missing}")
        field: dict[str, np.ndarray | float] = {
            "x_um": np.asarray(archive["x_coordinates_um"], dtype=float),
            "z_um": np.asarray(archive["z_coordinates_um"], dtype=float),
            "spacing_um": float(np.asarray(archive["spacing_um"]).reshape(-1)[0]),
            "fixed_y_um": float(np.asarray(archive["fixed_y_um"]).reshape(-1)[0]),
            "lumen_mask": np.asarray(archive["lumen_mask"], dtype=bool),
            "wall_mask": np.asarray(archive["wall_mask"], dtype=bool),
            "distance_to_wall_um": np.asarray(
                archive["distance_to_wall_um"], dtype=float
            ),
            "wall_start": np.asarray(
                archive["continuous_wall_start_xz_um"], dtype=float
            ),
            "wall_end": np.asarray(
                archive["continuous_wall_end_xz_um"], dtype=float
            ),
        }

    with np.load(vessel_archive, allow_pickle=False) as archive:
        vessel = {
            key: np.asarray(archive[key])
            for key in ("x_p", "x_d", "radius_um", "children_count")
        }
    return field, vessel


def _clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.5, width=0.7)
    ax.grid(False)
    ax.set_aspect("equal", adjustable="box")


def _tight_layout(
    fig: plt.Figure,
    *,
    pad: float = 0.65,
    rect: tuple[float, float, float, float] = (0.0, 0.0, 1.0, 1.0),
) -> None:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="This figure includes Axes that are not compatible with tight_layout",
            category=UserWarning,
        )
        fig.tight_layout(pad=pad, rect=rect)


def draw_two_dimensional_domain(
    field: dict[str, np.ndarray | float],
    vessel: dict[str, np.ndarray],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Visualize how endpoint bounds and padding define the 2-D X-Z domain."""

    _configure_style()
    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    fixed_y_um = float(field["fixed_y_um"])
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)

    vessel_starts = np.asarray(vessel["x_p"], dtype=float)[:, [0, 2]]
    vessel_ends = np.asarray(vessel["x_d"], dtype=float)[:, [0, 2]]
    vessel_segments = np.stack((vessel_starts, vessel_ends), axis=1)
    radii = np.asarray(vessel["radius_um"], dtype=float)
    points = np.vstack((vessel_starts, vessel_ends))
    anatomical_min = points.min(axis=0)
    anatomical_max = points.max(axis=0)
    domain_min = np.asarray([x_edges[0], z_edges[0]], dtype=float)
    domain_max = np.asarray([x_edges[-1], z_edges[-1]], dtype=float)

    radius_norm = Normalize(vmin=float(radii.min()), vmax=float(radii.max()))
    radius_cmap = mpl.colormaps["viridis"]
    line_widths = 0.7 + 3.0 * radii / float(radii.max())

    fig = plt.figure(figsize=(7.5, 5.5))
    spec = fig.add_gridspec(1, 2, width_ratios=(1.0, 0.035), wspace=0.12)
    ax = fig.add_subplot(spec[0, 0])
    cax = fig.add_subplot(spec[0, 1])

    collection = LineCollection(
        vessel_segments,
        array=radii,
        cmap=radius_cmap,
        norm=radius_norm,
        linewidths=line_widths,
        capstyle="round",
        joinstyle="round",
        zorder=4,
    )
    ax.add_collection(collection)

    anatomical_box = Rectangle(
        anatomical_min,
        *(anatomical_max - anatomical_min),
        fill=False,
        edgecolor=MUTED,
        linewidth=1.0,
        linestyle=(0, (4.0, 3.0)),
        zorder=2,
    )
    padded_box = Rectangle(
        padded_min,
        *(padded_max - padded_min),
        fill=False,
        edgecolor=ACCENT,
        linewidth=1.15,
        linestyle=(0, (2.5, 2.0)),
        zorder=2,
    )
    domain_box = Rectangle(
        domain_min,
        *(domain_max - domain_min),
        facecolor="#F4F4F4",
        edgecolor=INK,
        linewidth=1.15,
        zorder=0,
    )
    ax.add_patch(domain_box)
    ax.add_patch(padded_box)
    ax.add_patch(anatomical_box)

    arrow_y = anatomical_min[1] + 0.035 * (anatomical_max[1] - anatomical_min[1])
    ax.annotate(
        "",
        xy=(anatomical_min[0], arrow_y),
        xytext=(padded_min[0], arrow_y),
        arrowprops=dict(arrowstyle="<->", color=ACCENT, lw=1.0),
        zorder=6,
    )
    ax.text(
        0.5 * (padded_min[0] + anatomical_min[0]),
        arrow_y + 28.0,
        rf"$p={configured_padding:g}\ \mathrm{{\mu m}}$",
        ha="center",
        va="bottom",
        color=ACCENT,
        fontsize=10,
        zorder=6,
    )

    domain_width = float(domain_max[0] - domain_min[0])
    domain_height = float(domain_max[1] - domain_min[1])
    ax.annotate(
        "",
        xy=(domain_max[0], domain_max[1] + 34.0),
        xytext=(domain_min[0], domain_max[1] + 34.0),
        arrowprops=dict(arrowstyle="|-|", color=INK, lw=0.85),
        annotation_clip=False,
    )
    ax.text(
        0.5 * (domain_min[0] + domain_max[0]),
        domain_max[1] + 60.0,
        rf"$L_x={domain_width:g}\ \mathrm{{\mu m}}$  ($N_x={x_um.size}$)",
        ha="center",
        va="bottom",
        fontsize=10,
        clip_on=False,
    )
    ax.annotate(
        "",
        xy=(domain_max[0] + 34.0, domain_max[1]),
        xytext=(domain_max[0] + 34.0, domain_min[1]),
        arrowprops=dict(arrowstyle="|-|", color=INK, lw=0.85),
        annotation_clip=False,
    )
    ax.text(
        domain_max[0] + 62.0,
        0.5 * (domain_min[1] + domain_max[1]),
        rf"$L_z={domain_height:g}\ \mathrm{{\mu m}}$  ($N_z={z_um.size}$)",
        ha="left",
        va="center",
        rotation=90,
        fontsize=10,
        clip_on=False,
    )

    info = (
        rf"Projection plane: $y={fixed_y_um:g}\ \mathrm{{\mu m}}$"
        "\n"
        rf"Uniform spacing: $\Delta x=\Delta z={spacing_um:g}\ \mathrm{{\mu m}}$"
    )
    ax.text(
        0.035,
        0.96,
        info,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="#D0D0D0", lw=0.7),
        zorder=8,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=MUTED,
            lw=1.0,
            linestyle=(0, (4.0, 3.0)),
            label="Endpoint envelope",
        ),
        Line2D(
            [0],
            [0],
            color=ACCENT,
            lw=1.15,
            linestyle=(0, (2.5, 2.0)),
            label=rf"$+{configured_padding:g}\ \mathrm{{\mu m}}$ padded bounds",
        ),
        Patch(facecolor="#F4F4F4", edgecolor=INK, label="Final cell-edge domain"),
    ]
    legend = fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.47, -0.005),
        ncol=3,
        frameon=True,
        framealpha=0.95,
        facecolor="white",
        edgecolor="#D0D0D0",
        borderpad=0.45,
        handlelength=2.2,
        columnspacing=1.5,
        fontsize=9,
    )
    legend.get_frame().set_linewidth(0.7)

    margin_x = 0.06 * domain_width
    margin_z = 0.06 * domain_height
    ax.set(
        xlim=(domain_min[0] - margin_x, domain_max[0] + 1.25 * margin_x),
        ylim=(domain_min[1] - 0.08 * domain_height, domain_max[1] + margin_z),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Two-dimensional computational-domain construction",
    )
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    _clean_axes(ax)

    colorbar = fig.colorbar(collection, cax=cax, orientation="vertical")
    colorbar.set_label(r"Vessel radius ($\mathrm{\mu m}$)", fontsize=12)
    colorbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colorbar.outline.set_linewidth(0.7)

    _tight_layout(fig, rect=(0.0, 0.16, 1.0, 1.0))
    png_path = output_dir / f"{DOMAIN_BASENAME}.png"
    pdf_path = output_dir / f"{DOMAIN_BASENAME}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return png_path, pdf_path


def draw_final_cfd_grid(
    field: dict[str, np.ndarray | float],
    vessel: dict[str, np.ndarray],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Visualize the accepted Cartesian CFD mesh at global and cell scales."""

    _configure_style()
    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen_mask"], dtype=bool)
    wall_mask = np.asarray(field["wall_mask"], dtype=bool)
    wall_distance = np.asarray(field["distance_to_wall_um"], dtype=float)
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)
    x_limits = (float(x_edges[0]), float(x_edges[-1]))
    z_limits = (float(z_edges[0]), float(z_edges[-1]))

    vessel_starts = np.asarray(vessel["x_p"], dtype=float)[:, [0, 2]]
    vessel_ends = np.asarray(vessel["x_d"], dtype=float)[:, [0, 2]]
    branch_point, branch_radius = _representative_bifurcation(
        vessel, x_limits, z_limits
    )
    zoom_half_width = max(70.0, 4.0 * branch_radius)
    zoom_bounds = (
        float(branch_point[0] - zoom_half_width),
        float(branch_point[0] + zoom_half_width),
        float(branch_point[1] - zoom_half_width),
        float(branch_point[1] + zoom_half_width),
    )

    ix = np.flatnonzero(
        (x_um >= zoom_bounds[0] - spacing_um)
        & (x_um <= zoom_bounds[1] + spacing_um)
    )
    iz = np.flatnonzero(
        (z_um >= zoom_bounds[2] - spacing_um)
        & (z_um <= zoom_bounds[3] + spacing_um)
    )
    ix0, ix1 = int(ix[0]), int(ix[-1]) + 1
    iz0, iz1 = int(iz[0]), int(iz[-1]) + 1

    distance_values = wall_distance[lumen]
    distance_max = float(np.ceil(np.percentile(distance_values, 99.5) / 5.0) * 5.0)
    distance_norm = Normalize(vmin=0.0, vmax=distance_max, clip=True)
    distance_cmap = mpl.colormaps["viridis"].copy()
    distance_cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    masked_distance = np.ma.masked_where(~lumen, wall_distance)

    fig = plt.figure(figsize=(8.0, 4.6))
    spec = fig.add_gridspec(1, 3, width_ratios=(1.5, 1.0, 0.055), wspace=0.25)
    ax_global = fig.add_subplot(spec[0, 0])
    ax_zoom = fig.add_subplot(spec[0, 1])
    cax = fig.add_subplot(spec[0, 2])

    ax_global.set_facecolor("#F3F3F3")
    major_stride = 100
    ax_global.vlines(
        x_edges[::major_stride],
        z_edges[0],
        z_edges[-1],
        color="#D4D4D4",
        linewidth=0.35,
        zorder=0,
    )
    ax_global.hlines(
        z_edges[::major_stride],
        x_edges[0],
        x_edges[-1],
        color="#D4D4D4",
        linewidth=0.35,
        zorder=0,
    )
    overview = ax_global.imshow(
        masked_distance.T,
        origin="lower",
        extent=(x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]),
        cmap=distance_cmap,
        norm=distance_norm,
        interpolation="nearest",
        aspect="equal",
        rasterized=True,
        zorder=2,
    )
    ax_global.add_collection(
        LineCollection(
            np.stack((wall_start, wall_end), axis=1),
            colors=INK,
            linewidths=0.18,
            alpha=0.85,
            rasterized=True,
            zorder=3,
        )
    )
    ax_global.add_patch(
        Rectangle(
            (zoom_bounds[0], zoom_bounds[2]),
            zoom_bounds[1] - zoom_bounds[0],
            zoom_bounds[3] - zoom_bounds[2],
            fill=False,
            edgecolor=ACCENT,
            linewidth=1.2,
            zorder=5,
        )
    )
    ax_global.text(
        0.03,
        0.04,
        "Uniform Cartesian mesh\n"
        + rf"${x_um.size}\times{z_um.size}$ cells, $\Delta={spacing_um:g}\ \mathrm{{\mu m}}$",
        transform=ax_global.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        linespacing=1.35,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#D0D0D0", lw=0.7),
        zorder=7,
    )
    ax_global.set(
        xlim=x_limits,
        ylim=z_limits,
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Final Cartesian CFD mesh",
    )
    ax_global.xaxis.set_major_locator(MaxNLocator(5))
    ax_global.yaxis.set_major_locator(MaxNLocator(5))
    _clean_axes(ax_global)

    local_x_edges = x_edges[ix0 : ix1 + 1]
    local_z_edges = z_edges[iz0 : iz1 + 1]
    local_lumen = lumen[ix0:ix1, iz0:iz1]
    local_wall = wall_mask[ix0:ix1, iz0:iz1]
    local_distance = np.ma.masked_where(
        ~local_lumen, wall_distance[ix0:ix1, iz0:iz1]
    )
    background = np.zeros(local_lumen.shape, dtype=float)
    ax_zoom.pcolormesh(
        local_x_edges,
        local_z_edges,
        background.T,
        cmap=mpl.colors.ListedColormap(["#F0F0F0"]),
        vmin=0.0,
        vmax=1.0,
        shading="flat",
        edgecolors="#D0D0D0",
        linewidth=0.18,
        antialiased=True,
        rasterized=True,
        zorder=1,
    )
    ax_zoom.pcolormesh(
        local_x_edges,
        local_z_edges,
        local_distance.T,
        cmap=distance_cmap,
        norm=distance_norm,
        shading="flat",
        edgecolors=(1.0, 1.0, 1.0, 0.62),
        linewidth=0.18,
        antialiased=True,
        rasterized=True,
        zorder=2,
    )
    wall_overlay = np.ma.masked_where(~local_wall, np.ones(local_wall.shape))
    ax_zoom.pcolormesh(
        local_x_edges,
        local_z_edges,
        wall_overlay.T,
        cmap=mpl.colors.ListedColormap([(0.76, 0.23, 0.29, 0.28)]),
        vmin=0.0,
        vmax=1.0,
        shading="flat",
        edgecolors="none",
        rasterized=True,
        zorder=3,
    )

    local_wall_segments = _segments_in_box(wall_start, wall_end, zoom_bounds)
    local_centreline_segments = _segments_in_box(
        vessel_starts, vessel_ends, zoom_bounds
    )
    ax_zoom.add_collection(
        LineCollection(
            local_wall_segments,
            colors=INK,
            linewidths=0.8,
            zorder=5,
        )
    )
    ax_zoom.add_collection(
        LineCollection(
            local_centreline_segments,
            colors="white",
            linewidths=1.0,
            linestyles=(0, (3.0, 2.0)),
            zorder=6,
        )
    )
    ax_zoom.scatter(
        branch_point[0],
        branch_point[1],
        s=15,
        facecolor="white",
        edgecolor=INK,
        linewidth=0.7,
        zorder=7,
    )

    cell_i = int(np.argmin(np.abs(x_um - branch_point[0])))
    cell_j = int(np.argmin(np.abs(z_um - branch_point[1])))
    representative_x = float(x_um[cell_i])
    representative_z = float(z_um[cell_j])
    ax_zoom.scatter(
        representative_x,
        representative_z,
        s=9,
        marker="o",
        facecolor=INK,
        edgecolor="white",
        linewidth=0.45,
        zorder=8,
    )
    ax_zoom.annotate(
        "cell-centred variable",
        xy=(representative_x, representative_z),
        xytext=(0.04, 0.92),
        textcoords="axes fraction",
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9),
        arrowprops=dict(arrowstyle="-", color=INK, lw=0.7),
        zorder=9,
    )
    ax_zoom.set(
        xlim=(zoom_bounds[0], zoom_bounds[1]),
        ylim=(zoom_bounds[2], zoom_bounds[3]),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Cell-resolved CFD grid\nat a bifurcation",
    )
    ax_zoom.xaxis.set_major_locator(MaxNLocator(4))
    ax_zoom.yaxis.set_major_locator(MaxNLocator(4))
    _clean_axes(ax_zoom)

    legend_handles = [
        Patch(facecolor="#F0F0F0", edgecolor="#BEBEBE", label="Solid-domain cell"),
        Patch(facecolor=ACCENT, alpha=0.28, edgecolor=ACCENT, label="Wall-adjacent lumen cell"),
        Line2D([0], [0], color=INK, lw=1.0, label="Continuous vessel wall"),
        Line2D(
            [0],
            [0],
            color="#777777",
            lw=1.0,
            linestyle=(0, (3.0, 2.0)),
            label="Vessel centreline",
        ),
    ]
    legend = fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.48, 0.045),
        ncol=2,
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="#D0D0D0",
        borderpad=0.45,
        handlelength=1.8,
        columnspacing=1.8,
        fontsize=9,
    )
    legend.get_frame().set_linewidth(0.7)

    colorbar = fig.colorbar(overview, cax=cax, orientation="vertical")
    colorbar.set_label(r"Distance to wall ($\mathrm{\mu m}$)", fontsize=12)
    colorbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colorbar.outline.set_linewidth(0.7)

    for label, axis in (("a", ax_global), ("b", ax_zoom)):
        axis.text(
            -0.13,
            1.035,
            f"({label})",
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=12,
            fontweight="bold",
            clip_on=False,
        )

    _tight_layout(fig, pad=0.55, rect=(0.0, 0.04, 1.0, 1.0))
    png_path = output_dir / f"{GRID_BASENAME}.png"
    pdf_path = output_dir / f"{GRID_BASENAME}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return png_path, pdf_path


def draw_domain_with_local_cfd_grid(
    field: dict[str, np.ndarray | float],
    vessel: dict[str, np.ndarray],
    output_dir: Path,
) -> tuple[Path, Path]:
    """Draw the full 2-D domain with the CFD mesh shown only as an inset."""

    _configure_style()
    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    fixed_y_um = float(field["fixed_y_um"])
    lumen = np.asarray(field["lumen_mask"], dtype=bool)
    wall_mask = np.asarray(field["wall_mask"], dtype=bool)
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)

    vessel_starts = np.asarray(vessel["x_p"], dtype=float)[:, [0, 2]]
    vessel_ends = np.asarray(vessel["x_d"], dtype=float)[:, [0, 2]]
    vessel_segments = np.stack((vessel_starts, vessel_ends), axis=1)
    radii = np.asarray(vessel["radius_um"], dtype=float)
    points = np.vstack((vessel_starts, vessel_ends))
    anatomical_min = points.min(axis=0)
    anatomical_max = points.max(axis=0)
    configured_padding = float(anatomical_min[0] - x_um[0])
    padded_min = anatomical_min - configured_padding
    padded_max = anatomical_max + configured_padding
    domain_min = np.asarray([x_edges[0], z_edges[0]], dtype=float)
    domain_max = np.asarray([x_edges[-1], z_edges[-1]], dtype=float)
    x_limits = (float(domain_min[0]), float(domain_max[0]))
    z_limits = (float(domain_min[1]), float(domain_max[1]))

    branch_point, branch_radius = _representative_bifurcation(
        vessel,
        x_limits,
        z_limits,
    )
    zoom_half_width = max(70.0, 4.0 * branch_radius)
    zoom_bounds = (
        float(branch_point[0] - zoom_half_width),
        float(branch_point[0] + zoom_half_width),
        float(branch_point[1] - zoom_half_width),
        float(branch_point[1] + zoom_half_width),
    )
    ix = np.flatnonzero(
        (x_um >= zoom_bounds[0] - spacing_um)
        & (x_um <= zoom_bounds[1] + spacing_um)
    )
    iz = np.flatnonzero(
        (z_um >= zoom_bounds[2] - spacing_um)
        & (z_um <= zoom_bounds[3] + spacing_um)
    )
    if ix.size < 3 or iz.size < 3:
        raise ValueError("The selected bifurcation window falls outside the CFD grid.")
    ix0, ix1 = int(ix[0]), int(ix[-1]) + 1
    iz0, iz1 = int(iz[0]), int(iz[-1]) + 1

    radius_norm = Normalize(vmin=float(radii.min()), vmax=float(radii.max()))
    radius_cmap = mpl.colormaps["viridis"]
    line_widths = 0.7 + 3.0 * radii / float(radii.max())

    fig = plt.figure(figsize=(8.0, 5.9))
    spec = fig.add_gridspec(1, 2, width_ratios=(1.0, 0.035), wspace=0.11)
    ax = fig.add_subplot(spec[0, 0])
    cax = fig.add_subplot(spec[0, 1])

    domain_box = Rectangle(
        domain_min,
        *(domain_max - domain_min),
        facecolor="#F4F4F4",
        edgecolor=INK,
        linewidth=1.15,
        zorder=0,
    )
    anatomical_box = Rectangle(
        anatomical_min,
        *(anatomical_max - anatomical_min),
        fill=False,
        edgecolor=MUTED,
        linewidth=1.0,
        linestyle=(0, (4.0, 3.0)),
        zorder=2,
    )
    ax.add_patch(domain_box)
    ax.add_patch(anatomical_box)
    collection = LineCollection(
        vessel_segments,
        array=radii,
        cmap=radius_cmap,
        norm=radius_norm,
        linewidths=line_widths,
        capstyle="round",
        joinstyle="round",
        zorder=4,
    )
    ax.add_collection(collection)

    roi_box = Rectangle(
        (zoom_bounds[0], zoom_bounds[2]),
        zoom_bounds[1] - zoom_bounds[0],
        zoom_bounds[3] - zoom_bounds[2],
        fill=False,
        edgecolor=ZOOM_ACCENT,
        linewidth=1.35,
        zorder=7,
    )
    ax.add_patch(roi_box)
    ax.text(
        zoom_bounds[0],
        zoom_bounds[3] + 22.0,
        "Local mesh ROI",
        color=ZOOM_ACCENT,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="bottom",
        zorder=8,
    )

    domain_width = float(domain_max[0] - domain_min[0])
    domain_height = float(domain_max[1] - domain_min[1])
    info_box = AnchoredText(
        (
            rf"Grid resolution: $N_x={x_um.size}$, $N_z={z_um.size}$"
            "\n"
            rf"Uniform spacing: $\Delta x=\Delta z={spacing_um:g}\ "
            r"\mathrm{\mu m}$"
        ),
        loc="upper left",
        prop={"size": 10, "linespacing": 1.35},
        frameon=True,
        bbox_to_anchor=(domain_min[0], domain_max[1]),
        bbox_transform=(
            ax.transData
            + ScaledTranslation(
                3.2 / 72.0,
                -3.2 / 72.0,
                fig.dpi_scale_trans,
            )
        ),
        borderpad=0.0,
    )
    info_box.patch.set_boxstyle("round,pad=0.30")
    info_box.patch.set_facecolor("white")
    info_box.patch.set_edgecolor("#D0D0D0")
    info_box.patch.set_linewidth(0.7)
    info_box.set_zorder(9)
    ax.add_artist(info_box)

    # The mesh is intentionally an inset rather than a second global panel.
    ax_zoom = ax.inset_axes(
        [0.585, 0.000, 0.300, 0.360],
        facecolor="white",
        zorder=12,
    )
    local_lumen = lumen[ix0:ix1, iz0:iz1]
    cell_class = np.zeros(local_lumen.shape, dtype=int)
    cell_class[local_lumen] = 1
    cell_cmap = mpl.colors.ListedColormap(
        [
            "#F4F4F4",
            "#A8D3E8",
        ]
    )
    cell_norm = mpl.colors.BoundaryNorm((-0.5, 0.5, 1.5), 2)
    ax_zoom.pcolormesh(
        x_edges[ix0 : ix1 + 1],
        z_edges[iz0 : iz1 + 1],
        cell_class.T,
        cmap=cell_cmap,
        norm=cell_norm,
        shading="flat",
        edgecolors=(0.62, 0.62, 0.62, 0.34),
        linewidth=0.10,
        antialiased=True,
        rasterized=True,
        zorder=1,
    )
    local_walls = _segments_in_box(wall_start, wall_end, zoom_bounds)
    ax_zoom.add_collection(
        LineCollection(
            local_walls,
            colors=INK,
            linewidths=1.02,
            zorder=4,
        )
    )
    ax_zoom.scatter(
        branch_point[0],
        branch_point[1],
        s=13,
        facecolor="white",
        edgecolor=INK,
        linewidth=0.65,
        zorder=5,
    )
    ax_zoom.set(
        xlim=(zoom_bounds[0], zoom_bounds[1]),
        ylim=(zoom_bounds[2], zoom_bounds[3]),
    )
    ax_zoom.set_xticks([])
    ax_zoom.set_yticks([])
    ax_zoom.grid(False)
    ax_zoom.set_aspect("equal", adjustable="box")
    for spine in ax_zoom.spines.values():
        spine.set_visible(True)
        spine.set_color(ZOOM_ACCENT)
        spine.set_linewidth(1.0)
    inset_legend = ax_zoom.legend(
        handles=[
            Patch(facecolor="#F4F4F4", edgecolor="#B5B5B5", label="Solid"),
            Patch(facecolor="#A8D3E8", edgecolor="#7FA9C2", label="Lumen"),
        ],
        loc="lower right",
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="#C9C9C9",
        fontsize=7.6,
        borderpad=0.35,
        handlelength=1.1,
        handletextpad=0.45,
        labelspacing=0.30,
    )
    inset_legend.get_frame().set_linewidth(0.6)

    for x_value, inset_x in (
        (zoom_bounds[0], 0.02),
        (zoom_bounds[1], 0.98),
    ):
        fig.add_artist(
            ConnectionPatch(
                xyA=(x_value, zoom_bounds[2]),
                coordsA=ax.transData,
                xyB=(inset_x, 1.0),
                coordsB=ax_zoom.transAxes,
                color=ZOOM_ACCENT,
                linewidth=0.75,
                alpha=0.72,
                zorder=11,
                clip_on=False,
            )
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=MUTED,
            lw=1.0,
            linestyle=(0, (4.0, 3.0)),
            label="Endpoint envelope",
        ),
        Patch(
            facecolor="#F4F4F4",
            edgecolor=INK,
            label="Computational domain",
        ),
    ]
    legend = fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.47, 0.985),
        ncol=2,
        frameon=True,
        framealpha=0.95,
        facecolor="white",
        edgecolor="#D0D0D0",
        borderpad=0.42,
        handlelength=2.0,
        columnspacing=1.4,
        fontsize=9,
    )
    legend.get_frame().set_linewidth(0.7)

    ax.set(
        xlim=(domain_min[0], domain_max[0]),
        ylim=(domain_min[1], domain_max[1]),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Two-dimensional computational domain with local CFD mesh",
    )
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    _clean_axes(ax)

    colorbar = fig.colorbar(collection, cax=cax, orientation="vertical")
    colorbar.set_label(r"Vessel radius ($\mathrm{\mu m}$)", fontsize=12)
    colorbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colorbar.outline.set_linewidth(0.7)

    _tight_layout(fig, pad=0.60, rect=(0.0, 0.0, 1.0, 0.89))
    png_path = output_dir / f"{COMBINED_BASENAME}.png"
    pdf_path = output_dir / f"{COMBINED_BASENAME}.pdf"
    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.04,
        facecolor="white",
    )
    fig.savefig(
        pdf_path,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.04,
        facecolor="white",
    )
    plt.close(fig)
    return png_path, pdf_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw the full 2-D computational domain with the final CFD grid "
            "shown only as a local inset."
        )
    )
    parser.add_argument(
        "--field",
        type=Path,
        default=None,
        help="Field archive; defaults to the newest complete timestamped result.",
    )
    parser.add_argument(
        "--vessels",
        type=Path,
        default=DEFAULT_VESSEL_ARCHIVE,
        help="Exported vessel topology archive.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for 300-DPI PNG and PDF outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    field_archive = (
        args.field.resolve()
        if args.field is not None
        else _latest_field_archive(PACKAGE_DIR / "results")
    )
    vessel_archive = args.vessels.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    field, vessel = _load_inputs(field_archive, vessel_archive)
    combined_png, combined_pdf = draw_domain_with_local_cfd_grid(
        field,
        vessel,
        output_dir,
    )
    print(f"Field archive: {field_archive}")
    print(f"Saved combined PNG: {combined_png}")
    print(f"Saved combined PDF: {combined_pdf}")


if __name__ == "__main__":
    main()
