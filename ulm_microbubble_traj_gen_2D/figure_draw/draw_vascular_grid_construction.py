"""Publication-quality visualization of the vascular CFD-grid construction.

The figure uses the numerical artifacts written by
``generate_microbubble_trajectories.py``.  It combines the accepted full-field
velocity map with a cell-resolved view of one representative bifurcation,
showing the centreline, continuous wall, and Cartesian sampling grid together.

Example
-------
python ulm_microbubble_traj_gen/figure_draw/draw_vascular_grid_construction.py
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent
PROJECT_DIR = PACKAGE_DIR.parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "outputs"
DEFAULT_VESSEL_ARCHIVE = (
    PROJECT_DIR
    / "ulm_vascular_model_generator"
    / "vessel_swc_models"
    / "ulm_xz_planar_dcco_tree_seed_105.vessels.npz"
)


def _configure_style() -> None:
    """Apply a restrained, journal-compatible Matplotlib style."""

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "axes.titleweight": "normal",
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "legend.title_fontsize": 12,
            "axes.linewidth": 0.8,
            "axes.edgecolor": "#333333",
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _latest_field_archive(results_dir: Path) -> Path:
    """Return the newest complete velocity-field archive."""

    timestamped = [
        path
        for path in results_dir.glob("*/velocity_and_wall_shear_field.npz")
        if re.match(r"^\d{8}_\d{6}", path.parent.name)
    ]
    candidates = sorted(
        timestamped,
        key=lambda path: path.parent.name,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No velocity_and_wall_shear_field.npz found below {results_dir}."
        )
    continuous_keys = {
        "continuous_wall_start_xz_um",
        "continuous_wall_end_xz_um",
    }
    for candidate in candidates:
        with np.load(candidate, allow_pickle=False) as archive:
            if continuous_keys.issubset(archive.files):
                return candidate
    raise FileNotFoundError(
        "No timestamped field archive contains the continuous-wall arrays "
        f"required for this figure below {results_dir}."
    )


def _cell_edges(coordinates: np.ndarray, spacing_um: float) -> np.ndarray:
    """Convert regularly spaced cell centres to cell-edge coordinates."""

    centres = np.asarray(coordinates, dtype=float)
    return np.concatenate(
        (
            np.asarray([centres[0] - 0.5 * spacing_um]),
            centres + 0.5 * spacing_um,
        )
    )


def _representative_bifurcation(
    vessel_data: dict[str, np.ndarray],
    x_limits: tuple[float, float],
    z_limits: tuple[float, float],
) -> tuple[np.ndarray, float]:
    """Select a large, central bifurcation for the cell-scale panel."""

    is_branch = np.asarray(vessel_data["children_count"]) >= 2
    branch_indices = np.flatnonzero(is_branch)
    if branch_indices.size == 0:
        raise ValueError("The supplied vascular archive contains no bifurcation.")

    points = np.asarray(vessel_data["x_d"], dtype=float)[branch_indices][:, [0, 2]]
    radii = np.asarray(vessel_data["radius_um"], dtype=float)[branch_indices]
    centre = np.asarray(
        [0.5 * sum(x_limits), 0.5 * sum(z_limits)], dtype=float
    )
    half_extent = np.asarray(
        [0.5 * (x_limits[1] - x_limits[0]), 0.5 * (z_limits[1] - z_limits[0])],
        dtype=float,
    )
    normalized_distance = np.linalg.norm((points - centre) / half_extent, axis=1)
    centrality = np.clip(1.0 - normalized_distance / np.sqrt(2.0), 0.0, 1.0)
    radius_scale = radii / max(float(radii.max()), np.finfo(float).eps)
    score = 0.65 * radius_scale + 0.35 * centrality
    selected = int(np.argmax(score))
    return points[selected], float(radii[selected])


def _segments_in_box(
    starts: np.ndarray,
    ends: np.ndarray,
    bounds: tuple[float, float, float, float],
) -> np.ndarray:
    """Filter line segments by a conservative axis-aligned bounding box."""

    xmin, xmax, zmin, zmax = bounds
    segment_min = np.minimum(starts, ends)
    segment_max = np.maximum(starts, ends)
    keep = (
        (segment_max[:, 0] >= xmin)
        & (segment_min[:, 0] <= xmax)
        & (segment_max[:, 1] >= zmin)
        & (segment_min[:, 1] <= zmax)
    )
    return np.stack((starts[keep], ends[keep]), axis=1)


def _format_log_tick(value: float, _: int) -> str:
    return f"{value:g}"


def draw_figure(
    field_archive: Path,
    vessel_archive: Path,
    output_dir: Path,
    basename: str = "vascular_grid_construction",
) -> tuple[Path, Path]:
    """Draw and save the vascular-grid construction figure."""

    _configure_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    required_field_keys = {
        "x_coordinates_um",
        "z_coordinates_um",
        "spacing_um",
        "lumen_mask",
        "speed_um_s",
        "continuous_wall_start_xz_um",
        "continuous_wall_end_xz_um",
    }
    with np.load(field_archive, allow_pickle=False) as archive:
        missing = sorted(required_field_keys.difference(archive.files))
        if missing:
            raise KeyError(f"Field archive is missing required arrays: {missing}")
        x_um = np.asarray(archive["x_coordinates_um"], dtype=float)
        z_um = np.asarray(archive["z_coordinates_um"], dtype=float)
        spacing_um = float(np.asarray(archive["spacing_um"]).reshape(-1)[0])
        lumen_mask = np.asarray(archive["lumen_mask"], dtype=bool)
        speed_mm_s = np.asarray(archive["speed_um_s"], dtype=float) / 1000.0
        wall_start = np.asarray(
            archive["continuous_wall_start_xz_um"], dtype=float
        )
        wall_end = np.asarray(archive["continuous_wall_end_xz_um"], dtype=float)

    with np.load(vessel_archive, allow_pickle=False) as archive:
        vessel_data = {
            key: np.asarray(archive[key])
            for key in ("x_p", "x_d", "radius_um", "children_count")
        }

    if lumen_mask.shape != (x_um.size, z_um.size):
        raise ValueError("The lumen mask shape does not match the stored grid axes.")
    if speed_mm_s.shape != lumen_mask.shape:
        raise ValueError("The speed field shape does not match the lumen mask.")

    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)
    x_limits = (float(x_edges[0]), float(x_edges[-1]))
    z_limits = (float(z_edges[0]), float(z_edges[-1]))

    positive_speed = speed_mm_s[lumen_mask & (speed_mm_s > 0.0)]
    if positive_speed.size == 0:
        raise ValueError("The field archive contains no positive lumen velocity.")
    vmin = max(0.05, float(np.percentile(positive_speed, 0.5)))
    vmax = float(np.ceil(np.percentile(positive_speed, 99.8) * 2.0) / 2.0)
    if vmax <= vmin:
        vmax = float(positive_speed.max())
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=True)
    cmap = mpl.colormaps["cividis"].copy()
    cmap.set_bad("white")
    masked_speed = np.ma.masked_where(~lumen_mask, speed_mm_s)

    branch_point, branch_radius = _representative_bifurcation(
        vessel_data, x_limits, z_limits
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

    vessel_starts = np.asarray(vessel_data["x_p"], dtype=float)[:, [0, 2]]
    vessel_ends = np.asarray(vessel_data["x_d"], dtype=float)[:, [0, 2]]
    centreline_segments = np.stack((vessel_starts, vessel_ends), axis=1)
    wall_segments = np.stack((wall_start, wall_end), axis=1)

    fig = plt.figure(figsize=(7.8, 4.25))
    grid_spec = fig.add_gridspec(
        1,
        3,
        width_ratios=(1.55, 1.0, 0.055),
        wspace=0.26,
    )
    ax_network = fig.add_subplot(grid_spec[0, 0])
    ax_zoom = fig.add_subplot(grid_spec[0, 1])
    cax = fig.add_subplot(grid_spec[0, 2])

    image = ax_network.imshow(
        masked_speed.T,
        origin="lower",
        extent=(x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]),
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
        rasterized=True,
    )
    ax_network.add_collection(
        LineCollection(
            wall_segments,
            colors="#202020",
            linewidths=0.18,
            alpha=0.85,
            zorder=3,
            rasterized=True,
        )
    )
    ax_network.add_collection(
        LineCollection(
            centreline_segments,
            colors="white",
            linewidths=0.28,
            alpha=0.48,
            zorder=4,
        )
    )
    zoom_box = Rectangle(
        (zoom_bounds[0], zoom_bounds[2]),
        zoom_bounds[1] - zoom_bounds[0],
        zoom_bounds[3] - zoom_bounds[2],
        fill=False,
        edgecolor="#C23B4A",
        linewidth=1.25,
        zorder=6,
    )
    ax_network.add_patch(zoom_box)
    ax_network.set(
        xlim=x_limits,
        ylim=z_limits,
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Continuous vascular lumen and CFD field",
    )
    ax_network.set_aspect("equal", adjustable="box")
    ax_network.xaxis.set_major_locator(MaxNLocator(5))
    ax_network.yaxis.set_major_locator(MaxNLocator(5))

    zoom_x_edges = x_edges[ix0 : ix1 + 1]
    zoom_z_edges = z_edges[iz0 : iz1 + 1]
    zoom_speed = np.ma.masked_where(
        ~lumen_mask[ix0:ix1, iz0:iz1],
        speed_mm_s[ix0:ix1, iz0:iz1],
    )
    background = np.zeros((ix1 - ix0, iz1 - iz0), dtype=float)
    ax_zoom.pcolormesh(
        zoom_x_edges,
        zoom_z_edges,
        background.T,
        cmap=ListedColormap(["#F2F2F2"]),
        vmin=0.0,
        vmax=1.0,
        shading="flat",
        edgecolors="#D5D5D5",
        linewidth=0.18,
        antialiased=True,
        rasterized=True,
        zorder=1,
    )
    ax_zoom.pcolormesh(
        zoom_x_edges,
        zoom_z_edges,
        zoom_speed.T,
        cmap=cmap,
        norm=norm,
        shading="flat",
        edgecolors=(1.0, 1.0, 1.0, 0.58),
        linewidth=0.18,
        antialiased=True,
        rasterized=True,
        zorder=2,
    )

    local_wall_segments = _segments_in_box(wall_start, wall_end, zoom_bounds)
    local_centreline_segments = _segments_in_box(
        vessel_starts, vessel_ends, zoom_bounds
    )
    ax_zoom.add_collection(
        LineCollection(
            local_wall_segments,
            colors="#111111",
            linewidths=0.8,
            zorder=4,
        )
    )
    ax_zoom.add_collection(
        LineCollection(
            local_centreline_segments,
            colors="white",
            linewidths=1.0,
            linestyles=(0, (3.0, 2.0)),
            alpha=0.95,
            zorder=5,
        )
    )
    ax_zoom.scatter(
        branch_point[0],
        branch_point[1],
        s=16,
        facecolor="white",
        edgecolor="#111111",
        linewidth=0.7,
        zorder=6,
    )
    ax_zoom.set(
        xlim=(zoom_bounds[0], zoom_bounds[1]),
        ylim=(zoom_bounds[2], zoom_bounds[3]),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Cell-centred rasterization\n"
        + rf"($\Delta x=\Delta z={spacing_um:g}\ \mathrm{{\mu m}}$)",
    )
    ax_zoom.set_aspect("equal", adjustable="box")
    ax_zoom.xaxis.set_major_locator(MaxNLocator(4))
    ax_zoom.yaxis.set_major_locator(MaxNLocator(4))

    annotation_box = dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.88)
    local_wall_midpoints = local_wall_segments.mean(axis=1)
    wall_target = local_wall_midpoints[
        int(np.argmin(np.linalg.norm(local_wall_midpoints - branch_point, axis=1)))
    ]
    ax_zoom.annotate(
        "continuous wall",
        xy=wall_target,
        xytext=(0.04, 0.91),
        textcoords="axes fraction",
        ha="left",
        va="top",
        fontsize=9,
        bbox=annotation_box,
        arrowprops=dict(arrowstyle="-", color="#222222", lw=0.7),
        zorder=8,
    )
    ax_zoom.plot(
        [0.05, 0.13],
        [0.07, 0.07],
        transform=ax_zoom.transAxes,
        color="white",
        linewidth=1.3,
        linestyle=(0, (3.0, 2.0)),
        zorder=8,
    )
    ax_zoom.text(
        0.15,
        0.07,
        "vessel centreline",
        transform=ax_zoom.transAxes,
        ha="left",
        va="center",
        fontsize=9,
        color="#222222",
        bbox=annotation_box,
        zorder=8,
    )

    colorbar = fig.colorbar(image, cax=cax, orientation="vertical")
    tick_candidates = np.asarray([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    ticks = tick_candidates[(tick_candidates >= vmin) & (tick_candidates <= vmax)]
    if ticks.size >= 2:
        colorbar.locator = FixedLocator(ticks)
        colorbar.formatter = FuncFormatter(_format_log_tick)
        colorbar.update_ticks()
    colorbar.set_label(r"Velocity magnitude ($\mathrm{mm\ s^{-1}}$)", fontsize=12)
    colorbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colorbar.outline.set_linewidth(0.7)

    for label, axis in (("a", ax_network), ("b", ax_zoom)):
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
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(direction="out", length=3.5, width=0.7)
        axis.grid(False)

    for y_value in (zoom_bounds[2], zoom_bounds[3]):
        fig.add_artist(
            ConnectionPatch(
                xyA=(zoom_bounds[1], y_value),
                coordsA=ax_network.transData,
                xyB=(zoom_bounds[0], y_value),
                coordsB=ax_zoom.transData,
                color="#C23B4A",
                linewidth=0.65,
                alpha=0.65,
                zorder=0,
                clip_on=False,
            )
        )

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="This figure includes Axes that are not compatible with tight_layout",
            category=UserWarning,
        )
        fig.tight_layout(pad=0.55)
    png_path = output_dir / f"{basename}.png"
    pdf_path = output_dir / f"{basename}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return png_path, pdf_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw a publication-quality vascular CFD-grid construction figure."
    )
    parser.add_argument(
        "--field",
        type=Path,
        default=None,
        help=(
            "velocity_and_wall_shear_field.npz to visualize. By default, the "
            "newest archive in ulm_microbubble_traj_gen/results is used."
        ),
    )
    parser.add_argument(
        "--vessels",
        type=Path,
        default=DEFAULT_VESSEL_ARCHIVE,
        help="Exported .vessels.npz archive used to define centreline topology.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the 300-DPI PNG and PDF files.",
    )
    parser.add_argument(
        "--basename",
        default="vascular_grid_construction",
        help="Output file stem.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    field_archive = (
        args.field.resolve()
        if args.field is not None
        else _latest_field_archive(PACKAGE_DIR / "results")
    )
    png_path, pdf_path = draw_figure(
        field_archive=field_archive,
        vessel_archive=args.vessels.resolve(),
        output_dir=args.output_dir.resolve(),
        basename=str(args.basename),
    )
    print(f"Field archive: {field_archive}")
    print(f"Saved PNG: {png_path}")
    print(f"Saved PDF: {pdf_path}")


if __name__ == "__main__":
    main()
