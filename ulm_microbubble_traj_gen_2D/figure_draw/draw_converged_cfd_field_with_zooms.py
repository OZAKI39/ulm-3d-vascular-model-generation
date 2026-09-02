"""Plot the converged vascular CFD field with one representative enlargement.

The full-network panel and all zoom panels use the same logarithmic velocity
scale.  Local panels add sparse streamlines to reveal in-plane flow direction
without obscuring the underlying velocity-magnitude field.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator

from draw_representative_local_piv import (
    PivRegion,
    _load_field,
    _load_vessels,
    _local_indices,
    _select_regions,
)
from draw_vascular_grid_construction import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_VESSEL_ARCHIVE,
    PACKAGE_DIR,
    _cell_edges,
    _configure_style,
    _latest_field_archive,
    _segments_in_box,
)


OUTPUT_BASENAME = "converged_cfd_vascular_flow_with_local_zoom"


def _configure_cfd_style() -> None:
    _configure_style()
    mpl.rcParams.update(
        {
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "legend.fontsize": 9,
        }
    )


def _load_solver_status(field_archive: Path) -> dict[str, float | int | bool]:
    with np.load(field_archive, allow_pickle=False) as archive:
        required = {
            "solver_iterations",
            "solver_physical_converged",
            "solver_normalized_momentum_residual",
            "solver_final_normalized_divergence_error",
        }
        missing = sorted(required.difference(archive.files))
        if missing:
            raise KeyError(
                f"Field archive is missing convergence metadata: {missing}"
            )
        return {
            "iterations": int(
                np.asarray(archive["solver_iterations"]).reshape(-1)[0]
            ),
            "converged": bool(
                np.asarray(archive["solver_physical_converged"]).reshape(-1)[0]
            ),
            "momentum_residual": float(
                np.asarray(
                    archive["solver_normalized_momentum_residual"]
                ).reshape(-1)[0]
            ),
            "divergence_error": float(
                np.asarray(
                    archive["solver_final_normalized_divergence_error"]
                ).reshape(-1)[0]
            ),
        }


def _clean_spatial_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.0, width=0.7)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)


def _draw_local_field(
    ax: plt.Axes,
    region: PivRegion,
    field: dict[str, np.ndarray | float],
    cmap: mpl.colors.Colormap,
    norm: LogNorm,
    panel_label: str,
) -> None:
    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen"], dtype=bool)
    wall = np.asarray(field["wall"], dtype=bool)
    fluid = lumen & ~wall
    display_mask = fluid | wall
    velocity = np.asarray(field["velocity_um_s"], dtype=float)
    speed_mm_s = np.asarray(field["speed_um_s"], dtype=float) / 1000.0
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)

    ix0, ix1, iz0, iz1 = _local_indices(
        x_um,
        z_um,
        region.bounds,
        spacing_um,
    )
    local_display = display_mask[ix0:ix1, iz0:iz1]
    local_speed = np.ma.masked_where(
        ~local_display,
        speed_mm_s[ix0:ix1, iz0:iz1],
    )
    ax.pcolormesh(
        x_edges[ix0 : ix1 + 1],
        z_edges[iz0 : iz1 + 1],
        local_speed.T,
        cmap=cmap,
        norm=norm,
        shading="flat",
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )

    local_walls = _segments_in_box(wall_start, wall_end, region.bounds)
    ax.add_collection(
        LineCollection(
            local_walls,
            colors="#191919",
            linewidths=0.72,
            zorder=5,
        )
    )

    local_velocity = velocity[ix0:ix1, iz0:iz1]
    u = np.ma.masked_where(
        ~local_display.T,
        local_velocity[..., 0].T,
    )
    w = np.ma.masked_where(
        ~local_display.T,
        local_velocity[..., 1].T,
    )
    stream = ax.streamplot(
        x_um[ix0:ix1],
        z_um[iz0:iz1],
        u,
        w,
        density=0.88,
        color="white",
        linewidth=0.62,
        arrowsize=0.72,
        minlength=0.20,
        maxlength=8.0,
        integration_direction="both",
        broken_streamlines=False,
        zorder=4,
    )
    outline = [
        path_effects.Stroke(linewidth=1.45, foreground="#111111", alpha=0.78),
        path_effects.Normal(),
    ]
    stream.lines.set_path_effects(outline)
    stream.arrows.set_path_effects(outline)

    ax.set(
        xlim=(region.bounds[0], region.bounds[1]),
        ylim=(region.bounds[2], region.bounds[3]),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title=region.title.replace("\n", " "),
    )
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    _clean_spatial_axes(ax)
    ax.text(
        -0.16,
        1.035,
        f"({panel_label})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        clip_on=False,
    )
    ax.text(
        0.025,
        0.97,
        region.label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=region.color,
        bbox={
            "boxstyle": "round,pad=0.20",
            "facecolor": "white",
            "edgecolor": region.color,
            "linewidth": 0.85,
            "alpha": 0.94,
        },
        zorder=8,
    )


def _format_log_tick(value: float, _: int) -> str:
    return f"{value:g}"


def draw_converged_cfd_field(
    field_archive: Path,
    vessel_archive: Path,
    output_dir: Path,
) -> tuple[Path, Path, tuple[PivRegion, ...]]:
    _configure_cfd_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    field = _load_field(field_archive)
    with np.load(field_archive, allow_pickle=False) as archive:
        if "wall_mask" not in archive.files:
            raise KeyError("Field archive is missing required array: wall_mask")
        field["wall"] = np.asarray(archive["wall_mask"], dtype=bool)
    status = _load_solver_status(field_archive)
    if not bool(status["converged"]):
        raise ValueError(
            f"{field_archive} is not certified as physically converged."
        )
    vessel = _load_vessels(vessel_archive)
    candidate_regions = _select_regions(vessel)
    # The major central bifurcation is the most informative single enlargement:
    # it contains one parent and two resolved daughters while remaining large
    # enough for streamlines and the velocity gradient to be clearly visible.
    selected = candidate_regions[1]
    region = PivRegion(
        label="ROI",
        title="Major central bifurcation",
        centre_xz_um=selected.centre_xz_um,
        half_width_um=selected.half_width_um,
        color="#D55E00",
        vessel_id=selected.vessel_id,
        generation=selected.generation,
    )

    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen"], dtype=bool)
    wall = np.asarray(field["wall"], dtype=bool)
    fluid = lumen & ~wall
    display_mask = fluid | wall
    speed_mm_s = np.asarray(field["speed_um_s"], dtype=float) / 1000.0
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)

    positive_speed = speed_mm_s[
        display_mask & np.isfinite(speed_mm_s) & (speed_mm_s > 0)
    ]
    if positive_speed.size == 0:
        raise ValueError("The converged field contains no positive lumen speed.")
    vmin = max(0.05, float(np.percentile(positive_speed, 0.5)))
    vmax = float(np.ceil(np.percentile(positive_speed, 99.8) * 2.0) / 2.0)
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=True)
    # Viridis provides a perceptually uniform and visibly continuous speed
    # gradient from slow distal flow to fast proximal flow.
    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad("white")
    masked_speed = np.ma.masked_where(~display_mask, speed_mm_s)

    fig = plt.figure(figsize=(10.2, 5.45))
    spec = fig.add_gridspec(
        1,
        3,
        width_ratios=(1.28, 1.0, 0.055),
        wspace=0.30,
    )
    ax_overview = fig.add_subplot(spec[0, 0])
    ax_local = fig.add_subplot(spec[0, 1])
    cax = fig.add_subplot(spec[0, 2])

    overview = ax_overview.imshow(
        masked_speed.T,
        origin="lower",
        extent=(x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]),
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
        rasterized=True,
        zorder=1,
    )
    ax_overview.add_collection(
        LineCollection(
            np.stack((wall_start, wall_end), axis=1),
            colors="#222222",
            linewidths=0.10,
            alpha=0.45,
            rasterized=True,
            zorder=2,
        )
    )
    xmin, xmax, zmin, zmax = region.bounds
    ax_overview.add_patch(
        Rectangle(
            (xmin, zmin),
            xmax - xmin,
            zmax - zmin,
            fill=False,
            edgecolor=region.color,
            linewidth=1.15,
            zorder=5,
        )
    )
    ax_overview.text(
        xmin,
        zmax + 18.0,
        region.label,
        color=region.color,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="bottom",
        zorder=6,
    )

    ax_overview.set(
        xlim=(float(x_edges[0]), float(x_edges[-1])),
        ylim=(float(z_edges[0]), float(z_edges[-1])),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Converged CFD velocity field",
    )
    ax_overview.xaxis.set_major_locator(MaxNLocator(4))
    ax_overview.yaxis.set_major_locator(MaxNLocator(5))
    _clean_spatial_axes(ax_overview)
    ax_overview.text(
        -0.13,
        1.025,
        "(a)",
        transform=ax_overview.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        clip_on=False,
    )
    ax_overview.text(
        0.035,
        0.975,
        (
            f"Converged at iteration {int(status['iterations'])}\n"
            f"Momentum residual = {float(status['momentum_residual']):.2e}"
        ),
        transform=ax_overview.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#333333",
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "edgecolor": "#777777",
            "linewidth": 0.65,
            "alpha": 0.92,
        },
        zorder=8,
    )

    _draw_local_field(
        ax_local,
        region,
        field,
        cmap,
        norm,
        panel_label="b",
    )

    colorbar = fig.colorbar(overview, cax=cax, orientation="vertical")
    tick_candidates = np.asarray(
        [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    )
    ticks = tick_candidates[(tick_candidates >= vmin) & (tick_candidates <= vmax)]
    if ticks.size >= 2:
        colorbar.locator = FixedLocator(ticks)
        colorbar.formatter = FuncFormatter(_format_log_tick)
        colorbar.update_ticks()
    colorbar.set_label(
        r"Velocity magnitude ($\mathrm{mm\ s^{-1}}$)",
        fontsize=12,
    )
    colorbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colorbar.outline.set_linewidth(0.7)

    fig.text(
        0.69,
        0.010,
        "Streamlines indicate in-plane flow direction.",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#333333",
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="This figure includes Axes that are not compatible with tight_layout",
            category=UserWarning,
        )
        fig.tight_layout(pad=0.55, rect=(0.0, 0.045, 1.0, 1.0))

    # Match the colorbar to the physical plotting height of the full-network
    # panel; the dedicated GridSpec column would otherwise span the full figure.
    fig.canvas.draw()
    overview_box = ax_overview.get_position()
    colorbar_box = cax.get_position()
    cax.set_position(
        [
            colorbar_box.x0,
            overview_box.y0,
            colorbar_box.width,
            overview_box.height,
        ]
    )

    png_path = output_dir / f"{OUTPUT_BASENAME}.png"
    pdf_path = output_dir / f"{OUTPUT_BASENAME}.pdf"
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
    return png_path, pdf_path, (region,)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw the converged CFD vascular velocity field with local zooms."
        )
    )
    parser.add_argument(
        "--field",
        type=Path,
        default=None,
        help="Field archive; defaults to the newest complete result.",
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
        help="Directory for the 300-DPI PNG and vector PDF.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    field_archive = (
        args.field.resolve()
        if args.field is not None
        else _latest_field_archive(PACKAGE_DIR / "results")
    )
    png_path, pdf_path, regions = draw_converged_cfd_field(
        field_archive,
        args.vessels.resolve(),
        args.output_dir.resolve(),
    )
    print(f"Field archive: {field_archive}")
    print(
        "Zoom regions: "
        + ", ".join(
            f"{region.label}=vessel {region.vessel_id} (generation {region.generation})"
            for region in regions
        )
    )
    print(f"Saved 300-DPI PNG: {png_path}")
    print(f"Saved vector PDF: {pdf_path}")


if __name__ == "__main__":
    main()
