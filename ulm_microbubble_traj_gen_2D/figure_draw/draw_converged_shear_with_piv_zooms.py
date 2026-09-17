"""Plot the converged vascular shear field at the four PIV zoom locations.

The stored scalar is the solver's two-dimensional in-plane wall-shear-stress
proxy in pascals.  Values are not numerically smoothed; Gouraud interpolation
is used only for rendering the local cell-centred field without blocky pixels.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import PowerNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MaxNLocator
from scipy.ndimage import distance_transform_edt

from draw_converged_cfd_field_with_zooms import (
    _clean_spatial_axes,
    _configure_cfd_style,
    _load_solver_status,
)
from draw_representative_local_piv import (
    PivRegion,
    _load_vessels,
    _local_indices,
    _select_regions,
)
from draw_vascular_grid_construction import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_VESSEL_ARCHIVE,
    PACKAGE_DIR,
    _cell_edges,
    _latest_field_archive,
    _segments_in_box,
)


OUTPUT_BASENAME = "converged_shear_field_with_piv_matched_zooms"


def _load_shear_field(
    field_archive: Path,
) -> dict[str, np.ndarray | float]:
    required = {
        "x_coordinates_um",
        "z_coordinates_um",
        "spacing_um",
        "lumen_mask",
        "wall_shear_stress_pa",
        "continuous_wall_start_xz_um",
        "continuous_wall_end_xz_um",
    }
    with np.load(field_archive, allow_pickle=False) as archive:
        missing = sorted(required.difference(archive.files))
        if missing:
            raise KeyError(f"Field archive is missing required arrays: {missing}")
        return {
            "x_um": np.asarray(archive["x_coordinates_um"], dtype=float),
            "z_um": np.asarray(archive["z_coordinates_um"], dtype=float),
            "spacing_um": float(np.asarray(archive["spacing_um"]).reshape(-1)[0]),
            "lumen": np.asarray(archive["lumen_mask"], dtype=bool),
            "shear_pa": np.asarray(
                archive["wall_shear_stress_pa"], dtype=float
            ),
            "wall_start": np.asarray(
                archive["continuous_wall_start_xz_um"], dtype=float
            ),
            "wall_end": np.asarray(
                archive["continuous_wall_end_xz_um"], dtype=float
            ),
        }


def _panel_label(ax: plt.Axes, label: str, x: float = -0.17) -> None:
    ax.text(
        x,
        1.035,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        clip_on=False,
    )


def _draw_local_shear(
    ax: plt.Axes,
    region: PivRegion,
    field: dict[str, np.ndarray | float],
    cmap: mpl.colors.Colormap,
    norm: PowerNorm,
    panel_label: str,
) -> mpl.image.AxesImage:
    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen"], dtype=bool)
    shear_pa = np.asarray(field["shear_pa"], dtype=float)
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
    local_lumen = lumen[ix0:ix1, iz0:iz1]
    local_shear = shear_pa[ix0:ix1, iz0:iz1]
    # Extend the nearest valid lumen value into transparent pixels, then
    # interpolate the RGBA image.  The alpha channel yields a clean antialiased
    # lumen boundary without introducing zero-valued dark halos.  This is a
    # display operation only; the archived scalar values remain unchanged.
    nearest = distance_transform_edt(
        ~local_lumen,
        return_distances=False,
        return_indices=True,
    )
    filled_shear = local_shear[tuple(nearest)]
    rgba = cmap(norm(filled_shear))
    rgba[..., 3] = local_lumen.astype(float)
    image = ax.imshow(
        np.swapaxes(rgba, 0, 1),
        origin="lower",
        extent=(
            x_edges[ix0],
            x_edges[ix1],
            z_edges[iz0],
            z_edges[iz1],
        ),
        interpolation="bicubic",
        interpolation_stage="rgba",
        aspect="equal",
        rasterized=True,
        zorder=1,
    )
    local_walls = _segments_in_box(wall_start, wall_end, region.bounds)
    ax.add_collection(
        LineCollection(
            local_walls,
            colors="#181818",
            linewidths=0.68,
            alpha=0.88,
            zorder=4,
        )
    )
    ax.set(
        xlim=(region.bounds[0], region.bounds[1]),
        ylim=(region.bounds[2], region.bounds[3]),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title=region.title,
    )
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    _clean_spatial_axes(ax)
    _panel_label(ax, panel_label)
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
    return image


def draw_converged_shear_with_piv_zooms(
    field_archive: Path,
    vessel_archive: Path,
    output_dir: Path,
) -> tuple[Path, Path, tuple[PivRegion, ...]]:
    _configure_cfd_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    field = _load_shear_field(field_archive)
    status = _load_solver_status(field_archive)
    if not bool(status["converged"]):
        raise ValueError(
            f"{field_archive} is not certified as physically converged."
        )
    vessel = _load_vessels(vessel_archive)
    regions = _select_regions(vessel)

    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen"], dtype=bool)
    shear_pa = np.asarray(field["shear_pa"], dtype=float)
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)

    finite_shear = shear_pa[
        lumen & np.isfinite(shear_pa) & (shear_pa >= 0.0)
    ]
    if finite_shear.size == 0:
        raise ValueError("The converged field contains no finite shear values.")
    vmax = float(np.ceil(np.percentile(finite_shear, 99.8) * 10.0) / 10.0)
    if vmax <= 0.0:
        vmax = 1.0
    # A perceptual power law opens the low-to-intermediate range while retaining
    # a physically ordered, zero-based scale and a shared norm for every panel.
    norm = PowerNorm(gamma=0.62, vmin=0.0, vmax=vmax, clip=False)
    cmap = mpl.colormaps["magma"].copy()
    cmap.set_bad("white")
    masked_shear = np.ma.masked_where(~lumen, shear_pa)

    fig = plt.figure(figsize=(13.8, 5.75))
    spec = fig.add_gridspec(
        2,
        4,
        width_ratios=(3.0, 1.0, 1.0, 0.055),
        height_ratios=(1.0, 1.0),
        wspace=0.30,
        hspace=0.45,
    )
    ax_overview = fig.add_subplot(spec[:, 0])
    local_axes = (
        fig.add_subplot(spec[0, 1]),
        fig.add_subplot(spec[0, 2]),
        fig.add_subplot(spec[1, 1]),
        fig.add_subplot(spec[1, 2]),
    )
    cax = fig.add_subplot(spec[:, 3])

    overview = ax_overview.imshow(
        masked_shear.T,
        origin="lower",
        extent=(x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]),
        cmap=cmap,
        norm=norm,
        interpolation="bilinear",
        interpolation_stage="rgba",
        aspect="equal",
        rasterized=True,
        zorder=1,
    )
    ax_overview.add_collection(
        LineCollection(
            np.stack((wall_start, wall_end), axis=1),
            colors="#242424",
            linewidths=0.13,
            alpha=0.62,
            rasterized=True,
            zorder=3,
        )
    )
    for region in regions:
        xmin, xmax, zmin, zmax = region.bounds
        ax_overview.add_patch(
            Rectangle(
                (xmin, zmin),
                xmax - xmin,
                zmax - zmin,
                fill=False,
                edgecolor=region.color,
                linewidth=1.15,
                zorder=6,
            )
        )
        ax_overview.text(
            xmin,
            zmax + 20.0,
            region.label,
            color=region.color,
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
            zorder=7,
        )
    ax_overview.set(
        xlim=(float(x_edges[0]), float(x_edges[-1])),
        ylim=(float(z_edges[0]), float(z_edges[-1])),
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title="Converged shear-stress field",
    )
    ax_overview.title.set_fontsize(12)
    ax_overview.xaxis.set_major_locator(MaxNLocator(4))
    ax_overview.yaxis.set_major_locator(MaxNLocator(5))
    _clean_spatial_axes(ax_overview)
    _panel_label(ax_overview, "a", x=-0.13)
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
        fontsize=8.2,
        color="#333333",
        bbox={
            "boxstyle": "round,pad=0.24",
            "facecolor": "white",
            "edgecolor": "#777777",
            "linewidth": 0.62,
            "alpha": 0.92,
        },
        zorder=9,
    )

    for index, (axis, region) in enumerate(zip(local_axes, regions)):
        _draw_local_shear(
            axis,
            region,
            field,
            cmap,
            norm,
            panel_label=chr(ord("b") + index),
        )
        if index in (0, 1):
            axis.set_xlabel("")
        if index in (1, 3):
            axis.set_ylabel("")

    colourbar = fig.colorbar(
        overview,
        cax=cax,
        orientation="vertical",
        extend="max",
        extendfrac=0.025,
    )
    tick_step = 0.2 if vmax <= 1.6 else 0.5
    ticks = np.arange(0.0, vmax + 0.5 * tick_step, tick_step)
    colourbar.locator = FixedLocator(ticks)
    colourbar.formatter = FormatStrFormatter("%.1f")
    colourbar.update_ticks()
    colourbar.set_label(
        r"In-plane shear-stress proxy ($\mathrm{Pa}$)",
        fontsize=12,
    )
    colourbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colourbar.outline.set_linewidth(0.7)

    # Make the global overview and colour scale span the full height occupied
    # by the two rows of local panels while retaining physical equal aspect.
    fig.canvas.draw()
    local_bottom = min(axis.get_position().y0 for axis in local_axes)
    local_top = max(axis.get_position().y1 for axis in local_axes)
    target_height = local_top - local_bottom
    figure_width, figure_height = fig.get_size_inches()
    data_aspect = (x_edges[-1] - x_edges[0]) / (z_edges[-1] - z_edges[0])
    target_width = target_height * (figure_height / figure_width) * data_aspect
    overview_box = ax_overview.get_position()
    ax_overview.set_position(
        [
            overview_box.x0,
            local_bottom,
            target_width,
            target_height,
        ]
    )
    colorbar_box = cax.get_position()
    cax.set_position(
        [
            colorbar_box.x0,
            local_bottom,
            colorbar_box.width,
            target_height,
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
    return png_path, pdf_path, regions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw the converged shear-stress field at the four PIV-matched "
            "zoom locations."
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
        help="Directory for the 300-DPI PNG and PDF outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    field_archive = (
        args.field.resolve()
        if args.field is not None
        else _latest_field_archive(PACKAGE_DIR / "results")
    )
    png_path, pdf_path, regions = draw_converged_shear_with_piv_zooms(
        field_archive,
        args.vessels.resolve(),
        args.output_dir.resolve(),
    )
    print(f"Field archive: {field_archive}")
    for region in regions:
        x_um, z_um = region.centre_xz_um
        print(
            f"{region.label}: vessel={region.vessel_id}; "
            f"generation={region.generation}; "
            f"centre=({x_um:.3f}, {z_um:.3f}) um; "
            f"half-width={region.half_width_um:.3f} um"
        )
    print(f"Saved 300-DPI PNG: {png_path}")
    print(f"Saved PDF: {pdf_path}")


if __name__ == "__main__":
    main()
