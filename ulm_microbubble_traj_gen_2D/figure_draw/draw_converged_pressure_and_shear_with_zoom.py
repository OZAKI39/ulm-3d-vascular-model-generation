"""Draw converged projection-pressure and wall-shear fields with one zoom.

The pressure solution stored by the finite-volume projection is a gauge
potential rather than a calibrated thermodynamic pressure, so it is reported
in arbitrary units.  The shear result is the solver's two-dimensional
in-plane wall-shear-stress proxy and is reported in pascals.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

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


OUTPUT_BASENAME = "converged_pressure_and_wall_shear_with_local_zoom"
ROI_COLOUR = "#D55E00"


def _load_scalar_fields(
    field_archive: Path,
) -> dict[str, np.ndarray | float]:
    required = {
        "x_coordinates_um",
        "z_coordinates_um",
        "spacing_um",
        "lumen_mask",
        "pressure",
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
            "pressure_kau": np.asarray(archive["pressure"], dtype=float) / 1000.0,
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


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.145,
        1.025,
        f"({label})",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        clip_on=False,
    )


def _draw_roi(
    ax: plt.Axes,
    region: PivRegion,
    *,
    label_offset_um: float,
) -> None:
    xmin, xmax, zmin, zmax = region.bounds
    ax.add_patch(
        Rectangle(
            (xmin, zmin),
            xmax - xmin,
            zmax - zmin,
            fill=False,
            edgecolor=ROI_COLOUR,
            linewidth=1.15,
            zorder=6,
        )
    )
    ax.text(
        xmin,
        zmax + label_offset_um,
        "ROI",
        color=ROI_COLOUR,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="bottom",
        zorder=7,
    )


def _format_spatial_panel(
    ax: plt.Axes,
    *,
    title: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    panel_label: str,
    local: bool,
) -> None:
    ax.set(
        xlim=xlim,
        ylim=ylim,
        xlabel=r"$x$ ($\mathrm{\mu m}$)",
        ylabel=r"$z$ ($\mathrm{\mu m}$)",
        title=title,
    )
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4 if local else 5))
    _clean_spatial_axes(ax)
    _panel_label(ax, panel_label)


def _draw_global_scalar(
    ax: plt.Axes,
    values: np.ndarray,
    lumen: np.ndarray,
    extent: tuple[float, float, float, float],
    wall_segments: np.ndarray,
    cmap: mpl.colors.Colormap,
    norm: mpl.colors.Normalize,
    region: PivRegion,
    *,
    title: str,
    panel_label: str,
) -> mpl.image.AxesImage:
    masked = np.ma.masked_where(~lumen, values)
    image = ax.imshow(
        masked.T,
        origin="lower",
        extent=extent,
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
        rasterized=True,
        zorder=1,
    )
    ax.add_collection(
        LineCollection(
            wall_segments,
            colors="#202020",
            linewidths=0.10,
            alpha=0.52,
            rasterized=True,
            zorder=3,
        )
    )
    _draw_roi(ax, region, label_offset_um=18.0)
    _format_spatial_panel(
        ax,
        title=title,
        xlim=(extent[0], extent[1]),
        ylim=(extent[2], extent[3]),
        panel_label=panel_label,
        local=False,
    )
    return image


def _draw_local_scalar(
    ax: plt.Axes,
    values: np.ndarray,
    field: dict[str, np.ndarray | float],
    region: PivRegion,
    cmap: mpl.colors.Colormap,
    norm: mpl.colors.Normalize,
    *,
    title: str,
    panel_label: str,
    contour_count: int = 0,
) -> mpl.collections.QuadMesh:
    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen"], dtype=bool)
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
    local_values = np.ma.masked_where(
        ~local_lumen,
        values[ix0:ix1, iz0:iz1],
    )
    mesh = ax.pcolormesh(
        x_edges[ix0 : ix1 + 1],
        z_edges[iz0 : iz1 + 1],
        local_values.T,
        cmap=cmap,
        norm=norm,
        shading="flat",
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )
    if contour_count >= 2:
        finite_local = local_values.compressed()
        if finite_local.size:
            lower, upper = np.percentile(finite_local, (2.0, 98.0))
            if upper > lower:
                ax.contour(
                    x_um[ix0:ix1],
                    z_um[iz0:iz1],
                    local_values.T,
                    levels=np.linspace(lower, upper, contour_count),
                    colors="#4A4A4A",
                    linewidths=0.42,
                    alpha=0.58,
                    zorder=3,
                )
    local_walls = _segments_in_box(wall_start, wall_end, region.bounds)
    ax.add_collection(
        LineCollection(
            local_walls,
            colors="#171717",
            linewidths=0.74,
            zorder=4,
        )
    )
    _format_spatial_panel(
        ax,
        title=title,
        xlim=(region.bounds[0], region.bounds[1]),
        ylim=(region.bounds[2], region.bounds[3]),
        panel_label=panel_label,
        local=True,
    )
    ax.text(
        0.025,
        0.97,
        "ROI",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=ROI_COLOUR,
        bbox={
            "boxstyle": "round,pad=0.20",
            "facecolor": "white",
            "edgecolor": ROI_COLOUR,
            "linewidth": 0.85,
            "alpha": 0.94,
        },
        zorder=8,
    )
    return mesh


def _add_colourbar(
    fig: plt.Figure,
    artist: mpl.cm.ScalarMappable,
    cax: plt.Axes,
    label: str,
) -> None:
    colourbar = fig.colorbar(artist, cax=cax, orientation="vertical")
    colourbar.set_label(label, fontsize=12)
    colourbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colourbar.ax.yaxis.set_major_locator(MaxNLocator(6))
    colourbar.outline.set_linewidth(0.7)


def draw_converged_pressure_and_shear(
    field_archive: Path,
    vessel_archive: Path,
    output_dir: Path,
) -> tuple[Path, Path, PivRegion]:
    _configure_cfd_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    field = _load_scalar_fields(field_archive)
    status = _load_solver_status(field_archive)
    if not bool(status["converged"]):
        raise ValueError(
            f"{field_archive} is not certified as physically converged."
        )

    vessel = _load_vessels(vessel_archive)
    selected = _select_regions(vessel)[1]
    region = PivRegion(
        label="ROI",
        title="Major central bifurcation",
        centre_xz_um=selected.centre_xz_um,
        half_width_um=selected.half_width_um,
        color=ROI_COLOUR,
        vessel_id=selected.vessel_id,
        generation=selected.generation,
    )

    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen"], dtype=bool)
    pressure_kau = np.asarray(field["pressure_kau"], dtype=float)
    shear_pa = np.asarray(field["shear_pa"], dtype=float)
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)
    extent = (
        float(x_edges[0]),
        float(x_edges[-1]),
        float(z_edges[0]),
        float(z_edges[-1]),
    )
    wall_segments = np.stack((wall_start, wall_end), axis=1)

    pressure_values = pressure_kau[
        lumen & np.isfinite(pressure_kau)
    ]
    if pressure_values.size == 0:
        raise ValueError("The converged field contains no finite pressure values.")
    pressure_limit = float(
        np.ceil(np.percentile(np.abs(pressure_values), 99.8) / 5.0) * 5.0
    )
    if pressure_limit <= 0.0:
        pressure_limit = 1.0
    pressure_norm = TwoSlopeNorm(
        vmin=-pressure_limit,
        vcenter=0.0,
        vmax=pressure_limit,
    )
    pressure_cmap = mpl.colormaps["coolwarm"].copy()
    pressure_cmap.set_bad("white")

    shear_values = shear_pa[
        lumen & np.isfinite(shear_pa) & (shear_pa >= 0.0)
    ]
    if shear_values.size == 0:
        raise ValueError("The converged field contains no finite shear values.")
    shear_limit = float(
        np.ceil(np.percentile(shear_values, 99.8) * 10.0) / 10.0
    )
    if shear_limit <= 0.0:
        shear_limit = 1.0
    shear_norm = Normalize(vmin=0.0, vmax=shear_limit, clip=True)
    shear_cmap = mpl.colormaps["viridis"].copy()
    shear_cmap.set_bad("white")

    fig = plt.figure(figsize=(10.3, 9.6))
    spec = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.23, 1.0, 0.050),
        height_ratios=(1.0, 1.0),
        wspace=0.30,
        hspace=0.32,
    )
    ax_pressure_global = fig.add_subplot(spec[0, 0])
    ax_pressure_local = fig.add_subplot(spec[0, 1])
    cax_pressure = fig.add_subplot(spec[0, 2])
    ax_shear_global = fig.add_subplot(spec[1, 0])
    ax_shear_local = fig.add_subplot(spec[1, 1])
    cax_shear = fig.add_subplot(spec[1, 2])

    pressure_artist = _draw_global_scalar(
        ax_pressure_global,
        pressure_kau,
        lumen,
        extent,
        wall_segments,
        pressure_cmap,
        pressure_norm,
        region,
        title="Converged projection-pressure field",
        panel_label="a",
    )
    _draw_local_scalar(
        ax_pressure_local,
        pressure_kau,
        field,
        region,
        pressure_cmap,
        pressure_norm,
        title="Pressure: central bifurcation",
        panel_label="b",
        contour_count=6,
    )
    shear_artist = _draw_global_scalar(
        ax_shear_global,
        shear_pa,
        lumen,
        extent,
        wall_segments,
        shear_cmap,
        shear_norm,
        region,
        title="Converged in-plane wall-shear-stress proxy",
        panel_label="c",
    )
    _draw_local_scalar(
        ax_shear_local,
        shear_pa,
        field,
        region,
        shear_cmap,
        shear_norm,
        title="Wall shear stress: central bifurcation",
        panel_label="d",
    )

    ax_pressure_global.text(
        0.035,
        0.975,
        (
            f"Converged at iteration {int(status['iterations'])}\n"
            f"Momentum residual = {float(status['momentum_residual']):.2e}"
        ),
        transform=ax_pressure_global.transAxes,
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
        zorder=9,
    )

    _add_colourbar(
        fig,
        pressure_artist,
        cax_pressure,
        r"Projection pressure potential ($10^3$ a.u.)",
    )
    _add_colourbar(
        fig,
        shear_artist,
        cax_shear,
        r"Wall shear stress proxy ($\mathrm{Pa}$)",
    )

    fig.text(
        0.54,
        0.012,
        (
            "The same central-bifurcation ROI is used in both rows; "
            "pressure is a zero-centred projection potential."
        ),
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
        fig.tight_layout(pad=0.55, rect=(0.0, 0.035, 1.0, 1.0))

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
    return png_path, pdf_path, region


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw converged projection-pressure and wall-shear fields with "
            "one matched local enlargement."
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
    png_path, pdf_path, region = draw_converged_pressure_and_shear(
        field_archive,
        args.vessels.resolve(),
        args.output_dir.resolve(),
    )
    x_um, z_um = region.centre_xz_um
    print(f"Field archive: {field_archive}")
    print(
        f"ROI: vessel {region.vessel_id}, generation {region.generation}, "
        f"centre=({x_um:.3f}, {z_um:.3f}) um, "
        f"half-width={region.half_width_um:.3f} um"
    )
    print(f"Saved 300-DPI PNG: {png_path}")
    print(f"Saved vector PDF: {pdf_path}")


if __name__ == "__main__":
    main()
