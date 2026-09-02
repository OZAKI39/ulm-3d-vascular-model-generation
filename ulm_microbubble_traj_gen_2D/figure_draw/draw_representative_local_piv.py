"""Plot representative local PIV views of the converged vascular flow field.

The background encodes velocity magnitude.  Equal-length quiver vectors encode
direction using adaptive downsampling: smooth-flow regions remain sparse while
regions with stronger velocity or directional gradients retain more vectors.
A global context panel identifies all regions of interest.
All data are read from the accepted field archive produced by
``generate_microbubble_trajectories.py``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator

from draw_vascular_grid_construction import (
    DEFAULT_OUTPUT_DIR,
    DEFAULT_VESSEL_ARCHIVE,
    PACKAGE_DIR,
    _cell_edges,
    _configure_style,
    _latest_field_archive,
    _segments_in_box,
)


OUTPUT_BASENAME = "representative_local_piv_velocity_direction"


@dataclass(frozen=True)
class PivRegion:
    label: str
    title: str
    centre_xz_um: np.ndarray
    half_width_um: float
    color: str
    vessel_id: int
    generation: int

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        x, z = np.asarray(self.centre_xz_um, dtype=float)
        half = float(self.half_width_um)
        return x - half, x + half, z - half, z + half


def _load_field(field_archive: Path) -> dict[str, np.ndarray | float]:
    required = {
        "x_coordinates_um",
        "z_coordinates_um",
        "spacing_um",
        "lumen_mask",
        "velocity_xz_um_s",
        "speed_um_s",
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
            "velocity_um_s": np.asarray(archive["velocity_xz_um_s"], dtype=float),
            "speed_um_s": np.asarray(archive["speed_um_s"], dtype=float),
            "wall_start": np.asarray(
                archive["continuous_wall_start_xz_um"], dtype=float
            ),
            "wall_end": np.asarray(
                archive["continuous_wall_end_xz_um"], dtype=float
            ),
        }


def _load_vessels(vessel_archive: Path) -> dict[str, np.ndarray]:
    required = {
        "vessel_id",
        "parent_id",
        "children_count",
        "x_p",
        "x_d",
        "radius_um",
    }
    with np.load(vessel_archive, allow_pickle=False) as archive:
        missing = sorted(required.difference(archive.files))
        if missing:
            raise KeyError(f"Vessel archive is missing required arrays: {missing}")
        return {key: np.asarray(archive[key]) for key in required}


def _vessel_generations(vessel: dict[str, np.ndarray]) -> np.ndarray:
    vessel_ids = np.asarray(vessel["vessel_id"], dtype=int)
    parent_ids = np.asarray(vessel["parent_id"], dtype=int)
    id_to_row = {int(vessel_id): index for index, vessel_id in enumerate(vessel_ids)}
    generations = np.zeros(vessel_ids.size, dtype=int)
    for row in range(vessel_ids.size):
        parent = int(parent_ids[row])
        visited: set[int] = set()
        generation = 0
        while parent >= 0:
            if parent in visited or parent not in id_to_row:
                raise ValueError("The vessel parent graph is cyclic or incomplete.")
            visited.add(parent)
            generation += 1
            parent = int(parent_ids[id_to_row[parent]])
        generations[row] = generation
    return generations


def _select_regions(vessel: dict[str, np.ndarray]) -> tuple[PivRegion, ...]:
    """Choose proximal, central, intermediate, and distal bifurcations."""

    vessel_ids = np.asarray(vessel["vessel_id"], dtype=int)
    child_count = np.asarray(vessel["children_count"], dtype=int)
    endpoints = np.asarray(vessel["x_d"], dtype=float)[:, [0, 2]]
    radii = np.asarray(vessel["radius_um"], dtype=float)
    generations = _vessel_generations(vessel)
    branch_rows = np.flatnonzero(child_count >= 2)
    if branch_rows.size < 4:
        raise ValueError("At least four bifurcations are required for the PIV panel.")

    # The choices are based on topology and physical radius, not on the plotted
    # velocity outcome.  This avoids selecting regions because they happen to
    # look visually favourable in the converged field.
    proximal = int(branch_rows[np.argmin(generations[branch_rows])])

    network_centre = np.mean(
        np.vstack(
            (
                np.asarray(vessel["x_p"], dtype=float)[:, [0, 2]],
                endpoints,
            )
        ),
        axis=0,
    )
    remaining = branch_rows[branch_rows != proximal]
    extent = np.ptp(endpoints, axis=0)
    centrality = 1.0 - np.clip(
        np.linalg.norm((endpoints[remaining] - network_centre) / extent, axis=1),
        0.0,
        1.0,
    )
    radius_score = radii[remaining] / float(radii[remaining].max())
    central = int(remaining[np.argmax(0.65 * radius_score + 0.35 * centrality)])

    # Intermediate: the largest generation-7-to-9 branch separated from the
    # proximal and central views.
    intermediate_candidates = branch_rows[
        (generations[branch_rows] >= 7) & (generations[branch_rows] <= 9)
    ]
    if intermediate_candidates.size == 0:
        intermediate_candidates = remaining
    intermediate = int(
        intermediate_candidates[np.argmax(radii[intermediate_candidates])]
    )

    # Distal: deepest available bifurcation; break ties by choosing the one
    # furthest from the intermediate view so the panels sample distinct beds.
    deepest_generation = int(generations[branch_rows].max())
    distal_candidates = branch_rows[generations[branch_rows] == deepest_generation]
    separation = np.linalg.norm(
        endpoints[distal_candidates] - endpoints[intermediate], axis=1
    )
    distal = int(distal_candidates[np.argmax(separation)])

    rows = (proximal, central, intermediate, distal)
    titles = (
        "Proximal\nbifurcation",
        "Major central\nbifurcation",
        "Intermediate\nbranch bed",
        "Distal\nbifurcation",
    )
    half_widths = (
        max(105.0, 4.2 * float(radii[proximal])),
        max(88.0, 4.0 * float(radii[central])),
        max(68.0, 4.2 * float(radii[intermediate])),
        max(64.0, 4.6 * float(radii[distal])),
    )
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
    return tuple(
        PivRegion(
            label=f"R{index + 1}",
            title=titles[index],
            centre_xz_um=np.asarray(endpoints[row], dtype=float),
            half_width_um=float(half_widths[index]),
            color=colors[index],
            vessel_id=int(vessel_ids[row]),
            generation=int(generations[row]),
        )
        for index, row in enumerate(rows)
    )


def _clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.0, width=0.7)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)


def _local_indices(
    x_um: np.ndarray,
    z_um: np.ndarray,
    bounds: tuple[float, float, float, float],
    spacing_um: float,
) -> tuple[int, int, int, int]:
    xmin, xmax, zmin, zmax = bounds
    ix = np.flatnonzero(
        (x_um >= xmin - spacing_um) & (x_um <= xmax + spacing_um)
    )
    iz = np.flatnonzero(
        (z_um >= zmin - spacing_um) & (z_um <= zmax + spacing_um)
    )
    if ix.size < 3 or iz.size < 3:
        raise ValueError(f"PIV region {bounds} falls outside the field grid.")
    return int(ix[0]), int(ix[-1]) + 1, int(iz[0]), int(iz[-1]) + 1


def _draw_local_piv(
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
    velocity = np.asarray(field["velocity_um_s"], dtype=float)
    speed_um_s = np.asarray(field["speed_um_s"], dtype=float)
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)
    ix0, ix1, iz0, iz1 = _local_indices(
        x_um, z_um, region.bounds, spacing_um
    )

    local_lumen = lumen[ix0:ix1, iz0:iz1]
    local_speed_mm_s = speed_um_s[ix0:ix1, iz0:iz1] / 1000.0
    masked_speed = np.ma.masked_where(~local_lumen, local_speed_mm_s)
    ax.set_facecolor("white")
    ax.pcolormesh(
        x_edges[ix0 : ix1 + 1],
        z_edges[iz0 : iz1 + 1],
        masked_speed.T,
        cmap=cmap,
        norm=norm,
        shading="flat",
        edgecolors="none",
        alpha=0.78,
        rasterized=True,
        zorder=1,
    )

    local_walls = _segments_in_box(wall_start, wall_end, region.bounds)
    ax.add_collection(
        LineCollection(
            local_walls,
            colors="#171717",
            linewidths=0.75,
            zorder=4,
        )
    )

    local_velocity = velocity[ix0:ix1, iz0:iz1]
    valid = local_lumen & np.isfinite(local_speed_mm_s) & (
        local_speed_mm_s > 0.0
    )
    direction_x = np.zeros_like(local_speed_mm_s, dtype=float)
    direction_z = np.zeros_like(local_speed_mm_s, dtype=float)
    direction_x[valid] = (
        local_velocity[..., 0][valid] / speed_um_s[ix0:ix1, iz0:iz1][valid]
    )
    direction_z[valid] = (
        local_velocity[..., 1][valid] / speed_um_s[ix0:ix1, iz0:iz1][valid]
    )

    # Use only cells whose four neighbours also lie inside the lumen when
    # estimating gradients.  This avoids mistaking the binary vessel boundary
    # for a physical high-gradient flow feature.
    interior = valid.copy()
    interior[1:-1, 1:-1] &= (
        valid[:-2, 1:-1]
        & valid[2:, 1:-1]
        & valid[1:-1, :-2]
        & valid[1:-1, 2:]
    )
    interior[[0, -1], :] = False
    interior[:, [0, -1]] = False

    ddir_x_dx, ddir_x_dz = np.gradient(direction_x, spacing_um, spacing_um)
    ddir_z_dx, ddir_z_dz = np.gradient(direction_z, spacing_um, spacing_um)
    turning_strength = np.sqrt(
        ddir_x_dx**2
        + ddir_x_dz**2
        + ddir_z_dx**2
        + ddir_z_dz**2
    )
    log_speed = np.log(np.maximum(local_speed_mm_s, 1.0e-12))
    dspeed_dx, dspeed_dz = np.gradient(log_speed, spacing_um, spacing_um)
    speed_gradient = np.hypot(dspeed_dx, dspeed_dz)

    def robust_normalize(values: np.ndarray) -> np.ndarray:
        samples = values[interior & np.isfinite(values)]
        if samples.size == 0:
            return np.zeros_like(values)
        scale = float(np.percentile(samples, 90.0))
        if not np.isfinite(scale) or scale <= 0.0:
            return np.zeros_like(values)
        return np.clip(values / scale, 0.0, 2.0)

    feature_score = (
        0.65 * robust_normalize(turning_strength)
        + 0.35 * robust_normalize(speed_gradient)
    )
    feature_values = feature_score[interior]
    feature_threshold = (
        float(np.percentile(feature_values, 72.0))
        if feature_values.size
        else np.inf
    )

    cells_across = (2.0 * region.half_width_um) / spacing_um
    coarse_stride = max(8, int(round(cells_across / 12.0)))
    fine_stride = max(4, int(round(coarse_stride / 2.0)))
    grid_i, grid_j = np.indices(local_lumen.shape)
    coarse_grid = (
        (grid_i - coarse_stride // 2) % coarse_stride == 0
    ) & ((grid_j - coarse_stride // 2) % coarse_stride == 0)
    fine_grid = (
        (grid_i - fine_stride // 2) % fine_stride == 0
    ) & ((grid_j - fine_stride // 2) % fine_stride == 0)
    sample_mask = valid & (
        coarse_grid
        | (
            fine_grid
            & interior
            & (feature_score >= feature_threshold)
        )
    )

    arrow_length_um = 0.055 * (2.0 * region.half_width_um)
    local_x, local_z = np.meshgrid(
        x_um[ix0:ix1],
        z_um[iz0:iz1],
        indexing="ij",
    )
    ax.quiver(
        local_x[sample_mask],
        local_z[sample_mask],
        arrow_length_um * direction_x[sample_mask],
        arrow_length_um * direction_z[sample_mask],
        angles="xy",
        scale_units="xy",
        scale=1.0,
        pivot="middle",
        color="#000000",
        width=0.0026,
        headwidth=3.8,
        headlength=5.2,
        headaxislength=4.6,
        minlength=0.0,
        zorder=6,
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
    _clean_axes(ax)
    ax.text(
        -0.17,
        1.045,
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
        fontsize=10,
        fontweight="bold",
        color=region.color,
        bbox=dict(
            boxstyle="round,pad=0.22",
            facecolor="white",
            edgecolor=region.color,
            linewidth=0.9,
            alpha=0.94,
        ),
        zorder=8,
    )
def _format_log_tick(value: float, _: int) -> str:
    return f"{value:g}"


def draw_representative_piv(
    field_archive: Path,
    vessel_archive: Path,
    output_dir: Path,
) -> tuple[Path, Path, tuple[PivRegion, ...]]:
    _configure_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    field = _load_field(field_archive)
    vessel = _load_vessels(vessel_archive)
    regions = _select_regions(vessel)

    x_um = np.asarray(field["x_um"], dtype=float)
    z_um = np.asarray(field["z_um"], dtype=float)
    spacing_um = float(field["spacing_um"])
    lumen = np.asarray(field["lumen"], dtype=bool)
    speed_mm_s = np.asarray(field["speed_um_s"], dtype=float) / 1000.0
    wall_start = np.asarray(field["wall_start"], dtype=float)
    wall_end = np.asarray(field["wall_end"], dtype=float)
    x_edges = _cell_edges(x_um, spacing_um)
    z_edges = _cell_edges(z_um, spacing_um)
    positive_speed = speed_mm_s[lumen & (speed_mm_s > 0.0)]
    vmin = max(0.05, float(np.percentile(positive_speed, 0.5)))
    vmax = float(np.ceil(np.percentile(positive_speed, 99.8) * 2.0) / 2.0)
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=True)
    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad("white")
    masked_speed = np.ma.masked_where(~lumen, speed_mm_s)

    fig = plt.figure(figsize=(13.8, 5.7))
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
        masked_speed.T,
        origin="lower",
        extent=(x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]),
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="equal",
        alpha=0.78,
        rasterized=True,
        zorder=1,
    )
    ax_overview.add_collection(
        LineCollection(
            np.stack((wall_start, wall_end), axis=1),
            colors="#222222",
            linewidths=0.16,
            alpha=0.78,
            rasterized=True,
            zorder=2,
        )
    )
    for region in regions:
        bounds = region.bounds
        ax_overview.add_patch(
            Rectangle(
                (bounds[0], bounds[2]),
                bounds[1] - bounds[0],
                bounds[3] - bounds[2],
                fill=False,
                edgecolor=region.color,
                linewidth=1.15,
                zorder=5,
            )
        )
        ax_overview.text(
            bounds[0],
            bounds[3] + 20.0,
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
        title="PIV region locations",
    )
    ax_overview.xaxis.set_major_locator(MaxNLocator(4))
    ax_overview.yaxis.set_major_locator(MaxNLocator(5))
    _clean_axes(ax_overview)
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

    for index, (axis, region) in enumerate(zip(local_axes, regions)):
        _draw_local_piv(
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

    colorbar = fig.colorbar(overview, cax=cax, orientation="vertical")
    tick_candidates = np.asarray([0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    ticks = tick_candidates[(tick_candidates >= vmin) & (tick_candidates <= vmax)]
    if ticks.size >= 2:
        colorbar.locator = FixedLocator(ticks)
        colorbar.formatter = FuncFormatter(_format_log_tick)
        colorbar.update_ticks()
    colorbar.set_label(r"Velocity magnitude ($\mathrm{mm\ s^{-1}}$)", fontsize=12)
    colorbar.ax.tick_params(labelsize=10, length=3.0, width=0.7)
    colorbar.outline.set_linewidth(0.7)

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
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return png_path, pdf_path, regions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Draw representative local PIV maps of the converged velocity field."
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
    vessel_archive = args.vessels.resolve()
    output_dir = args.output_dir.resolve()
    png_path, pdf_path, regions = draw_representative_piv(
        field_archive,
        vessel_archive,
        output_dir,
    )
    print(f"Field archive: {field_archive}")
    for region in regions:
        x, z = region.centre_xz_um
        print(
            f"{region.label}: {region.title.replace(chr(10), ' ')}; "
            f"vessel={region.vessel_id}; "
            f"generation={region.generation}; centre=({x:.3f}, {z:.3f}) um; "
            f"half-width={region.half_width_um:.3f} um"
        )
    print(f"Saved PNG: {png_path}")
    print(f"Saved PDF: {pdf_path}")


if __name__ == "__main__":
    main()
