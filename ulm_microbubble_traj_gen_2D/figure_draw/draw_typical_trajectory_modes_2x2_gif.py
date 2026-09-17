"""Create a slow 2x2 GIF of three representative microbubble trajectories.

Layout
------
top-left
    Global continuous-wall context with three moving regions of interest.
top-right
    Normal bifurcation selection.
bottom-left
    Non-targeted near-wall motion.
bottom-right
    Targeted molecular binding and rolling.

Only the authoritative continuous-wall geometry is drawn.  No raster lumen or
wall mask is rendered, so lumen and exterior share one background while vessel
edges remain smooth.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle

try:
    from .draw_molecular_binding_global_local_gif import (
        BOND_COLOR,
        FOCAL_EDGE_COLOR,
        TARGET_COLOR,
        _bond_segments,
        _configure_style,
    )
    from .draw_typical_trajectory_modes_global_local_gif import (
        BINDING_COLOR,
        NORMAL_COLOR,
        OTHER_BUBBLE_COLOR,
        TRAIL_COLOR,
        WALL_COLOR,
        WALL_FOLLOWING_COLOR,
        LocalMode,
        _load_geometry,
        _load_trajectory,
        _record_map,
        _select_bifurcation_mode,
        _select_binding_mode,
        _select_wall_following_mode,
    )
except ImportError:
    from draw_molecular_binding_global_local_gif import (
        BOND_COLOR,
        FOCAL_EDGE_COLOR,
        TARGET_COLOR,
        _bond_segments,
        _configure_style,
    )
    from draw_typical_trajectory_modes_global_local_gif import (
        BINDING_COLOR,
        NORMAL_COLOR,
        OTHER_BUBBLE_COLOR,
        TRAIL_COLOR,
        WALL_COLOR,
        WALL_FOLLOWING_COLOR,
        LocalMode,
        _load_geometry,
        _load_trajectory,
        _record_map,
        _select_bifurcation_mode,
        _select_binding_mode,
        _select_wall_following_mode,
    )


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BINDING_RESULT = SCRIPT_DIR / "outputs" / "molecular_binding_demo_data"
DEFAULT_REFERENCE_RESULT = SCRIPT_DIR / "outputs" / "steady_reference_demo_data"
DEFAULT_OUTPUT = SCRIPT_DIR / "outputs" / "molecular_binding_global_local.gif"

BACKGROUND_COLOR = "#F7F8FA"
BUBBLE_FACE_COLOR = "#2A7AB0"
BUBBLE_EDGE_COLOR = "#173F5F"
BINDING_RING_COLOR = "#E4572E"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a slow 2x2 GIF with global context and three simultaneous "
            "representative local trajectory modes."
        )
    )
    parser.add_argument(
        "--binding-result", type=Path, default=DEFAULT_BINDING_RESULT
    )
    parser.add_argument(
        "--reference-result", type=Path, default=DEFAULT_REFERENCE_RESULT
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--frames",
        type=int,
        default=180,
        help="Animation frames. The default gives 18 s at 10 fps.",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=10,
        help="Playback rate; deliberately slower than the previous 15 fps.",
    )
    parser.add_argument(
        "--trail-duration-s",
        type=float,
        default=0.14,
        help="Physical trajectory history shown behind each focal bubble.",
    )
    parser.add_argument(
        "--preview",
        type=Path,
        default=None,
        help="Optional 2x2 PNG at the maximum-force binding frame.",
    )
    return parser.parse_args()


def _smooth_wall_collection(
    segments: np.ndarray,
    *,
    color: str,
    linewidth: float,
    alpha: float = 1.0,
    zorder: float = 1.0,
) -> LineCollection:
    collection = LineCollection(
        segments,
        colors=color,
        linewidths=linewidth,
        alpha=alpha,
        antialiaseds=True,
        zorder=zorder,
    )
    collection.set_capstyle("round")
    collection.set_joinstyle("round")
    return collection


def _mode_path(
    mode: LocalMode,
    record_map: dict[int, int],
    source_frames: np.ndarray,
) -> np.ndarray:
    lo = int(np.min(source_frames))
    hi = int(np.max(source_frames))
    rows = [
        record_map[frame]
        for frame in range(lo, hi + 1)
        if frame in record_map
    ]
    return mode.data.positions_xz_um[np.asarray(rows, dtype=np.int64)]


def _linear_source_frames(
    mode: LocalMode,
    output_count: int,
    physical_duration_s: float,
) -> np.ndarray:
    """Map equal physical time to every output frame without easing."""

    if output_count < 2:
        raise ValueError("At least two output frames are required.")
    start = int(np.min(mode.frames))
    available_end = int(np.max(mode.frames))
    requested_end = start + int(round(float(physical_duration_s) / mode.data.dt_s))
    end = min(available_end, requested_end)
    return np.rint(np.linspace(start, end, int(output_count))).astype(np.int64)


def _actual_speed_um_s(
    mode: LocalMode,
    record_map: dict[int, int],
    source_frame: int,
    *,
    averaging_duration_s: float = 0.010,
) -> float:
    """Measure centre speed from accepted positions over a short time window."""

    half_window = max(
        1,
        int(round(0.5 * float(averaging_duration_s) / mode.data.dt_s)),
    )
    first_available = int(np.min(mode.frames))
    last_available = int(np.max(mode.frames))
    first = max(first_available, int(source_frame) - half_window)
    last = min(last_available, int(source_frame) + half_window)
    while first < source_frame and first not in record_map:
        first += 1
    while last > source_frame and last not in record_map:
        last -= 1
    if last <= first:
        record = record_map[int(source_frame)]
        return float(
            np.linalg.norm(mode.data.realized_velocities_xz_um_s[record])
        )
    displacement_um = (
        mode.data.positions_xz_um[record_map[last]]
        - mode.data.positions_xz_um[record_map[first]]
    )
    duration_s = (last - first) * mode.data.dt_s
    return float(np.linalg.norm(displacement_um) / duration_s)


def _require_unmodulated_trajectory(result_dir: Path) -> None:
    """Reject a trajectory whose stored cardiac multiplier changes with time."""

    trajectory_path = (
        Path(result_dir).resolve() / "microbubble_field_trajectories.npz"
    )
    with np.load(trajectory_path, allow_pickle=False) as trajectory:
        if "record_cardiac_multiplier" not in trajectory:
            return
        multiplier = np.asarray(
            trajectory["record_cardiac_multiplier"], dtype=float
        )
    if multiplier.size and not np.allclose(
        multiplier, 1.0, rtol=0.0, atol=1.0e-6
    ):
        raise ValueError(
            f"Cardiac pulsatility is active in {Path(result_dir).resolve()}. "
            "Regenerate this trajectory with cardiac_pulsatility.enabled=false."
        )


def _segments_near_path(
    starts: np.ndarray,
    ends: np.ndarray,
    path: np.ndarray,
    padding_um: float,
) -> np.ndarray:
    lower = np.min(path, axis=0) - float(padding_um)
    upper = np.max(path, axis=0) + float(padding_um)
    centres = 0.5 * (starts + ends)
    keep = np.all((centres >= lower[None, :]) & (centres <= upper[None, :]), axis=1)
    return np.stack((starts[keep], ends[keep]), axis=1)


def _orientation_segment(
    axis: plt.Axes,
    centre_xz: np.ndarray,
    angle_rad: float,
    half_length_px: float = 9.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a fixed-screen-length orientation mark in physical coordinates."""

    display_centre = axis.transData.transform(centre_xz)
    direction = np.asarray([np.cos(angle_rad), np.sin(angle_rad)], dtype=float)
    display_start = display_centre - float(half_length_px) * direction
    display_end = display_centre + float(half_length_px) * direction
    inverse = axis.transData.inverted()
    return inverse.transform(display_start), inverse.transform(display_end)


def _update_connection_pair(
    connections: tuple[ConnectionPatch, ConnectionPatch],
    rectangle_bounds: tuple[float, float, float, float],
    mode_index: int,
) -> None:
    x0, x1, z0, z1 = rectangle_bounds
    first, second = connections
    if mode_index == 0:
        first.xy1 = (x1, z1)
        second.xy1 = (x1, z0)
    elif mode_index == 1:
        first.xy1 = (x0, z0)
        second.xy1 = (x1, z0)
    else:
        first.xy1 = (x0, z0)
        second.xy1 = (x1, z0)


def render_2x2_animation(
    binding_result: Path,
    reference_result: Path,
    output_path: Path,
    *,
    frame_count: int,
    fps: int,
    trail_duration_s: float,
    preview_path: Path | None,
) -> tuple[Path, Path, tuple[LocalMode, ...]]:
    if frame_count < 60 or fps <= 0:
        raise ValueError("At least 60 frames and a positive FPS are required.")
    if trail_duration_s <= 0.0:
        raise ValueError("trail_duration_s must be positive.")

    _require_unmodulated_trajectory(binding_result)
    _require_unmodulated_trajectory(reference_result)
    binding = _load_trajectory(binding_result)
    reference = _load_trajectory(reference_result)
    geometry = _load_geometry(binding_result)
    with np.load(
        reference.result_dir / "velocity_and_wall_shear_field.npz",
        allow_pickle=False,
    ) as reference_field:
        reference_hash = str(
            np.asarray(reference_field["continuous_geometry_hash_sha256"]).item()
        )
    binding_hash = str(np.asarray(geometry["geometry_hash"]).item())
    if reference_hash != binding_hash:
        raise ValueError("Reference and binding results use different vessel geometry.")

    modes = (
        _select_bifurcation_mode(reference, frame_count),
        _select_wall_following_mode(reference, frame_count),
        _select_binding_mode(binding, frame_count),
    )
    record_maps = tuple(_record_map(mode.data, mode.focal_id) for mode in modes)
    available_durations_s = tuple(
        int(np.ptp(mode.frames)) * mode.data.dt_s for mode in modes
    )
    shared_duration_s = min(available_durations_s)
    if shared_duration_s <= 0.0:
        raise ValueError("Every trajectory mode needs a positive physical duration.")
    source_timelines = tuple(
        _linear_source_frames(mode, frame_count, shared_duration_s)
        for mode in modes
    )
    paths = tuple(
        _mode_path(mode, record_map, source_frames)
        for mode, record_map, source_frames in zip(
            modes, record_maps, source_timelines, strict=True
        )
    )

    wall_start = geometry["wall_start"]
    wall_end = geometry["wall_end"]
    global_wall_segments = np.stack((wall_start, wall_end), axis=1)
    positive = geometry["target_positive"]
    target_start = geometry["target_start"][positive]
    target_end = geometry["target_end"][positive]
    global_target_segments = np.stack((target_start, target_end), axis=1)
    local_wall_segments = tuple(
        _segments_near_path(
            wall_start,
            wall_end,
            path,
            max(1.35 * mode.half_width_um, 20.0),
        )
        for mode, path in zip(modes, paths, strict=True)
    )
    local_target_segments = _segments_near_path(
        target_start,
        target_end,
        paths[2],
        max(1.35 * modes[2].half_width_um, 20.0),
    )

    _configure_style()
    mpl.rcParams.update(
        {
            "axes.facecolor": BACKGROUND_COLOR,
            "figure.facecolor": BACKGROUND_COLOR,
            "savefig.facecolor": BACKGROUND_COLOR,
        }
    )
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13.4, 10.0),
        dpi=105,
        gridspec_kw={"hspace": 0.27, "wspace": 0.20},
    )
    global_ax = axes[0, 0]
    local_axes = (axes[0, 1], axes[1, 0], axes[1, 1])
    for axis in (global_ax, *local_axes):
        axis.set_facecolor(BACKGROUND_COLOR)

    global_ax.add_collection(
        _smooth_wall_collection(
            global_wall_segments,
            color=WALL_COLOR,
            linewidth=0.45,
            alpha=0.9,
            zorder=1,
        )
    )
    global_ax.add_collection(
        _smooth_wall_collection(
            global_target_segments,
            color=TARGET_COLOR,
            linewidth=1.5,
            alpha=0.95,
            zorder=2,
        )
    )
    x_um = geometry["x_um"]
    z_um = geometry["z_um"]
    extent = (float(x_um[0]), float(x_um[-1]), float(z_um[0]), float(z_um[-1]))
    global_ax.set(
        xlim=(extent[0], extent[1]),
        ylim=(extent[2], extent[3]),
        xlabel="x (µm)",
        ylabel="z (µm)",
        title="Global context and local regions",
    )
    global_ax.set_aspect("equal", adjustable="box")

    mode_colors = (NORMAL_COLOR, WALL_FOLLOWING_COLOR, BINDING_COLOR)
    local_titles = (
        "Normal flow: bifurcation selection",
        "Non-targeted near-wall motion",
        "Targeted molecular binding",
    )
    global_full_paths = []
    global_trails = []
    global_bubbles = []
    rectangles = []
    for mode, path, color in zip(modes, paths, mode_colors, strict=True):
        full_path, = global_ax.plot(
            path[:, 0],
            path[:, 1],
            color=color,
            linewidth=0.9,
            linestyle=":",
            alpha=0.45,
            zorder=3,
        )
        trail, = global_ax.plot(
            [], [], color=color, linewidth=2.0, alpha=0.95, zorder=4
        )
        bubble = global_ax.scatter(
            [],
            [],
            s=58,
            marker="o",
            facecolors=BUBBLE_FACE_COLOR,
            edgecolors=BUBBLE_EDGE_COLOR,
            linewidths=1.1,
            zorder=6,
        )
        rectangle = Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            fill=False,
            edgecolor=color,
            linewidth=1.2,
            linestyle="--",
            zorder=5,
        )
        global_ax.add_patch(rectangle)
        global_full_paths.append(full_path)
        global_trails.append(trail)
        global_bubbles.append(bubble)
        rectangles.append(rectangle)

    local_full_paths = []
    local_trails = []
    local_other_bubbles = []
    local_focal_bubbles = []
    local_orientation_lines = []
    local_orientation_tips = []
    local_status = []
    for index, (axis, mode, path, color, wall_segments, title) in enumerate(
        zip(
            local_axes,
            modes,
            paths,
            mode_colors,
            local_wall_segments,
            local_titles,
            strict=True,
        )
    ):
        axis.add_collection(
            _smooth_wall_collection(
                wall_segments,
                color=WALL_COLOR,
                linewidth=1.25,
                alpha=0.95,
                zorder=1,
            )
        )
        if index == 2:
            axis.add_collection(
                _smooth_wall_collection(
                    local_target_segments,
                    color=TARGET_COLOR,
                    linewidth=3.0,
                    alpha=0.95,
                    zorder=2,
                )
            )
        full_path, = axis.plot(
            path[:, 0],
            path[:, 1],
            color=color,
            linewidth=1.0,
            linestyle=":",
            alpha=0.42,
            zorder=3,
        )
        trail, = axis.plot(
            [], [], color=color, linewidth=2.2, alpha=0.95, zorder=4
        )
        others = axis.scatter(
            [],
            [],
            s=34,
            marker="o",
            facecolors=OTHER_BUBBLE_COLOR,
            edgecolors="white",
            linewidths=0.55,
            zorder=5,
        )
        focal = axis.scatter(
            [],
            [],
            s=190,
            marker="o",
            facecolors=BUBBLE_FACE_COLOR,
            edgecolors=BUBBLE_EDGE_COLOR,
            linewidths=1.6,
            zorder=8,
        )
        orientation, = axis.plot(
            [],
            [],
            color="white",
            linewidth=2.4,
            solid_capstyle="round",
            zorder=9,
        )
        tip, = axis.plot(
            [],
            [],
            marker="o",
            markersize=3.3,
            color=BUBBLE_EDGE_COLOR,
            linestyle="none",
            zorder=10,
        )
        status = axis.text(
            0.02,
            0.98,
            "",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.2,
            color="#172B4D",
            bbox={
                "boxstyle": "round,pad=0.3",
                "facecolor": BACKGROUND_COLOR,
                "edgecolor": "#AEB8C2",
                "alpha": 0.94,
            },
            zorder=12,
        )
        axis.set(title=title, xlabel="x (µm)", ylabel="z (µm)")
        axis.set_aspect("equal", adjustable="box")
        local_full_paths.append(full_path)
        local_trails.append(trail)
        local_other_bubbles.append(others)
        local_focal_bubbles.append(focal)
        local_orientation_lines.append(orientation)
        local_orientation_tips.append(tip)
        local_status.append(status)

    binding_ring = local_axes[2].scatter(
        [],
        [],
        s=270,
        marker="o",
        facecolors="none",
        edgecolors=BINDING_RING_COLOR,
        linewidths=2.0,
        zorder=7,
    )
    bond_shadow = LineCollection(
        [], colors=FOCAL_EDGE_COLOR, linewidths=3.8, alpha=0.75, zorder=6
    )
    bond_lines = LineCollection(
        [], colors=BOND_COLOR, linewidths=2.2, alpha=1.0, zorder=7
    )
    local_axes[2].add_collection(bond_shadow)
    local_axes[2].add_collection(bond_lines)
    branch_marker, = local_axes[0].plot(
        [],
        [],
        marker="X",
        markersize=7,
        markerfacecolor=NORMAL_COLOR,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linestyle="none",
        zorder=7,
    )
    branch_text = local_axes[0].text(
        0.0, 0.0, "branch decision", fontsize=8, color="#172B4D", zorder=8
    )

    connection_endpoints = (
        ((0.0, 1.0), (0.0, 0.0)),
        ((0.0, 1.0), (1.0, 1.0)),
        ((0.0, 1.0), (1.0, 1.0)),
    )
    connection_pairs: list[tuple[ConnectionPatch, ConnectionPatch]] = []
    for index, (axis, color, endpoints) in enumerate(
        zip(local_axes, mode_colors, connection_endpoints, strict=True)
    ):
        pair = []
        for endpoint in endpoints:
            connection = ConnectionPatch(
                xyA=(0.0, 0.0),
                xyB=endpoint,
                coordsA="data",
                coordsB="axes fraction",
                axesA=global_ax,
                axesB=axis,
                arrowstyle="-",
                linewidth=0.85,
                linestyle="--",
                color=color,
                alpha=0.62,
                clip_on=False,
                zorder=0.5,
            )
            figure.add_artist(connection)
            pair.append(connection)
        connection_pairs.append((pair[0], pair[1]))

    figure.legend(
        handles=[
            Line2D([0], [0], color=WALL_COLOR, linewidth=1.5, label="smooth vessel edge"),
            Line2D([0], [0], color=NORMAL_COLOR, linewidth=2.0, label="normal branch"),
            Line2D(
                [0], [0], color=WALL_FOLLOWING_COLOR, linewidth=2.0, label="wall following"
            ),
            Line2D([0], [0], color=BINDING_COLOR, linewidth=2.0, label="binding path"),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=BUBBLE_FACE_COLOR,
                markeredgecolor=BUBBLE_EDGE_COLOR,
                markersize=8,
                label="microbubble",
            ),
            Line2D([0], [0], color=BOND_COLOR, linewidth=2.2, label="molecular bonds"),
        ],
        loc="lower center",
        ncol=6,
        frameon=False,
        bbox_to_anchor=(0.5, 0.012),
    )
    figure.suptitle("Typical microbubble trajectory modes", fontsize=14, y=0.985)
    figure.subplots_adjust(left=0.065, right=0.97, bottom=0.085, top=0.945)
    figure.canvas.draw()

    def update(output_index: int):
        for mode_index, mode in enumerate(modes):
            data = mode.data
            record_map = record_maps[mode_index]
            source_frame = int(source_timelines[mode_index][int(output_index)])
            focal_record = record_map.get(source_frame)
            if focal_record is None:
                continue
            centre = data.positions_xz_um[focal_record]
            angle = float(data.rotation_angle_rad[focal_record])
            gap = float(data.wall_gap_um[focal_record])
            speed = _actual_speed_um_s(mode, record_map, source_frame)
            expected_bonds = max(float(data.bond_count[focal_record]), 0.0)
            radius = 0.5 * float(data.diameters_um[focal_record])

            half_width = float(mode.half_width_um)
            view_centre = (
                centre
                if mode.fixed_centre_xz_um is None
                else np.asarray(mode.fixed_centre_xz_um, dtype=float)
            )
            axis = local_axes[mode_index]
            axis.set_xlim(view_centre[0] - half_width, view_centre[0] + half_width)
            axis.set_ylim(view_centre[1] - half_width, view_centre[1] + half_width)

            rectangle = rectangles[mode_index]
            rectangle.set_xy(
                (view_centre[0] - half_width, view_centre[1] - half_width)
            )
            rectangle.set_width(2.0 * half_width)
            rectangle.set_height(2.0 * half_width)
            rectangle_bounds = (
                view_centre[0] - half_width,
                view_centre[0] + half_width,
                view_centre[1] - half_width,
                view_centre[1] + half_width,
            )
            _update_connection_pair(
                connection_pairs[mode_index], rectangle_bounds, mode_index
            )

            local_focal_bubbles[mode_index].set_offsets(centre[None, :])
            global_bubbles[mode_index].set_offsets(centre[None, :])

            frame_start = int(data.offsets[source_frame])
            frame_end = int(data.offsets[source_frame + 1])
            rows = np.arange(frame_start, frame_end, dtype=np.int64)
            other = rows[data.bubble_ids[rows] != mode.focal_id]
            local_other_bubbles[mode_index].set_offsets(
                data.positions_xz_um[other]
            )

            trail_frames = max(2, int(round(trail_duration_s / data.dt_s)))
            lower_frame = max(
                int(np.min(mode.frames)), source_frame - trail_frames
            )
            trail_rows = [
                record_map[frame]
                for frame in range(lower_frame, source_frame + 1)
                if frame in record_map
            ]
            trail = data.positions_xz_um[np.asarray(trail_rows, dtype=np.int64)]
            local_trails[mode_index].set_data(trail[:, 0], trail[:, 1])
            global_trails[mode_index].set_data(trail[:, 0], trail[:, 1])

            orientation_start, orientation_end = _orientation_segment(
                axis, centre, angle
            )
            local_orientation_lines[mode_index].set_data(
                [orientation_start[0], orientation_end[0]],
                [orientation_start[1], orientation_end[1]],
            )
            local_orientation_tips[mode_index].set_data(
                [orientation_end[0]], [orientation_end[1]]
            )

            if mode_index == 0:
                branch = np.asarray(mode.branch_position_xz_um, dtype=float)
                branch_marker.set_data([branch[0]], [branch[1]])
                branch_text.set_position((branch[0] + 3.0, branch[1] + 3.0))
                local_status[mode_index].set_text(
                    f"t = {source_frame * data.dt_s:.3f} s | type: normal bifurcation\n"
                    f"speed = {speed:.1f} µm/s | angle = {np.degrees(angle):.1f}°"
                )
            elif mode_index == 1:
                local_status[mode_index].set_text(
                    f"t = {source_frame * data.dt_s:.3f} s | type: non-targeted wall motion\n"
                    f"speed = {speed:.1f} µm/s | angle = {np.degrees(angle):.1f}°"
                )
            else:
                binding_ring.set_offsets(centre[None, :])
                displayed_bonds = _bond_segments(
                    centre,
                    radius,
                    gap,
                    data.wall_normal_xz[focal_record],
                    expected_bonds,
                    float(data.bond_mean_extension_um[focal_record]),
                )
                bond_shadow.set_segments(displayed_bonds)
                bond_lines.set_segments(displayed_bonds)
                local_status[mode_index].set_text(
                    f"t = {source_frame * data.dt_s:.3f} s | type: targeted binding\n"
                    f"speed = {speed:.1f} µm/s | angle = {np.degrees(angle):.1f}°"
                )
        return ()

    animation = FuncAnimation(
        figure,
        update,
        frames=np.arange(frame_count, dtype=np.int64),
        interval=1000.0 / float(fps),
        blit=False,
        repeat=True,
    )
    destination = Path(output_path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    animation.save(
        destination,
        writer=PillowWriter(
            fps=int(fps),
            metadata={
                "title": "Typical microbubble trajectory modes in a 2x2 layout",
                "artist": "ULM microbubble trajectory generator",
            },
        ),
        dpi=105,
    )

    binding_records = np.asarray(
        [record_maps[2][int(frame)] for frame in source_timelines[2]],
        dtype=np.int64,
    )
    preview_index = int(np.argmax(binding.bond_count[binding_records]))
    update(preview_index)
    preview = (
        destination.with_name(destination.stem + "_peak.png")
        if preview_path is None
        else Path(preview_path).resolve()
    )
    figure.savefig(preview, dpi=180, bbox_inches="tight")
    plt.close(figure)
    return destination, preview, modes


def main() -> None:
    args = _parse_args()
    output, preview, modes = render_2x2_animation(
        args.binding_result,
        args.reference_result,
        args.output,
        frame_count=int(args.frames),
        fps=int(args.fps),
        trail_duration_s=float(args.trail_duration_s),
        preview_path=args.preview,
    )
    print(f"Saved GIF: {output}")
    print(f"Saved 2x2 preview: {preview}")
    print(
        f"Modes: {modes[0].detail}; {modes[1].detail}; {modes[2].detail}"
    )


if __name__ == "__main__":
    main()
