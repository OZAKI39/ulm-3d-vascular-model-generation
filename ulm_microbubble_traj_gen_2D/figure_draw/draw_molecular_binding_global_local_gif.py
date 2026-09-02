"""Render a global + tracking-zoom GIF of microbubble molecular binding.

The animation reads the append-only trajectory archive written by
``generate_microbubble_trajectories.py``.  It automatically selects the bubble
with the largest expected bond count and shows its translation, accumulated
rotation, target-wall contact, and mean-field bonds in one synchronized view.
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
from matplotlib.patches import Circle, Rectangle


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULT_DIR = SCRIPT_DIR / "outputs" / "molecular_binding_demo_data"
DEFAULT_OUTPUT = SCRIPT_DIR / "outputs" / "molecular_binding_global_local.gif"

WALL_COLOR = "#7B8794"
LUMEN_COLOR = "#EAF4F8"
TARGET_COLOR = "#C23B72"
FREE_BUBBLE_COLOR = "#2A7AB0"
BOUND_BUBBLE_COLOR = "#E4572E"
FOCAL_EDGE_COLOR = "#173F5F"
BOND_COLOR = "#FFD000"
TRAIL_COLOR = "#195E8B"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a synchronized global and local molecular-binding GIF from "
            "one completed trajectory result."
        )
    )
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=DEFAULT_RESULT_DIR,
        help="Completed result directory.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="GIF output path.",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=110,
        help="Number of sampled animation frames.",
    )
    parser.add_argument("--fps", type=int, default=15, help="GIF playback rate.")
    parser.add_argument(
        "--zoom-half-width-um",
        type=float,
        default=12.0,
        help="Half-width of the tracking local view.",
    )
    parser.add_argument(
        "--trail-duration-s",
        type=float,
        default=0.10,
        help="Focal-bubble trail duration.",
    )
    parser.add_argument(
        "--preview",
        type=Path,
        default=None,
        help=(
            "Optional peak-binding PNG. Defaults to the GIF name with "
            "'_peak.png' appended."
        ),
    )
    return parser.parse_args()


def _metadata_float(
    trajectory: np.lib.npyio.NpzFile, key: str, default: float
) -> float:
    if "metadata_keys" not in trajectory or "metadata_values" not in trajectory:
        return float(default)
    keys = np.asarray(trajectory["metadata_keys"], dtype=str)
    values = np.asarray(trajectory["metadata_values"], dtype=str)
    match = np.flatnonzero(keys == key)
    return float(values[int(match[0])]) if match.size else float(default)


def _require(archive: np.lib.npyio.NpzFile, keys: set[str], label: str) -> None:
    missing = sorted(keys.difference(archive.files))
    if missing:
        raise KeyError(f"{label} is missing arrays: {missing}")


def _frame_for_records(offsets: np.ndarray, record_indices: np.ndarray) -> np.ndarray:
    return np.searchsorted(offsets, record_indices, side="right") - 1


def _select_event(
    offsets: np.ndarray,
    bubble_ids: np.ndarray,
    bond_count: np.ndarray,
    requested_frames: int,
) -> tuple[int, np.ndarray, int]:
    finite_bonds = np.nan_to_num(bond_count, nan=0.0, posinf=0.0, neginf=0.0)
    peak_record = int(np.argmax(finite_bonds))
    peak_bonds = float(finite_bonds[peak_record])
    if peak_bonds <= 1.0e-8:
        raise ValueError(
            "No molecular-binding event was recorded: maximum expected bond "
            f"count is {peak_bonds:.6g}."
        )
    focal_id = int(bubble_ids[peak_record])
    focal_records = np.flatnonzero(bubble_ids == focal_id)
    focal_frames = _frame_for_records(offsets, focal_records)
    focal_bonds = finite_bonds[focal_records]
    positive = focal_bonds > max(1.0e-6, peak_bonds * 1.0e-4)
    positive_frames = focal_frames[positive]
    first_binding = int(positive_frames[0])
    last_binding = int(positive_frames[-1])
    total_frames = offsets.size - 1
    pre_roll = max(25, int(round(0.08 * (last_binding - first_binding + 1))))
    post_roll = max(15, pre_roll // 2)
    start_frame = max(0, first_binding - pre_roll)
    end_frame = min(total_frames - 1, last_binding + post_roll)
    count = min(max(int(requested_frames), 2), end_frame - start_frame + 1)
    sampled_frames = np.unique(
        np.rint(np.linspace(start_frame, end_frame, count)).astype(np.int64)
    )
    peak_frame = int(
        _frame_for_records(offsets, np.asarray([peak_record], dtype=np.int64))[0]
    )
    return focal_id, sampled_frames, peak_frame


def _bond_segments(
    centre_xz: np.ndarray,
    radius_um: float,
    gap_um: float,
    inward_normal_xz: np.ndarray,
    expected_count: float,
    mean_tangential_extension_um: float,
) -> list[np.ndarray]:
    if expected_count <= 1.0e-6:
        return []
    normal = np.asarray(inward_normal_xz, dtype=float)
    normal_norm = float(np.linalg.norm(normal))
    if not np.isfinite(normal_norm) or normal_norm <= 0.0:
        return []
    normal /= normal_norm
    tangent = np.asarray([-normal[1], normal[0]], dtype=float)
    count = min(9, max(1, int(np.ceil(expected_count))))
    span_um = min(0.62 * radius_um, 0.75)
    offsets_um = (
        np.asarray([0.0])
        if count == 1
        else np.linspace(-span_um, span_um, count)
    )
    wall_anchor = centre_xz - normal * (radius_um + max(gap_um, 0.0))
    tangential_extension_um = float(
        np.clip(mean_tangential_extension_um, -1.5 * radius_um, 1.5 * radius_um)
    )
    segments: list[np.ndarray] = []
    for offset_um in offsets_um:
        normal_radius_um = np.sqrt(
            max(radius_um * radius_um - float(offset_um) ** 2, 0.0)
        )
        bubble_anchor = (
            centre_xz
            - normal * normal_radius_um
            + tangent * float(offset_um)
        )
        target_anchor = (
            wall_anchor
            + tangent
            * (float(offset_um) - tangential_extension_um)
        )
        segments.append(np.vstack((bubble_anchor, target_anchor)))
    return segments


def _configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.linewidth": 0.8,
            "axes.edgecolor": "#344054",
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def render_binding_animation(
    result_dir: Path,
    output_path: Path,
    *,
    frame_count: int,
    fps: int,
    zoom_half_width_um: float,
    trail_duration_s: float,
    preview_path: Path | None,
) -> tuple[Path, Path, dict[str, float | int]]:
    result = Path(result_dir).resolve()
    field_path = result / "velocity_and_wall_shear_field.npz"
    trajectory_path = result / "microbubble_field_trajectories.npz"
    target_path = result / "molecular_target_field.npz"
    if not field_path.is_file() or not trajectory_path.is_file():
        raise FileNotFoundError(
            "Result directory must contain the accepted field and trajectory NPZ files."
        )
    if not target_path.is_file():
        raise FileNotFoundError(f"Missing molecular target field: {target_path}")
    if frame_count < 2 or fps <= 0:
        raise ValueError("frames must be at least 2 and fps must be positive.")
    if zoom_half_width_um <= 0.0 or trail_duration_s <= 0.0:
        raise ValueError("Zoom width and trail duration must be positive.")

    with np.load(field_path, allow_pickle=False) as field:
        _require(
            field,
            {
                "x_coordinates_um",
                "z_coordinates_um",
                "lumen_mask",
                "continuous_wall_start_xz_um",
                "continuous_wall_end_xz_um",
            },
            "Accepted field",
        )
        x_um = np.asarray(field["x_coordinates_um"], dtype=float)
        z_um = np.asarray(field["z_coordinates_um"], dtype=float)
        lumen = np.asarray(field["lumen_mask"], dtype=bool)
        wall_start = np.asarray(
            field["continuous_wall_start_xz_um"], dtype=float
        )
        wall_end = np.asarray(field["continuous_wall_end_xz_um"], dtype=float)

    with np.load(target_path, allow_pickle=False) as target:
        _require(
            target,
            {
                "wall_start_xz_um",
                "wall_end_xz_um",
                "wall_target_positive",
            },
            "Molecular target field",
        )
        target_start = np.asarray(target["wall_start_xz_um"], dtype=float)
        target_end = np.asarray(target["wall_end_xz_um"], dtype=float)
        target_positive = np.asarray(target["wall_target_positive"], dtype=bool)

    with np.load(trajectory_path, allow_pickle=False) as trajectory:
        required_trajectory = {
            "frame_offsets",
            "record_bubble_id",
            "record_positions_um",
            "record_velocities_um_s",
            "record_diameter_um",
            "record_wall_gap_um",
            "record_wall_normal_xz",
            "record_rotation_angle_rad",
            "record_bond_count_expected",
            "record_bond_force_xz_pn",
            "record_bond_mean_tangential_extension_um",
        }
        _require(trajectory, required_trajectory, "Trajectory")
        offsets = np.asarray(trajectory["frame_offsets"], dtype=np.int64)
        bubble_ids = np.asarray(trajectory["record_bubble_id"], dtype=np.int64)
        positions_xz = np.asarray(
            trajectory["record_positions_um"][:, [0, 2]], dtype=float
        )
        velocities_xz = np.asarray(
            trajectory["record_velocities_um_s"][:, [0, 2]], dtype=float
        )
        diameters_um = np.asarray(trajectory["record_diameter_um"], dtype=float)
        wall_gap_um = np.asarray(trajectory["record_wall_gap_um"], dtype=float)
        wall_normal = np.asarray(trajectory["record_wall_normal_xz"], dtype=float)
        rotation_rad = np.asarray(
            trajectory["record_rotation_angle_rad"], dtype=float
        )
        bond_count = np.asarray(
            trajectory["record_bond_count_expected"], dtype=float
        )
        bond_force = np.asarray(
            trajectory["record_bond_force_xz_pn"], dtype=float
        )
        bond_mean_extension_um = np.asarray(
            trajectory["record_bond_mean_tangential_extension_um"], dtype=float
        )
        dt_s = _metadata_float(trajectory, "dt_s", 1.0)

    focal_id, animation_frames, peak_frame = _select_event(
        offsets, bubble_ids, bond_count, frame_count
    )
    focal_records = np.flatnonzero(bubble_ids == focal_id)
    focal_frames = _frame_for_records(offsets, focal_records)
    focal_record_by_frame = {
        int(frame): int(record)
        for frame, record in zip(focal_frames, focal_records, strict=True)
    }
    peak_bond_record = int(
        focal_records[np.argmax(np.nan_to_num(bond_count[focal_records]))]
    )
    focal_force_magnitude = np.linalg.norm(bond_force[focal_records], axis=1)
    representative_record = int(
        focal_records[np.argmax(np.nan_to_num(focal_force_magnitude))]
    )
    representative_frame = int(
        _frame_for_records(
            offsets, np.asarray([representative_record], dtype=np.int64)
        )[0]
    )
    trail_frames = max(2, int(round(trail_duration_s / dt_s)))

    _configure_style()
    figure, (global_ax, local_ax) = plt.subplots(
        1,
        2,
        figsize=(13.2, 6.5),
        dpi=105,
        gridspec_kw={"width_ratios": (1.28, 1.0), "wspace": 0.15},
    )
    wall_segments = np.stack((wall_start, wall_end), axis=1)
    target_segments = np.stack(
        (target_start[target_positive], target_end[target_positive]), axis=1
    )

    extent = (float(x_um[0]), float(x_um[-1]), float(z_um[0]), float(z_um[-1]))
    global_ax.imshow(
        lumen.T,
        origin="lower",
        extent=extent,
        cmap=mpl.colors.ListedColormap(("white", LUMEN_COLOR)),
        interpolation="nearest",
        alpha=0.95,
        zorder=0,
    )
    global_ax.add_collection(
        LineCollection(wall_segments, colors=WALL_COLOR, linewidths=0.26, zorder=1)
    )
    global_ax.add_collection(
        LineCollection(
            target_segments,
            colors=TARGET_COLOR,
            linewidths=1.35,
            alpha=0.95,
            zorder=2,
        )
    )
    global_bubbles = global_ax.scatter(
        [],
        [],
        s=[],
        c=[],
        edgecolors="white",
        linewidths=0.25,
        zorder=4,
    )
    global_trail, = global_ax.plot(
        [], [], color=TRAIL_COLOR, linewidth=1.25, alpha=0.9, zorder=3
    )
    global_focus, = global_ax.plot(
        [],
        [],
        marker="o",
        markersize=7.0,
        markerfacecolor="none",
        markeredgecolor=FOCAL_EDGE_COLOR,
        markeredgewidth=1.2,
        linestyle="none",
        zorder=5,
    )
    zoom_box = Rectangle(
        (0.0, 0.0),
        2.0 * zoom_half_width_um,
        2.0 * zoom_half_width_um,
        fill=False,
        edgecolor=FOCAL_EDGE_COLOR,
        linewidth=1.0,
        linestyle="--",
        zorder=5,
    )
    global_ax.add_patch(zoom_box)
    global_ax.set(
        xlim=(extent[0], extent[1]),
        ylim=(extent[2], extent[3]),
        xlabel="x (µm)",
        ylabel="z (µm)",
        title="Global vascular context",
    )
    global_ax.set_aspect("equal", adjustable="box")

    local_ax.set_facecolor(LUMEN_COLOR)
    local_ax.add_collection(
        LineCollection(wall_segments, colors=WALL_COLOR, linewidths=1.0, zorder=1)
    )
    local_ax.add_collection(
        LineCollection(
            target_segments,
            colors=TARGET_COLOR,
            linewidths=3.0,
            alpha=0.95,
            zorder=2,
        )
    )
    local_other_bubbles = local_ax.scatter(
        [],
        [],
        s=[],
        c=[],
        edgecolors="white",
        linewidths=0.5,
        zorder=4,
    )
    local_trail, = local_ax.plot(
        [], [], color=TRAIL_COLOR, linewidth=2.0, alpha=0.85, zorder=3
    )
    focal_circle = Circle(
        (0.0, 0.0),
        radius=1.0,
        facecolor=FREE_BUBBLE_COLOR,
        edgecolor=FOCAL_EDGE_COLOR,
        linewidth=1.4,
        zorder=7,
    )
    local_ax.add_patch(focal_circle)
    orientation_line, = local_ax.plot(
        [], [], color="white", linewidth=2.2, solid_capstyle="round", zorder=8
    )
    orientation_tip, = local_ax.plot(
        [],
        [],
        marker="o",
        markersize=3.2,
        color=FOCAL_EDGE_COLOR,
        linestyle="none",
        zorder=9,
    )
    bond_shadow_collection = LineCollection(
        [],
        colors=FOCAL_EDGE_COLOR,
        linewidths=3.8,
        alpha=0.75,
        zorder=5,
    )
    bond_collection = LineCollection(
        [], colors=BOND_COLOR, linewidths=2.2, alpha=1.0, zorder=6
    )
    local_ax.add_collection(bond_shadow_collection)
    local_ax.add_collection(bond_collection)
    status_text = local_ax.text(
        0.02,
        0.98,
        "",
        transform=local_ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="#172B4D",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#B8C4CE",
            "alpha": 0.92,
        },
        zorder=10,
    )
    local_ax.set(
        xlabel="x (µm)",
        ylabel="z (µm)",
        title=f"Tracking zoom — focal microbubble ID {focal_id}",
    )
    local_ax.set_aspect("equal", adjustable="box")

    legend = [
        Line2D(
            [0],
            [0],
            color=TARGET_COLOR,
            linewidth=3,
            label="target-positive wall",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=FREE_BUBBLE_COLOR,
            markeredgecolor="white",
            markersize=7,
            label="moving microbubble",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=BOUND_BUBBLE_COLOR,
            markeredgecolor=FOCAL_EDGE_COLOR,
            markersize=7,
            label="bonded focal microbubble",
        ),
        Line2D(
            [0],
            [0],
            color=BOND_COLOR,
            markeredgecolor=FOCAL_EDGE_COLOR,
            linewidth=2.2,
            label="mean-field bonds",
        ),
    ]
    figure.legend(
        handles=legend,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )
    figure.suptitle(
        "Microbubble translation, rotation and molecular binding",
        fontsize=13,
        y=0.985,
    )
    figure.subplots_adjust(bottom=0.10, top=0.93)

    def update(frame_index: int):
        frame = int(frame_index)
        start = int(offsets[frame])
        end = int(offsets[frame + 1])
        rows = np.arange(start, end, dtype=np.int64)
        current_positions = positions_xz[rows]
        current_bonds = bond_count[rows]
        current_colors = np.where(
            current_bonds > 1.0e-6, BOUND_BUBBLE_COLOR, FREE_BUBBLE_COLOR
        )
        global_bubbles.set_offsets(current_positions)
        global_bubbles.set_sizes(np.clip(5.0 * diameters_um[rows], 5.0, 17.0))
        global_bubbles.set_color(current_colors)

        focal_record = focal_record_by_frame.get(frame)
        if focal_record is None:
            global_focus.set_data([], [])
            local_other_bubbles.set_offsets(np.empty((0, 2)))
            focal_circle.set_visible(False)
            orientation_line.set_data([], [])
            orientation_tip.set_data([], [])
            bond_shadow_collection.set_segments([])
            bond_collection.set_segments([])
            status_text.set_text(f"t = {frame * dt_s:.3f} s\nfocal bubble inactive")
            return ()

        focal_circle.set_visible(True)
        centre = positions_xz[focal_record]
        radius = 0.5 * float(diameters_um[focal_record])
        bonds = max(float(bond_count[focal_record]), 0.0)
        angle = float(rotation_rad[focal_record])
        gap = float(wall_gap_um[focal_record])
        force_pn = float(np.linalg.norm(bond_force[focal_record]))
        speed_um_s = float(np.linalg.norm(velocities_xz[focal_record]))
        is_bound = bonds > 1.0e-6

        global_focus.set_data([centre[0]], [centre[1]])
        zoom_box.set_xy(
            (
                centre[0] - zoom_half_width_um,
                centre[1] - zoom_half_width_um,
            )
        )
        local_ax.set_xlim(
            centre[0] - zoom_half_width_um, centre[0] + zoom_half_width_um
        )
        local_ax.set_ylim(
            centre[1] - zoom_half_width_um, centre[1] + zoom_half_width_um
        )

        other = rows[bubble_ids[rows] != focal_id]
        local_other_bubbles.set_offsets(positions_xz[other])
        local_other_bubbles.set_sizes(
            np.clip(18.0 * diameters_um[other] ** 2, 15.0, 95.0)
        )
        local_other_bubbles.set_color(
            np.where(
                bond_count[other] > 1.0e-6,
                BOUND_BUBBLE_COLOR,
                FREE_BUBBLE_COLOR,
            )
        )
        focal_circle.center = (float(centre[0]), float(centre[1]))
        focal_circle.set_radius(radius)
        focal_circle.set_facecolor(
            BOUND_BUBBLE_COLOR if is_bound else FREE_BUBBLE_COLOR
        )

        orientation = np.asarray([np.cos(angle), np.sin(angle)], dtype=float)
        line_start = centre - 0.72 * radius * orientation
        line_end = centre + 0.72 * radius * orientation
        orientation_line.set_data(
            [line_start[0], line_end[0]], [line_start[1], line_end[1]]
        )
        orientation_tip.set_data([line_end[0]], [line_end[1]])
        displayed_bonds = _bond_segments(
            centre,
            radius,
            gap,
            wall_normal[focal_record],
            bonds,
            float(bond_mean_extension_um[focal_record]),
        )
        bond_shadow_collection.set_segments(displayed_bonds)
        bond_collection.set_segments(displayed_bonds)

        trail_start_frame = max(int(focal_frames[0]), frame - trail_frames)
        trail_rows = [
            focal_record_by_frame[f]
            for f in range(trail_start_frame, frame + 1)
            if f in focal_record_by_frame
        ]
        trail_positions = positions_xz[np.asarray(trail_rows, dtype=np.int64)]
        global_trail.set_data(trail_positions[:, 0], trail_positions[:, 1])
        local_trail.set_data(trail_positions[:, 0], trail_positions[:, 1])

        state = "BONDED" if is_bound else "approaching target"
        status_text.set_text(
            f"t = {frame * dt_s:.3f} s   |   {state}\n"
            f"speed = {speed_um_s:.1f} µm/s   "
            f"rotation = {np.degrees(angle):.1f}°\n"
            f"expected bonds = {bonds:.2f}   "
            f"|F_bond| = {force_pn:.2f} pN   "
            f"gap = {gap:.3f} µm"
        )
        return ()

    animation = FuncAnimation(
        figure,
        update,
        frames=animation_frames,
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
                "title": "Microbubble molecular binding: global and local views",
                "artist": "ULM microbubble trajectory generator",
            },
        ),
        dpi=105,
    )

    peak_preview = (
        destination.with_name(destination.stem + "_peak.png")
        if preview_path is None
        else Path(preview_path).resolve()
    )
    update(representative_frame)
    figure.savefig(peak_preview, dpi=180, bbox_inches="tight")
    plt.close(figure)
    summary: dict[str, float | int] = {
        "focal_bubble_id": focal_id,
        "peak_frame": peak_frame,
        "peak_time_s": peak_frame * dt_s,
        "peak_expected_bonds": float(bond_count[peak_bond_record]),
        "representative_frame": representative_frame,
        "representative_time_s": representative_frame * dt_s,
        "peak_bond_force_pn": float(
            np.linalg.norm(bond_force[representative_record])
        ),
        "animation_frames": int(animation_frames.size),
    }
    return destination, peak_preview, summary


def main() -> None:
    args = _parse_args()
    output, preview, summary = render_binding_animation(
        args.result_dir,
        args.output,
        frame_count=int(args.frames),
        fps=int(args.fps),
        zoom_half_width_um=float(args.zoom_half_width_um),
        trail_duration_s=float(args.trail_duration_s),
        preview_path=args.preview,
    )
    print(f"Saved GIF: {output}")
    print(f"Saved peak preview: {preview}")
    print(
        "Selected event: "
        f"bubble={summary['focal_bubble_id']}, "
        f"maximum-force t={summary['representative_time_s']:.3f} s, "
        f"peak expected bonds={summary['peak_expected_bonds']:.3f}, "
        f"peak |F_bond|={summary['peak_bond_force_pn']:.3f} pN"
    )


if __name__ == "__main__":
    main()
