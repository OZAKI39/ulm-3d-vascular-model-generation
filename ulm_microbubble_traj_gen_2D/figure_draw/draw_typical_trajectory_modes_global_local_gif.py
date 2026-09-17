"""Render three typical microbubble trajectory modes in one global/local GIF.

The three sequential local scenes are selected from recorded numerical data:

1. normal flow through a bifurcation,
2. non-targeted near-wall motion with zero molecular bonds,
3. target-positive molecular binding and rolling.

The left panel provides global context.  The right panel switches between a
fixed bifurcation view and tracking views for the two wall-interaction modes.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle

try:
    from .draw_molecular_binding_global_local_gif import (
        BOND_COLOR,
        FOCAL_EDGE_COLOR,
        TARGET_COLOR,
        _bond_segments,
        _configure_style,
        _frame_for_records,
        _metadata_float,
        _require,
    )
except ImportError:
    from draw_molecular_binding_global_local_gif import (
        BOND_COLOR,
        FOCAL_EDGE_COLOR,
        TARGET_COLOR,
        _bond_segments,
        _configure_style,
        _frame_for_records,
        _metadata_float,
        _require,
    )


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent
DEFAULT_BINDING_RESULT = SCRIPT_DIR / "outputs" / "molecular_binding_demo_data"
DEFAULT_REFERENCE_RESULT = PACKAGE_DIR / "results" / "20260721_150913"
DEFAULT_OUTPUT = SCRIPT_DIR / "outputs" / "molecular_binding_global_local.gif"

SOLID_COLOR = "#EEE9E2"
LUMEN_COLOR = "#DCEFF7"
WALL_COLOR = "#65717C"
OTHER_BUBBLE_COLOR = "#2A7AB0"
NORMAL_COLOR = "#16877A"
WALL_FOLLOWING_COLOR = "#6C5CE7"
BINDING_COLOR = "#E4572E"
TRAIL_COLOR = "#174E6F"
GAP_COLOR = "#4E5968"


@dataclass(frozen=True)
class TrajectoryData:
    result_dir: Path
    offsets: np.ndarray
    bubble_ids: np.ndarray
    positions_xz_um: np.ndarray
    velocities_xz_um_s: np.ndarray
    realized_velocities_xz_um_s: np.ndarray
    diameters_um: np.ndarray
    wall_gap_um: np.ndarray
    wall_normal_xz: np.ndarray
    rotation_angle_rad: np.ndarray
    vessel_id: np.ndarray
    bond_count: np.ndarray
    bond_force_xz_pn: np.ndarray
    bond_mean_extension_um: np.ndarray
    target_reaction_area_um2: np.ndarray
    dt_s: float
    event_bubble_id: np.ndarray
    event_time_s: np.ndarray
    event_from_vessel_id: np.ndarray
    event_to_vessel_id: np.ndarray
    event_position_xz_um: np.ndarray


@dataclass(frozen=True)
class LocalMode:
    key: str
    title: str
    color: str
    data: TrajectoryData
    focal_id: int
    frames: np.ndarray
    half_width_um: float
    fixed_centre_xz_um: np.ndarray | None
    detail: str
    show_target: bool = False
    show_bonds: bool = False
    show_gap: bool = False
    branch_position_xz_um: np.ndarray | None = None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a global/local GIF containing normal bifurcation selection, "
            "non-targeted wall following, and molecular binding."
        )
    )
    parser.add_argument(
        "--binding-result",
        type=Path,
        default=DEFAULT_BINDING_RESULT,
        help="Result containing the target-positive binding trajectory.",
    )
    parser.add_argument(
        "--reference-result",
        type=Path,
        default=DEFAULT_REFERENCE_RESULT,
        help=(
            "Full-diagnostics non-binding result used for the normal-flow and "
            "non-targeted near-wall scenes."
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--frames",
        type=int,
        default=165,
        help="Total GIF frames divided between the three modes.",
    )
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument(
        "--trail-duration-s",
        type=float,
        default=0.10,
        help="Visible trajectory history in tracking scenes.",
    )
    parser.add_argument(
        "--preview",
        type=Path,
        default=None,
        help="Optional PNG path for the maximum-force binding frame.",
    )
    return parser.parse_args()


def _optional_array(
    archive: np.lib.npyio.NpzFile,
    key: str,
    shape: tuple[int, ...],
    dtype: np.dtype | type,
) -> np.ndarray:
    if key in archive:
        return np.asarray(archive[key], dtype=dtype)
    return np.zeros(shape, dtype=dtype)


def _load_trajectory(result_dir: Path) -> TrajectoryData:
    result = Path(result_dir).resolve()
    path = result / "microbubble_field_trajectories.npz"
    if not path.is_file():
        raise FileNotFoundError(f"Missing trajectory archive: {path}")
    with np.load(path, allow_pickle=False) as trajectory:
        required = {
            "frame_offsets",
            "record_bubble_id",
            "record_positions_um",
            "record_velocities_um_s",
            "record_realized_velocities_um_s",
            "record_diameter_um",
            "record_wall_gap_um",
            "record_wall_normal_xz",
            "record_rotation_angle_rad",
            "record_vessel_id",
        }
        _require(trajectory, required, f"Trajectory {path}")
        record_count = int(trajectory["record_bubble_id"].size)
        event_count = int(
            trajectory["topological_event_bubble_id"].size
            if "topological_event_bubble_id" in trajectory
            else 0
        )
        return TrajectoryData(
            result_dir=result,
            offsets=np.asarray(trajectory["frame_offsets"], dtype=np.int64),
            bubble_ids=np.asarray(trajectory["record_bubble_id"], dtype=np.int64),
            positions_xz_um=np.asarray(
                trajectory["record_positions_um"][:, [0, 2]], dtype=float
            ),
            velocities_xz_um_s=np.asarray(
                trajectory["record_velocities_um_s"][:, [0, 2]], dtype=float
            ),
            realized_velocities_xz_um_s=np.asarray(
                trajectory["record_realized_velocities_um_s"][:, [0, 2]],
                dtype=float,
            ),
            diameters_um=np.asarray(
                trajectory["record_diameter_um"], dtype=float
            ),
            wall_gap_um=np.asarray(
                trajectory["record_wall_gap_um"], dtype=float
            ),
            wall_normal_xz=np.asarray(
                trajectory["record_wall_normal_xz"], dtype=float
            ),
            rotation_angle_rad=np.asarray(
                trajectory["record_rotation_angle_rad"], dtype=float
            ),
            vessel_id=np.asarray(trajectory["record_vessel_id"], dtype=np.int64),
            bond_count=_optional_array(
                trajectory, "record_bond_count_expected", (record_count,), float
            ),
            bond_force_xz_pn=_optional_array(
                trajectory, "record_bond_force_xz_pn", (record_count, 2), float
            ),
            bond_mean_extension_um=_optional_array(
                trajectory,
                "record_bond_mean_tangential_extension_um",
                (record_count,),
                float,
            ),
            target_reaction_area_um2=_optional_array(
                trajectory,
                "record_target_reaction_area_um2",
                (record_count,),
                float,
            ),
            dt_s=_metadata_float(trajectory, "dt_s", 1.0),
            event_bubble_id=_optional_array(
                trajectory, "topological_event_bubble_id", (event_count,), np.int64
            ),
            event_time_s=_optional_array(
                trajectory, "topological_event_time_s", (event_count,), float
            ),
            event_from_vessel_id=_optional_array(
                trajectory,
                "topological_event_from_vessel_id",
                (event_count,),
                np.int64,
            ),
            event_to_vessel_id=_optional_array(
                trajectory,
                "topological_event_to_vessel_id",
                (event_count,),
                np.int64,
            ),
            event_position_xz_um=_optional_array(
                trajectory,
                "topological_event_position_xz_um",
                (event_count, 2),
                float,
            ),
        )


def _sample_frames(start: int, end: int, count: int) -> np.ndarray:
    if end < start:
        raise ValueError("Animation frame interval is reversed.")
    available = end - start + 1
    core_count = min(max(count - 6, 2), available)
    core = np.unique(
        np.rint(np.linspace(start, end, core_count)).astype(np.int64)
    )
    return np.concatenate(
        (
            np.repeat(core[0], 3),
            core,
            np.repeat(core[-1], 3),
        )
    )[:count]


def _record_map(data: TrajectoryData, bubble_id: int) -> dict[int, int]:
    records = np.flatnonzero(data.bubble_ids == int(bubble_id))
    frames = _frame_for_records(data.offsets, records)
    return {
        int(frame): int(record)
        for frame, record in zip(frames, records, strict=True)
    }


def _select_bifurcation_mode(
    data: TrajectoryData, frame_count: int
) -> LocalMode:
    if data.event_bubble_id.size == 0:
        raise ValueError("Reference trajectory has no topological transition events.")
    parent_destinations: dict[int, set[int]] = {}
    for source, destination in zip(
        data.event_from_vessel_id, data.event_to_vessel_id, strict=True
    ):
        parent_destinations.setdefault(int(source), set()).add(int(destination))
    branch_parents = [
        parent
        for parent, destinations in parent_destinations.items()
        if len(destinations) >= 2
    ]
    if not branch_parents:
        raise ValueError("No parent vessel with multiple observed choices was found.")
    root_like_parent = int(np.min(data.event_from_vessel_id))
    non_root_branch_parents = [
        candidate for candidate in branch_parents if candidate != root_like_parent
    ]
    if non_root_branch_parents:
        branch_parents = non_root_branch_parents
    parent = max(
        branch_parents,
        key=lambda candidate: int(
            np.count_nonzero(data.event_from_vessel_id == candidate)
        ),
    )
    destinations = sorted(parent_destinations[parent])
    destination_counts = {
        destination: int(
            np.count_nonzero(
                (data.event_from_vessel_id == parent)
                & (data.event_to_vessel_id == destination)
            )
        )
        for destination in destinations
    }
    chosen_destination = min(
        destinations, key=lambda destination: (destination_counts[destination], destination)
    )
    candidates = np.flatnonzero(
        (data.event_from_vessel_id == parent)
        & (data.event_to_vessel_id == chosen_destination)
    )
    event_index = int(candidates[0])
    focal_id = int(data.event_bubble_id[event_index])
    event_frame = int(round(float(data.event_time_s[event_index]) / data.dt_s))
    record_frames = np.asarray(sorted(_record_map(data, focal_id)), dtype=np.int64)
    start = max(int(record_frames[0]), event_frame - 60)
    end = min(int(record_frames[-1]), event_frame + 60)
    return LocalMode(
        key="bifurcation",
        title="Normal flow: bifurcation selection",
        color=NORMAL_COLOR,
        data=data,
        focal_id=focal_id,
        frames=_sample_frames(start, end, frame_count),
        half_width_um=68.0,
        fixed_centre_xz_um=np.asarray(
            data.event_position_xz_um[event_index], dtype=float
        ),
        detail=f"vessel {parent} -> {chosen_destination}",
        branch_position_xz_um=np.asarray(
            data.event_position_xz_um[event_index], dtype=float
        ),
    )


def _longest_true_run(frames: np.ndarray, condition: np.ndarray) -> tuple[int, int]:
    best_start = -1
    best_end = -2
    start = -1
    for index, selected in enumerate(condition):
        if selected and start < 0:
            start = index
        run_ends = start >= 0 and (
            not selected
            or index == condition.size - 1
            or (
                index + 1 < condition.size
                and frames[index + 1] != frames[index] + 1
            )
        )
        if not run_ends:
            continue
        end = index if selected else index - 1
        if end - start > best_end - best_start:
            best_start, best_end = start, end
        start = -1
    return best_start, best_end


def _select_wall_following_mode(
    data: TrajectoryData, frame_count: int
) -> LocalMode:
    best: tuple[int, int, int, int] | None = None
    for bubble_id in np.unique(data.bubble_ids):
        records = np.flatnonzero(data.bubble_ids == bubble_id)
        frames = _frame_for_records(data.offsets, records)
        near_wall = (
            (data.wall_gap_um[records] <= 0.30)
            & (data.bond_count[records] <= 1.0e-10)
            & (data.target_reaction_area_um2[records] <= 1.0e-10)
        )
        local_start, local_end = _longest_true_run(frames, near_wall)
        if local_start < 0:
            continue
        length = local_end - local_start + 1
        candidate = (
            length,
            int(bubble_id),
            int(frames[local_start]),
            int(frames[local_end]),
        )
        if best is None or candidate[0] > best[0]:
            best = candidate
    if best is None or best[0] < 10:
        raise ValueError("No sustained non-targeted near-wall trajectory was found.")
    length, focal_id, start, end = best
    return LocalMode(
        key="wall_following",
        title="Non-targeted near-wall motion",
        color=WALL_FOLLOWING_COLOR,
        data=data,
        focal_id=focal_id,
        frames=_sample_frames(start, end, frame_count),
        half_width_um=12.0,
        fixed_centre_xz_um=None,
        detail=f"gap <= 0.30 um for {length * data.dt_s:.3f} s; bonds = 0",
        show_gap=True,
    )


def _select_binding_mode(
    data: TrajectoryData, frame_count: int
) -> LocalMode:
    peak_record = int(np.argmax(np.nan_to_num(data.bond_count)))
    if float(data.bond_count[peak_record]) <= 1.0e-8:
        raise ValueError("Binding result contains no nonzero molecular bonds.")
    focal_id = int(data.bubble_ids[peak_record])
    records = np.flatnonzero(data.bubble_ids == focal_id)
    frames = _frame_for_records(data.offsets, records)
    positive = data.bond_count[records] > 1.0e-6
    positive_frames = frames[positive]
    start = max(int(frames[0]), int(positive_frames[0]) - 35)
    end = min(int(frames[-1]), int(positive_frames[-1]) + 15)
    return LocalMode(
        key="binding",
        title="Targeted molecular binding and rolling",
        color=BINDING_COLOR,
        data=data,
        focal_id=focal_id,
        frames=_sample_frames(start, end, frame_count),
        half_width_um=12.0,
        fixed_centre_xz_um=None,
        detail="target-positive wall; deterministic mean-field bonds",
        show_target=True,
        show_bonds=True,
    )


def _load_geometry(binding_result: Path) -> dict[str, np.ndarray]:
    field_path = Path(binding_result).resolve() / "velocity_and_wall_shear_field.npz"
    target_path = Path(binding_result).resolve() / "molecular_target_field.npz"
    with np.load(field_path, allow_pickle=False) as field:
        _require(
            field,
            {
                "x_coordinates_um",
                "z_coordinates_um",
                "lumen_mask",
                "continuous_geometry_hash_sha256",
                "continuous_wall_start_xz_um",
                "continuous_wall_end_xz_um",
            },
            "Accepted field",
        )
        geometry = {
            "x_um": np.asarray(field["x_coordinates_um"], dtype=float),
            "z_um": np.asarray(field["z_coordinates_um"], dtype=float),
            "lumen": np.asarray(field["lumen_mask"], dtype=bool),
            "geometry_hash": np.asarray(field["continuous_geometry_hash_sha256"]),
            "wall_start": np.asarray(
                field["continuous_wall_start_xz_um"], dtype=float
            ),
            "wall_end": np.asarray(
                field["continuous_wall_end_xz_um"], dtype=float
            ),
        }
    with np.load(target_path, allow_pickle=False) as target:
        _require(
            target,
            {"wall_start_xz_um", "wall_end_xz_um", "wall_target_positive"},
            "Molecular target field",
        )
        geometry["target_start"] = np.asarray(
            target["wall_start_xz_um"], dtype=float
        )
        geometry["target_end"] = np.asarray(
            target["wall_end_xz_um"], dtype=float
        )
        geometry["target_positive"] = np.asarray(
            target["wall_target_positive"], dtype=bool
        )
    return geometry


def render_typical_modes_animation(
    binding_result: Path,
    reference_result: Path,
    output_path: Path,
    *,
    total_frames: int,
    fps: int,
    trail_duration_s: float,
    preview_path: Path | None,
) -> tuple[Path, Path, tuple[LocalMode, ...]]:
    if total_frames < 45 or fps <= 0:
        raise ValueError("At least 45 total frames and a positive FPS are required.")
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

    counts = [total_frames // 3] * 3
    for index in range(total_frames % 3):
        counts[index] += 1
    modes = (
        _select_bifurcation_mode(reference, counts[0]),
        _select_wall_following_mode(reference, counts[1]),
        _select_binding_mode(binding, counts[2]),
    )
    sequence = [
        (mode_index, int(frame))
        for mode_index, mode in enumerate(modes)
        for frame in mode.frames
    ]
    mode_record_maps = tuple(_record_map(mode.data, mode.focal_id) for mode in modes)
    mode_paths = []
    for mode, record_map in zip(modes, mode_record_maps, strict=True):
        lo = int(np.min(mode.frames))
        hi = int(np.max(mode.frames))
        path_rows = [
            record_map[frame]
            for frame in range(lo, hi + 1)
            if frame in record_map
        ]
        mode_paths.append(mode.data.positions_xz_um[np.asarray(path_rows, dtype=int)])

    _configure_style()
    mpl.rcParams["axes.facecolor"] = SOLID_COLOR
    figure, (global_ax, local_ax) = plt.subplots(
        1,
        2,
        figsize=(13.2, 6.5),
        dpi=105,
        gridspec_kw={"width_ratios": (1.28, 1.0), "wspace": 0.15},
    )
    x_um = geometry["x_um"]
    z_um = geometry["z_um"]
    lumen = geometry["lumen"]
    extent = (float(x_um[0]), float(x_um[-1]), float(z_um[0]), float(z_um[-1]))
    region_cmap = mpl.colors.ListedColormap((SOLID_COLOR, LUMEN_COLOR))
    for axis in (global_ax, local_ax):
        axis.imshow(
            lumen.T,
            origin="lower",
            extent=extent,
            cmap=region_cmap,
            interpolation="nearest",
            zorder=0,
        )
    wall_segments = np.stack(
        (geometry["wall_start"], geometry["wall_end"]), axis=1
    )
    positive = geometry["target_positive"]
    target_segments = np.stack(
        (
            geometry["target_start"][positive],
            geometry["target_end"][positive],
        ),
        axis=1,
    )
    global_ax.add_collection(
        LineCollection(wall_segments, colors=WALL_COLOR, linewidths=0.30, zorder=1)
    )
    local_ax.add_collection(
        LineCollection(wall_segments, colors=WALL_COLOR, linewidths=1.0, zorder=1)
    )
    global_target = LineCollection(
        target_segments,
        colors=TARGET_COLOR,
        linewidths=1.35,
        alpha=0.95,
        zorder=2,
    )
    local_target = LineCollection(
        target_segments,
        colors=TARGET_COLOR,
        linewidths=3.0,
        alpha=0.95,
        zorder=2,
    )
    global_ax.add_collection(global_target)
    local_ax.add_collection(local_target)

    global_bubbles = global_ax.scatter(
        [], [], s=[], c=[], edgecolors="white", linewidths=0.25, zorder=4
    )
    global_path, = global_ax.plot(
        [], [], color=TRAIL_COLOR, linewidth=1.3, alpha=0.9, zorder=3
    )
    global_focus, = global_ax.plot(
        [],
        [],
        marker="o",
        markersize=7,
        markerfacecolor="none",
        markeredgecolor=FOCAL_EDGE_COLOR,
        markeredgewidth=1.2,
        linestyle="none",
        zorder=5,
    )
    zoom_box = Rectangle(
        (0.0, 0.0),
        1.0,
        1.0,
        fill=False,
        edgecolor=FOCAL_EDGE_COLOR,
        linewidth=1.0,
        linestyle="--",
        zorder=5,
    )
    global_ax.add_patch(zoom_box)
    mode_label = global_ax.text(
        0.02,
        0.98,
        "",
        transform=global_ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="#172B4D",
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "edgecolor": "#B8C4CE",
            "alpha": 0.92,
        },
        zorder=8,
    )
    global_ax.set(
        xlim=(extent[0], extent[1]),
        ylim=(extent[2], extent[3]),
        xlabel="x (µm)",
        ylabel="z (µm)",
        title="Global vascular context",
    )
    global_ax.set_aspect("equal", adjustable="box")

    local_other_bubbles = local_ax.scatter(
        [], [], s=[], c=[], edgecolors="white", linewidths=0.5, zorder=4
    )
    local_full_path, = local_ax.plot(
        [],
        [],
        color=TRAIL_COLOR,
        linewidth=1.0,
        linestyle=":",
        alpha=0.45,
        zorder=3,
    )
    local_trail, = local_ax.plot(
        [], [], color=TRAIL_COLOR, linewidth=2.0, alpha=0.9, zorder=4
    )
    focal_circle = Circle(
        (0.0, 0.0),
        radius=1.0,
        facecolor=NORMAL_COLOR,
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
    gap_line, = local_ax.plot(
        [],
        [],
        color=GAP_COLOR,
        linewidth=1.8,
        linestyle="--",
        zorder=6,
    )
    bond_shadow = LineCollection(
        [], colors=FOCAL_EDGE_COLOR, linewidths=3.8, alpha=0.75, zorder=5
    )
    bonds = LineCollection(
        [], colors=BOND_COLOR, linewidths=2.2, alpha=1.0, zorder=6
    )
    local_ax.add_collection(bond_shadow)
    local_ax.add_collection(bonds)
    branch_marker, = local_ax.plot(
        [],
        [],
        marker="X",
        markersize=7,
        markerfacecolor=NORMAL_COLOR,
        markeredgecolor="white",
        markeredgewidth=0.8,
        linestyle="none",
        zorder=6,
    )
    branch_annotation = local_ax.text(
        0.0,
        0.0,
        "",
        color="#172B4D",
        fontsize=8,
        ha="left",
        va="bottom",
        zorder=7,
    )
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
    local_ax.set(xlabel="x (µm)", ylabel="z (µm)")
    local_ax.set_aspect("equal", adjustable="box")

    figure.legend(
        handles=[
            Patch(facecolor=LUMEN_COLOR, edgecolor=WALL_COLOR, label="lumen"),
            Patch(facecolor=SOLID_COLOR, edgecolor=WALL_COLOR, label="external solid"),
            Line2D(
                [0], [0], color=TARGET_COLOR, linewidth=3, label="target-positive wall"
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=OTHER_BUBBLE_COLOR,
                markeredgecolor="white",
                markersize=7,
                label="microbubble",
            ),
            Line2D([0], [0], color=BOND_COLOR, linewidth=2.2, label="molecular bonds"),
        ],
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )
    figure.suptitle(
        "Typical microbubble trajectory modes",
        fontsize=13,
        y=0.985,
    )
    figure.subplots_adjust(bottom=0.10, top=0.93)

    def update(item: tuple[int, int]):
        mode_index, frame = item
        mode = modes[int(mode_index)]
        data = mode.data
        record_map = mode_record_maps[int(mode_index)]
        focal_record = record_map.get(int(frame))
        if focal_record is None:
            return ()
        start = int(data.offsets[frame])
        end = int(data.offsets[frame + 1])
        rows = np.arange(start, end, dtype=np.int64)
        centre = data.positions_xz_um[focal_record]
        radius = 0.5 * float(data.diameters_um[focal_record])
        angle = float(data.rotation_angle_rad[focal_record])
        gap = float(data.wall_gap_um[focal_record])
        speed = float(np.linalg.norm(data.velocities_xz_um_s[focal_record]))
        expected_bonds = max(float(data.bond_count[focal_record]), 0.0)
        force = float(np.linalg.norm(data.bond_force_xz_pn[focal_record]))

        global_bubbles.set_offsets(data.positions_xz_um[rows])
        global_bubbles.set_sizes(
            np.clip(5.0 * data.diameters_um[rows], 5.0, 17.0)
        )
        bubble_colors = np.full(rows.size, OTHER_BUBBLE_COLOR, dtype=object)
        bubble_colors[data.bond_count[rows] > 1.0e-6] = BINDING_COLOR
        global_bubbles.set_color(bubble_colors)
        global_focus.set_data([centre[0]], [centre[1]])
        global_path.set_data(mode_paths[mode_index][:, 0], mode_paths[mode_index][:, 1])
        global_path.set_color(mode.color)
        mode_label.set_text(f"Mode {mode_index + 1}/3\n{mode.title}")

        half_width = float(mode.half_width_um)
        view_centre = (
            centre
            if mode.fixed_centre_xz_um is None
            else np.asarray(mode.fixed_centre_xz_um, dtype=float)
        )
        zoom_box.set_xy(
            (view_centre[0] - half_width, view_centre[1] - half_width)
        )
        zoom_box.set_width(2.0 * half_width)
        zoom_box.set_height(2.0 * half_width)
        local_ax.set_xlim(view_centre[0] - half_width, view_centre[0] + half_width)
        local_ax.set_ylim(view_centre[1] - half_width, view_centre[1] + half_width)
        local_ax.set_title(mode.title)
        global_target.set_visible(mode.show_target)
        local_target.set_visible(mode.show_target)

        other = rows[data.bubble_ids[rows] != mode.focal_id]
        local_other_bubbles.set_offsets(data.positions_xz_um[other])
        local_other_bubbles.set_sizes(
            np.clip(18.0 * data.diameters_um[other] ** 2, 15.0, 95.0)
        )
        local_other_bubbles.set_color(OTHER_BUBBLE_COLOR)
        local_full_path.set_data(
            mode_paths[mode_index][:, 0], mode_paths[mode_index][:, 1]
        )
        local_full_path.set_color(mode.color)

        trail_frames = max(2, int(round(trail_duration_s / data.dt_s)))
        trail_start = max(int(np.min(mode.frames)), int(frame) - trail_frames)
        trail_rows = [
            record_map[source_frame]
            for source_frame in range(trail_start, int(frame) + 1)
            if source_frame in record_map
        ]
        trail_positions = data.positions_xz_um[np.asarray(trail_rows, dtype=int)]
        local_trail.set_data(trail_positions[:, 0], trail_positions[:, 1])
        local_trail.set_color(mode.color)

        focal_circle.center = (float(centre[0]), float(centre[1]))
        focal_circle.set_radius(radius)
        focal_circle.set_facecolor(mode.color)
        orientation = np.asarray([np.cos(angle), np.sin(angle)])
        line_start = centre - 0.72 * radius * orientation
        line_end = centre + 0.72 * radius * orientation
        orientation_line.set_data(
            [line_start[0], line_end[0]], [line_start[1], line_end[1]]
        )
        orientation_tip.set_data([line_end[0]], [line_end[1]])

        normal = np.asarray(data.wall_normal_xz[focal_record], dtype=float)
        normal_norm = float(np.linalg.norm(normal))
        if mode.show_gap and normal_norm > 0.0:
            normal /= normal_norm
            surface = centre - normal * radius
            wall = centre - normal * (radius + max(gap, 0.0))
            gap_line.set_data([surface[0], wall[0]], [surface[1], wall[1]])
        else:
            gap_line.set_data([], [])

        displayed_bonds = (
            _bond_segments(
                centre,
                radius,
                gap,
                normal,
                expected_bonds,
                float(data.bond_mean_extension_um[focal_record]),
            )
            if mode.show_bonds
            else []
        )
        bond_shadow.set_segments(displayed_bonds)
        bonds.set_segments(displayed_bonds)

        if mode.branch_position_xz_um is not None:
            branch = np.asarray(mode.branch_position_xz_um, dtype=float)
            branch_marker.set_data([branch[0]], [branch[1]])
            branch_annotation.set_position((branch[0] + 3.0, branch[1] + 3.0))
            branch_annotation.set_text("branch decision")
        else:
            branch_marker.set_data([], [])
            branch_annotation.set_text("")

        if mode.key == "bifurcation":
            status = (
                f"t = {frame * data.dt_s:.3f} s   |   NORMAL FLOW\n"
                f"{mode.detail}   speed = {speed:.1f} µm/s\n"
                f"rotation = {np.degrees(angle):.1f}°   "
                f"wall gap = {gap:.2f} µm"
            )
        elif mode.key == "wall_following":
            status = (
                f"t = {frame * data.dt_s:.3f} s   |   NON-TARGETED\n"
                f"near-wall speed = {speed:.1f} µm/s   "
                f"rotation = {np.degrees(angle):.1f}°\n"
                f"wall gap = {gap:.3f} µm   expected bonds = 0"
            )
        else:
            status = (
                f"t = {frame * data.dt_s:.3f} s   |   TARGETED BINDING\n"
                f"speed = {speed:.1f} µm/s   "
                f"rotation = {np.degrees(angle):.1f}°\n"
                f"expected bonds = {expected_bonds:.2f}   "
                f"|F_bond| = {force:.2f} pN   gap = {gap:.3f} µm"
            )
        status_text.set_text(status)
        return ()

    animation = FuncAnimation(
        figure,
        update,
        frames=sequence,
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
                "title": "Typical microbubble trajectory modes",
                "artist": "ULM microbubble trajectory generator",
            },
        ),
        dpi=105,
    )

    binding_mode_index = 2
    binding_force = np.linalg.norm(binding.bond_force_xz_pn, axis=1)
    preview_record = int(np.argmax(binding_force))
    preview_frame = int(
        _frame_for_records(
            binding.offsets, np.asarray([preview_record], dtype=np.int64)
        )[0]
    )
    update((binding_mode_index, preview_frame))
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
    output, preview, modes = render_typical_modes_animation(
        args.binding_result,
        args.reference_result,
        args.output,
        total_frames=int(args.frames),
        fps=int(args.fps),
        trail_duration_s=float(args.trail_duration_s),
        preview_path=args.preview,
    )
    print(f"Saved GIF: {output}")
    print(f"Saved binding preview: {preview}")
    for index, mode in enumerate(modes, start=1):
        print(
            f"Mode {index}: {mode.title}; bubble={mode.focal_id}; "
            f"frames={int(np.min(mode.frames))}..{int(np.max(mode.frames))}; "
            f"{mode.detail}"
        )


if __name__ == "__main__":
    main()
