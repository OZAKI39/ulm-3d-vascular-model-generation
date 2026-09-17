"""Load field-based microbubble result folders for visualization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

from ...core.types import FiniteElementVelocityField
from ...flow.hybrid_velocity import sample_finite_element_velocity
from ...geometry.lumen_width import classify_narrow_lumen_cells

from .molecular_target_overlay import MolecularTargetOverlay, load_molecular_target_overlay


PACKAGE_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULTS_ROOT = PACKAGE_ROOT / "results"
_NUMERICAL_WALL_LOCK_SPEED_ATOL_UM_S = 1.0e-12


def _assign_display_lanes_impl(
    offsets: np.ndarray,
    compact_bubble_ids: np.ndarray,
    target: int,
) -> np.ndarray:
    """Map ragged observations onto stable, recyclable display lanes."""

    record_lanes = np.full(compact_bubble_ids.size, -1, dtype=np.int64)
    if target <= 0 or compact_bubble_ids.size == 0:
        return record_lanes

    valid_ids = compact_bubble_ids[compact_bubble_ids >= 0]
    id_count = int(np.max(valid_ids)) + 1 if valid_ids.size else 0
    id_to_lane = np.full(id_count, -1, dtype=np.int64)
    lane_owner = np.full(int(target), -1, dtype=np.int64)
    last_seen_frame = np.full(id_count, -1, dtype=np.int64)

    for frame in range(offsets.size - 1):
        start = int(offsets[frame])
        end = int(offsets[frame + 1])
        for record in range(start, end):
            bubble = int(compact_bubble_ids[record])
            if bubble >= 0:
                last_seen_frame[bubble] = frame

        for lane in range(int(target)):
            bubble = int(lane_owner[lane])
            if bubble >= 0 and last_seen_frame[bubble] != frame:
                id_to_lane[bubble] = -1
                lane_owner[lane] = -1

        for record in range(start, end):
            bubble = int(compact_bubble_ids[record])
            if bubble < 0:
                continue
            lane = int(id_to_lane[bubble])
            if lane < 0:
                for candidate in range(int(target)):
                    if lane_owner[candidate] < 0:
                        lane = candidate
                        break
                if lane < 0:
                    raise ValueError("Trajectory frame population exceeds allocated display lanes.")
                id_to_lane[bubble] = lane
                lane_owner[lane] = bubble
            record_lanes[record] = lane
    return record_lanes


try:
    from numba import njit as _numba_njit
except ImportError:  # pragma: no cover - exercised on installations without Numba.
    _assign_display_lanes = _assign_display_lanes_impl
else:
    _assign_display_lanes = _numba_njit(cache=True)(_assign_display_lanes_impl)


@dataclass(frozen=True)
class FieldVisualizationData:
    """Arrays required to draw vessel background fields and bubble trajectories."""

    result_dir: Path
    x_coordinates_um: np.ndarray
    z_coordinates_um: np.ndarray
    grid_spacing_um: float
    lumen_mask: np.ndarray
    wall_mask: np.ndarray
    wall_distance_um: np.ndarray
    effective_lumen_diameter_px: np.ndarray
    graph_lumen_diameter_px: np.ndarray
    mask_lumen_diameter_px: np.ndarray
    hole_mask: np.ndarray
    narrow_lumen_mask: np.ndarray
    vessel_id_grid: np.ndarray
    initial_velocity_xz_um_s: np.ndarray
    initial_speed_um_s: np.ndarray
    final_velocity_xz_um_s: np.ndarray
    speed_um_s: np.ndarray
    pressure: np.ndarray
    wall_shear_grid_pa: np.ndarray
    wall_shear_display_mask: np.ndarray
    local_shear_grid_pa: np.ndarray
    local_shear_source: str
    positions_um: np.ndarray
    velocities_um_s: np.ndarray
    realized_velocities_um_s: np.ndarray
    bubble_wall_shear_pa: np.ndarray
    bubble_vessel_id: np.ndarray
    active: np.ndarray
    bubble_id: np.ndarray
    bubble_diameter_um: np.ndarray
    bubble_rotation_angle_rad: np.ndarray
    bubble_wall_gap_um: np.ndarray
    bubble_near_wall_weight: np.ndarray
    bubble_wall_contact: np.ndarray
    bubble_bond_count_expected: np.ndarray
    bubble_contact_constraint_active: np.ndarray
    bubble_contact_reaction_force_pn: np.ndarray
    bubble_contact_free_normal_velocity_um_s: np.ndarray
    bubble_contact_constrained_normal_velocity_um_s: np.ndarray
    bubble_numerical_wall_lock: np.ndarray
    numerical_wall_lock_source: str
    active_count_per_frame: np.ndarray
    injected_count_per_frame: np.ndarray
    terminated_count_per_frame: np.ndarray
    peak_active_bubbles: int
    contact_geometry_tolerance_um: float
    wall_gap_invalid_below_um: float
    strict_nonnegative_wall_gap: bool
    near_wall_xi_far: float
    molecular_target: MolecularTargetOverlay


def find_latest_result(results_root: Path = DEFAULT_RESULTS_ROOT) -> Path:
    """Return the newest result folder containing both field and trajectory files."""

    root = Path(results_root)
    if not root.exists():
        raise FileNotFoundError(f"Result root does not exist: {root}")

    candidates = []
    for item in root.iterdir():
        if not item.is_dir():
            continue
        field_path = item / "velocity_and_wall_shear_field.npz"
        traj_path = item / "microbubble_field_trajectories.npz"
        if field_path.exists() and traj_path.exists():
            candidates.append(item)

    if not candidates:
        raise FileNotFoundError(f"No complete results folder found under {root}")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def load_result(result_dir: Path | None = None) -> FieldVisualizationData:
    """Load the selected result folder or the newest complete results folder."""

    selected_dir = find_latest_result() if result_dir is None else Path(result_dir).resolve()
    field_path = selected_dir / "velocity_and_wall_shear_field.npz"
    traj_path = selected_dir / "microbubble_field_trajectories.npz"

    if not field_path.exists():
        raise FileNotFoundError(f"Missing velocity field file: {field_path}")
    if not traj_path.exists():
        raise FileNotFoundError(f"Missing trajectory file: {traj_path}")

    with np.load(field_path) as field, np.load(traj_path) as traj:
        x_coordinates_um = np.asarray(field["x_coordinates_um"], dtype=float)
        z_coordinates_um = np.asarray(field["z_coordinates_um"], dtype=float)
        lumen_mask = np.asarray(field["lumen_mask"], dtype=bool)
        spacing_um = _field_spacing_um(field, x_coordinates_um, z_coordinates_um)
        wall_distance_um = _field_wall_distance_um(field, lumen_mask)
        hole_mask = _hole_mask(lumen_mask)
        radius_um = _field_radius_um(field, lumen_mask)
        junction_mask = _field_junction_mask(field, lumen_mask)
        open_boundary_mask = _field_open_boundary_mask(field, lumen_mask)
        final_velocity = _field_velocity(field, "velocity_xz_um_s", lumen_mask)
        final_speed = np.asarray(field["speed_um_s"], dtype=float)
        initial_velocity = _field_velocity(field, "initial_velocity_xz_um_s", lumen_mask, fallback=final_velocity)
        initial_speed = _field_speed(field, "initial_speed_um_s", initial_velocity, lumen_mask)
        pressure = _field_scalar(field, "pressure", lumen_mask, default_value=0.0)
        local_shear, local_shear_source = _field_local_shear_stress_pa(
            field,
            final_velocity,
            lumen_mask,
            spacing_um,
        )
        narrow_result = classify_narrow_lumen_cells(
            lumen_mask,
            spacing_um,
            radius_um,
            min_diameter_px=8.0,
            junction_mask=junction_mask,
            open_boundary_mask=open_boundary_mask,
        )
        trajectory = _trajectory_arrays(traj)
        molecular_target = load_molecular_target_overlay(
            selected_dir,
            x_coordinates_um,
            z_coordinates_um,
        )
        return FieldVisualizationData(
            result_dir=selected_dir,
            x_coordinates_um=x_coordinates_um,
            z_coordinates_um=z_coordinates_um,
            grid_spacing_um=float(spacing_um),
            lumen_mask=lumen_mask,
            wall_mask=np.asarray(field["wall_mask"], dtype=bool),
            wall_distance_um=wall_distance_um,
            effective_lumen_diameter_px=narrow_result.effective_diameter_px,
            graph_lumen_diameter_px=narrow_result.graph_diameter_px,
            mask_lumen_diameter_px=narrow_result.mask_diameter_px,
            hole_mask=hole_mask,
            narrow_lumen_mask=narrow_result.narrow_mask,
            vessel_id_grid=np.asarray(field["vessel_id"], dtype=int),
            initial_velocity_xz_um_s=initial_velocity,
            initial_speed_um_s=initial_speed,
            final_velocity_xz_um_s=final_velocity,
            speed_um_s=final_speed,
            pressure=pressure,
            wall_shear_grid_pa=np.asarray(field["wall_shear_stress_pa"], dtype=float),
            wall_shear_display_mask=(
                lumen_mask
                & np.asarray(field["wall_mask"], dtype=bool)
                & ~open_boundary_mask
            ),
            local_shear_grid_pa=local_shear,
            local_shear_source=local_shear_source,
            positions_um=trajectory["positions_um"],
            velocities_um_s=trajectory["velocities_um_s"],
            realized_velocities_um_s=trajectory["realized_velocities_um_s"],
            bubble_wall_shear_pa=trajectory["wall_shear_stress_pa"],
            bubble_vessel_id=trajectory["vessel_id"],
            active=trajectory["active"],
            bubble_id=trajectory["bubble_id"],
            bubble_diameter_um=trajectory["bubble_diameter_um"],
            bubble_rotation_angle_rad=trajectory["rotation_angle_rad"],
            bubble_wall_gap_um=trajectory["bubble_wall_gap_um"],
            bubble_near_wall_weight=trajectory["bubble_near_wall_weight"],
            bubble_wall_contact=trajectory["bubble_wall_contact"],
            bubble_bond_count_expected=trajectory["bubble_bond_count_expected"],
            bubble_contact_constraint_active=trajectory["contact_constraint_active"],
            bubble_contact_reaction_force_pn=trajectory["contact_reaction_force_pn"],
            bubble_contact_free_normal_velocity_um_s=trajectory[
                "contact_free_normal_velocity_um_s"
            ],
            bubble_contact_constrained_normal_velocity_um_s=trajectory[
                "contact_constrained_normal_velocity_um_s"
            ],
            bubble_numerical_wall_lock=trajectory["numerical_wall_lock"],
            numerical_wall_lock_source=str(trajectory["numerical_wall_lock_source"]),
            active_count_per_frame=trajectory["active_count_per_frame"],
            injected_count_per_frame=trajectory["injected_count_per_frame"],
            terminated_count_per_frame=trajectory["terminated_count_per_frame"],
            peak_active_bubbles=int(trajectory["peak_active_bubbles"]),
            contact_geometry_tolerance_um=float(
                trajectory["contact_geometry_tolerance_um"]
            ),
            wall_gap_invalid_below_um=float(trajectory["wall_gap_invalid_below_um"]),
            strict_nonnegative_wall_gap=bool(trajectory["strict_nonnegative_wall_gap"]),
            near_wall_xi_far=float(trajectory["near_wall_xi_far"]),
            molecular_target=molecular_target,
        )


def _trajectory_arrays(
    traj: np.lib.npyio.NpzFile,
) -> dict[str, np.ndarray | int | float | bool]:
    """Load the append-only schema or adapt a legacy dense trajectory file."""

    if "frame_offsets" not in traj:
        positions = np.asarray(traj["positions_um"], dtype=float)
        frames, bubbles = positions.shape[:2]
        active = np.asarray(traj["active"], dtype=bool)
        bubble_id = np.broadcast_to(np.arange(bubbles, dtype=np.int64), (frames, bubbles)).copy()
        dt_s = _trajectory_time_step_s(traj)
        rotation_angle_rad = _legacy_rotation_angles(traj, frames, bubbles)
        fixed_diameter_um = _trajectory_metadata_float(traj, "bubble_diameter_um", 2.0)
        legacy_failure_duration_s = (
            np.asarray(traj["wall_slide_failure_duration_s"], dtype=float)
            if "wall_slide_failure_duration_s" in traj
            else None
        )
        if (
            legacy_failure_duration_s is not None
            and legacy_failure_duration_s.shape != (frames, bubbles)
        ):
            raise ValueError(
                "Legacy wall_slide_failure_duration_s must match the dense trajectory shape."
            )
        numerical_wall_lock = (
            legacy_failure_duration_s >= dt_s * (1.0 - 1.0e-9)
            if legacy_failure_duration_s is not None
            else np.zeros((frames, bubbles), dtype=bool)
        )
        return {
            "positions_um": positions,
            "velocities_um_s": np.asarray(traj["velocities_um_s"], dtype=float),
            "realized_velocities_um_s": _realized_center_velocities_um_s(
                positions,
                bubble_id,
                active,
                dt_s,
            ),
            "wall_shear_stress_pa": np.asarray(traj["wall_shear_stress_pa"], dtype=float),
            "vessel_id": np.asarray(traj["vessel_id"], dtype=int),
            "active": active,
            "bubble_id": bubble_id,
            "bubble_diameter_um": np.full((frames, bubbles), fixed_diameter_um, dtype=float),
            "rotation_angle_rad": rotation_angle_rad,
            "bubble_wall_gap_um": np.full((frames, bubbles), np.nan, dtype=float),
            "bubble_near_wall_weight": np.full((frames, bubbles), np.nan, dtype=float),
            "bubble_wall_contact": np.zeros((frames, bubbles), dtype=bool),
            "bubble_bond_count_expected": np.zeros((frames, bubbles), dtype=float),
            "contact_constraint_active": np.zeros((frames, bubbles), dtype=bool),
            "contact_reaction_force_pn": np.zeros((frames, bubbles), dtype=float),
            "contact_free_normal_velocity_um_s": np.full(
                (frames, bubbles),
                np.nan,
                dtype=float,
            ),
            "contact_constrained_normal_velocity_um_s": np.full(
                (frames, bubbles),
                np.nan,
                dtype=float,
            ),
            "numerical_wall_lock": numerical_wall_lock,
            "numerical_wall_lock_source": (
                "legacy_wall_slide_failure_duration"
                if legacy_failure_duration_s is not None
                else "unavailable"
            ),
            "active_count_per_frame": np.count_nonzero(active, axis=1).astype(np.int32),
            "injected_count_per_frame": np.zeros(frames, dtype=np.int32),
            "terminated_count_per_frame": np.zeros(frames, dtype=np.int32),
            "peak_active_bubbles": bubbles,
            "near_wall_xi_far": _trajectory_metadata_float(
                traj, "near_wall_xi_far", 1.0
            ),
            **_wall_gap_metadata(traj),
        }

    offsets = np.asarray(traj["frame_offsets"], dtype=np.int64)
    if offsets.ndim != 1 or offsets.size < 2 or offsets[0] != 0 or np.any(np.diff(offsets) < 0):
        raise ValueError("Trajectory frame_offsets must be a monotone array starting at zero.")
    counts = np.diff(offsets)
    frames = counts.size
    target = int(np.max(counts)) if counts.size else 0
    positions = np.full((frames, target, 3), np.nan, dtype=float)
    velocities = np.zeros((frames, target, 3), dtype=float)
    realized_velocities = np.zeros((frames, target, 3), dtype=float)
    wall_shear = np.zeros((frames, target), dtype=float)
    vessel_id = np.full((frames, target), -1, dtype=int)
    active = np.zeros((frames, target), dtype=bool)
    bubble_id = np.full((frames, target), -1, dtype=np.int64)
    bubble_diameter_um = np.full((frames, target), np.nan, dtype=float)
    rotation_angle_rad = np.zeros((frames, target), dtype=float)
    bubble_wall_gap_um = np.full((frames, target), np.nan, dtype=float)
    bubble_near_wall_weight = np.full((frames, target), np.nan, dtype=float)
    bubble_wall_contact = np.zeros((frames, target), dtype=bool)
    bubble_bond_count_expected = np.zeros((frames, target), dtype=float)
    contact_constraint_active = np.zeros((frames, target), dtype=bool)
    contact_reaction_force_pn = np.zeros((frames, target), dtype=float)
    contact_free_normal_velocity_um_s = np.full((frames, target), np.nan, dtype=float)
    contact_constrained_normal_velocity_um_s = np.full(
        (frames, target),
        np.nan,
        dtype=float,
    )
    legacy_numerical_wall_lock = np.zeros((frames, target), dtype=bool)

    flat_positions = np.asarray(traj["record_positions_um"], dtype=float)
    flat_velocities = np.asarray(traj["record_velocities_um_s"], dtype=float)
    flat_realized_velocities = (
        np.asarray(traj["record_realized_velocities_um_s"], dtype=float)
        if "record_realized_velocities_um_s" in traj
        else None
    )
    flat_shear = np.asarray(traj["record_wall_shear_stress_pa"], dtype=float)
    flat_vessel = np.asarray(traj["record_vessel_id"], dtype=int)
    flat_active = np.asarray(traj["record_active"], dtype=bool)
    flat_bubble_id = np.asarray(traj["record_bubble_id"], dtype=np.int64)
    flat_rotation_angle_rad = (
        np.asarray(traj["record_rotation_angle_rad"], dtype=float)
        if "record_rotation_angle_rad" in traj
        else np.zeros(flat_bubble_id.shape, dtype=float)
    )
    if "record_diameter_um" in traj:
        flat_diameter_um = np.asarray(traj["record_diameter_um"], dtype=float)
    else:
        flat_diameter_um = np.full(
            flat_bubble_id.shape,
            _trajectory_metadata_float(traj, "bubble_diameter_um", 2.0),
            dtype=float,
        )
    flat_wall_gap_um = (
        np.asarray(traj["record_wall_gap_um"], dtype=float)
        if "record_wall_gap_um" in traj
        else np.full(flat_bubble_id.shape, np.nan, dtype=float)
    )
    flat_wall_contact = (
        np.asarray(traj["record_wall_contact"], dtype=bool)
        if "record_wall_contact" in traj
        else np.zeros(flat_bubble_id.shape, dtype=bool)
    )
    flat_near_wall_weight = _optional_record_array(
        traj,
        "record_near_wall_weight",
        flat_bubble_id.shape,
        dtype=float,
        fill_value=np.nan,
    )
    flat_bond_count_expected = _optional_record_array(
        traj,
        "record_bond_count_expected",
        flat_bubble_id.shape,
        dtype=float,
        fill_value=0.0,
    )
    flat_contact_constraint_active = _optional_record_array(
        traj,
        "record_contact_constraint_active",
        flat_bubble_id.shape,
        dtype=bool,
        fill_value=False,
    )
    flat_contact_reaction_force_pn = _optional_record_array(
        traj,
        "record_contact_reaction_force_pn",
        flat_bubble_id.shape,
        dtype=float,
        fill_value=0.0,
    )
    flat_contact_free_normal_velocity_um_s = _optional_record_array(
        traj,
        "record_contact_free_normal_velocity_um_s",
        flat_bubble_id.shape,
        dtype=float,
        fill_value=np.nan,
    )
    flat_contact_constrained_normal_velocity_um_s = _optional_record_array(
        traj,
        "record_contact_constrained_normal_velocity_um_s",
        flat_bubble_id.shape,
        dtype=float,
        fill_value=np.nan,
    )
    flat_legacy_wall_slide_failure_duration_s = (
        _optional_record_array(
            traj,
            "record_wall_slide_failure_duration_s",
            flat_bubble_id.shape,
            dtype=float,
            fill_value=0.0,
        )
        if "record_wall_slide_failure_duration_s" in traj
        else None
    )
    if int(offsets[-1]) != int(flat_positions.shape[0]):
        raise ValueError("Trajectory frame_offsets do not match the number of flat observation records.")
    dt_s = _trajectory_time_step_s(traj)

    # Ragged records are sorted by permanent ID inside each frame, so writing
    # them into columns by row number would shift every later bubble whenever a
    # smaller ID exits.  Keep a live ID-to-display-lane map instead.  Display
    # lanes may be recycled only after an ID disappears; simulation state slots
    # and permanent IDs remain untouched.
    valid_flat_ids = flat_bubble_id >= 0
    compact_bubble_ids = np.full(flat_bubble_id.shape, -1, dtype=np.int64)
    if np.any(valid_flat_ids):
        unique_bubble_ids = np.unique(flat_bubble_id[valid_flat_ids])
        compact_bubble_ids[valid_flat_ids] = np.searchsorted(
            unique_bubble_ids,
            flat_bubble_id[valid_flat_ids],
        )
    record_lanes = _assign_display_lanes(offsets, compact_bubble_ids, target)
    valid_records = np.flatnonzero(record_lanes >= 0)
    record_frames = np.repeat(np.arange(frames, dtype=np.int64), counts)[valid_records]
    display_lanes = record_lanes[valid_records]

    positions[record_frames, display_lanes] = flat_positions[valid_records]
    velocities[record_frames, display_lanes] = flat_velocities[valid_records]
    if flat_realized_velocities is not None:
        realized_velocities[record_frames, display_lanes] = flat_realized_velocities[valid_records]
    wall_shear[record_frames, display_lanes] = flat_shear[valid_records]
    vessel_id[record_frames, display_lanes] = flat_vessel[valid_records]
    active[record_frames, display_lanes] = flat_active[valid_records]
    bubble_id[record_frames, display_lanes] = flat_bubble_id[valid_records]
    bubble_diameter_um[record_frames, display_lanes] = flat_diameter_um[valid_records]
    rotation_angle_rad[record_frames, display_lanes] = flat_rotation_angle_rad[valid_records]
    bubble_wall_gap_um[record_frames, display_lanes] = flat_wall_gap_um[valid_records]
    bubble_near_wall_weight[record_frames, display_lanes] = flat_near_wall_weight[
        valid_records
    ]
    bubble_wall_contact[record_frames, display_lanes] = flat_wall_contact[valid_records]
    bubble_bond_count_expected[record_frames, display_lanes] = (
        flat_bond_count_expected[valid_records]
    )
    contact_constraint_active[record_frames, display_lanes] = flat_contact_constraint_active[valid_records]
    contact_reaction_force_pn[record_frames, display_lanes] = flat_contact_reaction_force_pn[valid_records]
    contact_free_normal_velocity_um_s[record_frames, display_lanes] = (
        flat_contact_free_normal_velocity_um_s[valid_records]
    )
    contact_constrained_normal_velocity_um_s[record_frames, display_lanes] = (
        flat_contact_constrained_normal_velocity_um_s[valid_records]
    )
    if flat_legacy_wall_slide_failure_duration_s is not None:
        legacy_numerical_wall_lock[record_frames, display_lanes] = (
            flat_legacy_wall_slide_failure_duration_s[valid_records]
            >= dt_s * (1.0 - 1.0e-9)
        )

    display_realized_velocities = (
        realized_velocities
        if flat_realized_velocities is not None
        else _realized_center_velocities_um_s(
            positions,
            bubble_id,
            active,
            dt_s,
        )
    )
    if "record_contact_constraint_active" in traj:
        constrained_rhs_speed = np.linalg.norm(velocities[..., (0, 2)], axis=-1)
        realized_speed = np.linalg.norm(display_realized_velocities[..., (0, 2)], axis=-1)
        finite_speed = np.isfinite(constrained_rhs_speed) & np.isfinite(realized_speed)
        numerical_wall_lock = (
            contact_constraint_active
            & active
            & finite_speed
            & (constrained_rhs_speed > _NUMERICAL_WALL_LOCK_SPEED_ATOL_UM_S)
            & (realized_speed <= _NUMERICAL_WALL_LOCK_SPEED_ATOL_UM_S)
        )
        numerical_wall_lock_source = (
            "v16_accepted_internal_velocity_vs_same_id_realized_velocity"
        )
    elif flat_legacy_wall_slide_failure_duration_s is not None:
        numerical_wall_lock = legacy_numerical_wall_lock
        numerical_wall_lock_source = "legacy_wall_slide_failure_duration"
    else:
        numerical_wall_lock = np.zeros((frames, target), dtype=bool)
        numerical_wall_lock_source = "unavailable"
    return {
        "positions_um": positions,
        "velocities_um_s": velocities,
        "realized_velocities_um_s": display_realized_velocities,
        "wall_shear_stress_pa": wall_shear,
        "vessel_id": vessel_id,
        "active": active,
        "bubble_id": bubble_id,
        "bubble_diameter_um": bubble_diameter_um,
        "rotation_angle_rad": rotation_angle_rad,
        "bubble_wall_gap_um": bubble_wall_gap_um,
        "bubble_near_wall_weight": bubble_near_wall_weight,
        "bubble_wall_contact": bubble_wall_contact,
        "bubble_bond_count_expected": bubble_bond_count_expected,
        "contact_constraint_active": contact_constraint_active,
        "contact_reaction_force_pn": contact_reaction_force_pn,
        "contact_free_normal_velocity_um_s": contact_free_normal_velocity_um_s,
        "contact_constrained_normal_velocity_um_s": (
            contact_constrained_normal_velocity_um_s
        ),
        "numerical_wall_lock": numerical_wall_lock,
        "numerical_wall_lock_source": numerical_wall_lock_source,
        "active_count_per_frame": _optional_frame_array(
            traj, "active_count_per_frame", np.count_nonzero(active, axis=1), frames
        ),
        "injected_count_per_frame": _optional_frame_array(
            traj, "injected_count_per_frame", np.zeros(frames, dtype=np.int32), frames
        ),
        "terminated_count_per_frame": _optional_frame_array(
            traj, "terminated_count_per_frame", np.zeros(frames, dtype=np.int32), frames
        ),
        "peak_active_bubbles": int(np.max(np.count_nonzero(active, axis=1))) if frames else 0,
        "near_wall_xi_far": _trajectory_metadata_float(
            traj, "near_wall_xi_far", 1.0
        ),
        **_wall_gap_metadata(traj),
    }


def _realized_center_velocities_um_s(
    positions_um: np.ndarray,
    bubble_id: np.ndarray,
    active: np.ndarray,
    dt_s: float,
) -> np.ndarray:
    """Derive display velocity from same-ID center displacement after geometry constraints."""

    positions = np.asarray(positions_um, dtype=float)
    ids = np.asarray(bubble_id, dtype=np.int64)
    is_active = np.asarray(active, dtype=bool)
    if positions.ndim != 3 or positions.shape[-1] != 3:
        raise ValueError("Trajectory positions must have shape (frames, lanes, 3).")
    if ids.shape != positions.shape[:2] or is_active.shape != positions.shape[:2]:
        raise ValueError("Trajectory IDs and activity must match positions.shape[:2].")
    step_s = float(dt_s)
    if not np.isfinite(step_s) or step_s <= 0.0:
        raise ValueError("Trajectory dt_s must be finite and positive for realized-speed visualization.")

    realized = np.zeros_like(positions, dtype=float)
    if positions.shape[0] < 2:
        return realized

    finite = np.all(np.isfinite(positions), axis=-1)
    same_forward = (
        is_active[:-1]
        & is_active[1:]
        & finite[:-1]
        & finite[1:]
        & (ids[:-1] >= 0)
        & (ids[:-1] == ids[1:])
    )
    forward_velocity = (positions[1:] - positions[:-1]) / step_s
    realized[:-1][same_forward] = forward_velocity[same_forward]

    # The last observation of a permanent ID has no same-ID forward record.
    # Reuse its immediately preceding realized displacement without ever
    # differencing across a replacement ID in the same display lane.
    has_forward = np.zeros_like(is_active)
    has_forward[:-1] = same_forward
    has_backward = np.zeros_like(is_active)
    has_backward[1:] = same_forward
    use_backward = is_active & finite & ~has_forward & has_backward
    realized[1:][use_backward[1:]] = forward_velocity[use_backward[1:]]
    realized[~is_active] = 0.0
    return realized


def _trajectory_time_step_s(traj: np.lib.npyio.NpzFile) -> float:
    """Return a valid saved time step, using one second only for legacy files without metadata."""

    dt_s = _trajectory_metadata_float(traj, "dt_s", 1.0)
    return float(dt_s) if np.isfinite(dt_s) and dt_s > 0.0 else 1.0


def _legacy_rotation_angles(traj: np.lib.npyio.NpzFile, frames: int, bubbles: int) -> np.ndarray:
    """Load an optional dense legacy rotation array or provide a stationary visual marker."""

    if "rotation_angle_rad" in traj:
        values = np.asarray(traj["rotation_angle_rad"], dtype=float)
        if values.shape == (frames, bubbles):
            return values
    return np.zeros((frames, bubbles), dtype=float)


def _optional_frame_array(
    traj: np.lib.npyio.NpzFile,
    key: str,
    fallback: np.ndarray,
    frames: int,
) -> np.ndarray:
    if key in traj:
        values = np.asarray(traj[key], dtype=np.int32)
        if values.shape == (frames,):
            return values
    return np.asarray(fallback, dtype=np.int32)


def _optional_record_array(
    traj: np.lib.npyio.NpzFile,
    key: str,
    expected_shape: tuple[int, ...],
    *,
    dtype,
    fill_value,
) -> np.ndarray:
    """Load one optional flat diagnostic and reject a mismatched saved schema."""

    if key in traj:
        values = np.asarray(traj[key], dtype=dtype)
        if values.shape != expected_shape:
            raise ValueError(
                f"Trajectory {key} must have shape {expected_shape}, got {values.shape}."
            )
        return values
    return np.full(expected_shape, fill_value, dtype=dtype)


def _wall_gap_metadata(traj: np.lib.npyio.NpzFile) -> dict[str, float | bool]:
    """Return current strict feasibility or the historical penetration convention.

    The contact tolerance controls how accurately the feasible-side contact
    point is solved; it never licenses a negative accepted gap.  Older files
    used a penetration allowance, so their visual warning threshold remains
    ``-penetration_tolerance_um`` for faithful retrospective inspection.
    """

    if _trajectory_metadata_has_key(traj, "contact_geometry_tolerance_um"):
        tolerance_um = max(
            0.0,
            _trajectory_metadata_float(traj, "contact_geometry_tolerance_um", 0.0),
        )
        return {
            "contact_geometry_tolerance_um": tolerance_um,
            "wall_gap_invalid_below_um": 0.0,
            "strict_nonnegative_wall_gap": True,
        }

    legacy_tolerance_um = max(
        0.0,
        _trajectory_metadata_float(traj, "penetration_tolerance_um", 0.0),
    )
    return {
        "contact_geometry_tolerance_um": legacy_tolerance_um,
        "wall_gap_invalid_below_um": -legacy_tolerance_um,
        "strict_nonnegative_wall_gap": False,
    }


def _trajectory_metadata_has_key(traj: np.lib.npyio.NpzFile, key: str) -> bool:
    if "metadata_keys" not in traj or "metadata_values" not in traj:
        return False
    keys = np.asarray(traj["metadata_keys"], dtype=str)
    return bool(np.any(keys == str(key)))


def _trajectory_metadata_text(
    traj: np.lib.npyio.NpzFile,
    key: str,
    fallback: str,
) -> str:
    if "metadata_keys" not in traj or "metadata_values" not in traj:
        return str(fallback)
    keys = np.asarray(traj["metadata_keys"], dtype=str)
    values = np.asarray(traj["metadata_values"], dtype=str)
    matches = np.flatnonzero(keys == str(key))
    return str(values[int(matches[0])]) if matches.size else str(fallback)


def _trajectory_metadata_float(traj: np.lib.npyio.NpzFile, key: str, fallback: float) -> float:
    if "metadata_keys" not in traj or "metadata_values" not in traj:
        return float(fallback)
    keys = np.asarray(traj["metadata_keys"], dtype=str)
    values = np.asarray(traj["metadata_values"], dtype=str)
    matches = np.flatnonzero(keys == str(key))
    if matches.size == 0:
        return float(fallback)
    try:
        return float(values[int(matches[0])])
    except (TypeError, ValueError):
        return float(fallback)


def _field_spacing_um(field: np.lib.npyio.NpzFile, x_coordinates_um: np.ndarray, z_coordinates_um: np.ndarray) -> float:
    if "spacing_um" in field:
        spacing = np.asarray(field["spacing_um"], dtype=float).ravel()
        if spacing.size > 0 and float(spacing[0]) > 0.0:
            return float(spacing[0])
    dx = np.diff(np.asarray(x_coordinates_um, dtype=float))
    dz = np.diff(np.asarray(z_coordinates_um, dtype=float))
    candidates = np.concatenate([dx[dx > 0.0], dz[dz > 0.0]])
    if candidates.size == 0:
        return 1.0
    return float(np.median(candidates))


def _field_wall_distance_um(
    field: np.lib.npyio.NpzFile,
    lumen_mask: np.ndarray,
) -> np.ndarray:
    if "distance_to_wall_um" not in field:
        raise ValueError(
            "The field archive is missing the continuous-wall distance array."
        )
    distance = np.asarray(field["distance_to_wall_um"], dtype=float)
    if distance.shape != np.asarray(lumen_mask).shape:
        raise ValueError(
            "distance_to_wall_um does not match the saved lumen grid."
        )
    return distance


def _field_radius_um(field: np.lib.npyio.NpzFile, lumen_mask: np.ndarray) -> np.ndarray:
    if "radius_um" in field:
        radius = np.asarray(field["radius_um"], dtype=float)
        if radius.shape == lumen_mask.shape:
            return radius
    return np.full(lumen_mask.shape, np.nan, dtype=float)


def _field_velocity(
    field: np.lib.npyio.NpzFile,
    key: str,
    lumen_mask: np.ndarray,
    fallback: np.ndarray | None = None,
) -> np.ndarray:
    if key in field:
        velocity = np.asarray(field[key], dtype=float)
        if velocity.shape == (*lumen_mask.shape, 2):
            velocity = velocity.copy()
            velocity[~lumen_mask] = 0.0
            return velocity
    if fallback is not None and fallback.shape == (*lumen_mask.shape, 2):
        velocity = np.asarray(fallback, dtype=float).copy()
        velocity[~lumen_mask] = 0.0
        return velocity
    return np.zeros((*lumen_mask.shape, 2), dtype=float)


def _field_speed(field: np.lib.npyio.NpzFile, key: str, velocity_xz_um_s: np.ndarray, lumen_mask: np.ndarray) -> np.ndarray:
    if key in field:
        speed = np.asarray(field[key], dtype=float)
        if speed.shape == lumen_mask.shape:
            speed = speed.copy()
            speed[~lumen_mask] = 0.0
            return speed
    speed = np.linalg.norm(np.asarray(velocity_xz_um_s, dtype=float), axis=-1)
    speed[~lumen_mask] = 0.0
    return speed


def _field_scalar(field: np.lib.npyio.NpzFile, key: str, lumen_mask: np.ndarray, *, default_value: float) -> np.ndarray:
    if key in field:
        scalar = np.asarray(field[key], dtype=float)
        if scalar.shape == lumen_mask.shape:
            scalar = scalar.copy()
            scalar[~lumen_mask] = np.nan
            return scalar
    scalar = np.full(lumen_mask.shape, float(default_value), dtype=float)
    scalar[~lumen_mask] = np.nan
    return scalar


def _field_local_shear_stress_pa(
    field: np.lib.npyio.NpzFile,
    velocity_xz_um_s: np.ndarray,
    lumen_mask: np.ndarray,
    spacing_um: float,
) -> tuple[np.ndarray, str]:
    """Load whole-lumen local shear, deriving a legacy diagnostic if absent."""

    lumen = np.asarray(lumen_mask, dtype=bool)
    if "local_shear_stress_pa" in field:
        values = np.asarray(field["local_shear_stress_pa"], dtype=float)
        if values.shape != lumen.shape:
            raise ValueError(
                "local_shear_stress_pa does not match the saved lumen grid."
            )
        values = values.copy()
        invalid = lumen & (~np.isfinite(values) | (values < 0.0))
        if np.any(invalid):
            raise ValueError(
                "local_shear_stress_pa contains invalid lumen values."
            )
        values[~lumen] = np.nan
        return values, "saved FEM local velocity gradient"

    fem_values = _legacy_fem_local_shear_stress_pa(field, lumen)
    if fem_values is not None:
        return fem_values, "reconstructed FEM local velocity gradient"

    # Legacy v18 fields saved the velocity and viscosity but not the local FEM
    # coefficients only in still older archives. A Cartesian derivative keeps
    # those result folders visualizable; new v19 fields save the authoritative
    # local-shear array directly.
    velocity = np.asarray(velocity_xz_um_s, dtype=float)
    if velocity.shape != (*lumen.shape, 2):
        raise ValueError("velocity_xz_um_s does not match the saved lumen grid.")
    step = float(spacing_um)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("Grid spacing must be finite and positive.")

    dux_dx, dux_dz = np.gradient(velocity[..., 0], step, step, edge_order=2)
    duz_dx, duz_dz = np.gradient(velocity[..., 1], step, step, edge_order=2)
    d_xz = 0.5 * (dux_dz + duz_dx)
    strain_double_contraction = dux_dx**2 + duz_dz**2 + 2.0 * d_xz**2
    equivalent_shear_rate = np.sqrt(
        np.maximum(2.0 * strain_double_contraction, 0.0)
    )

    if "viscosity_mpas" in field:
        viscosity_mpas = np.asarray(field["viscosity_mpas"], dtype=float)
    else:
        viscosity_mpas = np.full(lumen.shape, 3.0, dtype=float)
    if viscosity_mpas.shape != lumen.shape:
        raise ValueError("viscosity_mpas does not match the saved lumen grid.")
    positive = viscosity_mpas[
        lumen & np.isfinite(viscosity_mpas) & (viscosity_mpas > 0.0)
    ]
    default_viscosity_mpas = float(np.median(positive)) if positive.size else 3.0
    viscosity_pa_s = (
        np.where(
            np.isfinite(viscosity_mpas) & (viscosity_mpas > 0.0),
            viscosity_mpas,
            default_viscosity_mpas,
        )
        * 1.0e-3
    )
    values = viscosity_pa_s * equivalent_shear_rate
    values[~lumen] = np.nan
    return values, "legacy Cartesian velocity-gradient estimate"


def _legacy_fem_local_shear_stress_pa(
    field: np.lib.npyio.NpzFile,
    lumen_mask: np.ndarray,
) -> np.ndarray | None:
    """Reconstruct exact local FEM gradients from a legacy exported cache."""

    keys = (
        "fem_velocity_degree",
        "fem_cell_vertices_xz_um",
        "fem_polynomial_exponents",
        "fem_velocity_coefficients_um_s",
        "fem_bin_origin_xz_um",
        "fem_bin_size_um",
        "fem_bin_shape",
        "fem_bin_offsets",
        "fem_bin_cell_indices",
        "x_coordinates_um",
        "z_coordinates_um",
        "viscosity_mpas",
    )
    if any(key not in field for key in keys):
        return None

    finite_element = FiniteElementVelocityField(
        degree=int(np.asarray(field["fem_velocity_degree"]).reshape(-1)[0]),
        cell_vertices_xz_um=np.asarray(
            field["fem_cell_vertices_xz_um"], dtype=np.float64
        ),
        polynomial_exponents=np.asarray(
            field["fem_polynomial_exponents"], dtype=np.int16
        ),
        velocity_coefficients_um_s=np.asarray(
            field["fem_velocity_coefficients_um_s"], dtype=np.float64
        ),
        bin_origin_xz_um=np.asarray(field["fem_bin_origin_xz_um"], dtype=np.float64),
        bin_size_um=float(np.asarray(field["fem_bin_size_um"]).reshape(-1)[0]),
        bin_shape=tuple(
            int(value)
            for value in np.asarray(field["fem_bin_shape"]).reshape(-1)[:2]
        ),
        bin_offsets=np.asarray(field["fem_bin_offsets"], dtype=np.int64),
        bin_cell_indices=np.asarray(field["fem_bin_cell_indices"], dtype=np.int32),
    )
    lumen = np.asarray(lumen_mask, dtype=bool)
    indices = np.argwhere(lumen)
    values = np.full(lumen.shape, np.nan, dtype=float)
    if indices.size == 0:
        return values
    x = np.asarray(field["x_coordinates_um"], dtype=float)
    z = np.asarray(field["z_coordinates_um"], dtype=float)
    points = np.column_stack((x[indices[:, 0]], z[indices[:, 1]]))
    _, gradient, _ = sample_finite_element_velocity(
        finite_element,
        points,
        np.ones(points.shape[0], dtype=bool),
        use_numba=True,
    )

    strain_rate = 0.5 * (gradient + np.swapaxes(gradient, 1, 2))
    strain_invariant = np.einsum("nij,nij->n", strain_rate, strain_rate)
    equivalent_shear_rate = np.sqrt(np.maximum(2.0 * strain_invariant, 0.0))
    viscosity_mpas = np.asarray(field["viscosity_mpas"], dtype=float)[lumen]
    positive = viscosity_mpas[
        np.isfinite(viscosity_mpas) & (viscosity_mpas > 0.0)
    ]
    default_viscosity_mpas = float(np.median(positive)) if positive.size else 3.0
    viscosity_pa_s = (
        np.where(
            np.isfinite(viscosity_mpas) & (viscosity_mpas > 0.0),
            viscosity_mpas,
            default_viscosity_mpas,
        )
        * 1.0e-3
    )
    values[lumen] = viscosity_pa_s * equivalent_shear_rate
    return values


def _field_junction_mask(field: np.lib.npyio.NpzFile, lumen_mask: np.ndarray) -> np.ndarray:
    if "junction_core_mask" in field:
        junction = np.asarray(field["junction_core_mask"], dtype=bool)
        if junction.shape == lumen_mask.shape:
            return junction
    return np.zeros(lumen_mask.shape, dtype=bool)


def _field_open_boundary_mask(field: np.lib.npyio.NpzFile, lumen_mask: np.ndarray) -> np.ndarray:
    open_boundary = np.zeros(lumen_mask.shape, dtype=bool)
    if "inlet_label" in field:
        open_boundary |= np.asarray(field["inlet_label"] > 0, dtype=bool)
    if "outlet_label" in field:
        open_boundary |= np.asarray(field["outlet_label"] > 0, dtype=bool)
    return open_boundary & np.asarray(lumen_mask, dtype=bool)


def _hole_mask(lumen_mask: np.ndarray) -> np.ndarray:
    lumen = np.asarray(lumen_mask, dtype=bool)
    return ndimage.binary_fill_holes(lumen) & ~lumen
