"""Create a local continuous-wall target at the earliest observable contact.

The target location is selected from an already generated, non-binding baseline
trajectory.  The first saved microbubble observation inside the configured
capture shell identifies the wall anchor; a contiguous arclength patch around
that anchor is then written in the v16 continuous-wall target format.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent
DEFAULT_BASELINE_RESULT = PACKAGE_DIR / "results" / "20260723_082654"
DEFAULT_OUTPUT = (
    PACKAGE_DIR
    / "study_inputs"
    / "fast_binding_visualization"
    / "continuous_wall_target.npz"
)
TARGET_SCHEMA = "v16_continuous_wall_arclength_target"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Place a local molecular target where a baseline microbubble first "
            "enters the requested capture shell."
        )
    )
    parser.add_argument(
        "--baseline-result",
        type=Path,
        default=DEFAULT_BASELINE_RESULT,
        help="Result directory containing the accepted field and baseline trajectory.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Destination v16 continuous-wall target NPZ.",
    )
    parser.add_argument(
        "--capture-distance-um",
        type=float,
        default=0.50,
        help="Capture-shell distance used to select the earliest observation.",
    )
    parser.add_argument(
        "--patch-half-length-um",
        type=float,
        default=120.0,
        help="Target-positive wall arclength on each side of the anchor.",
    )
    return parser.parse_args()


def _metadata_float(archive: np.lib.npyio.NpzFile, key: str, default: float) -> float:
    if "metadata_keys" not in archive or "metadata_values" not in archive:
        return float(default)
    keys = np.asarray(archive["metadata_keys"], dtype=str)
    values = np.asarray(archive["metadata_values"], dtype=str)
    match = np.flatnonzero(keys == key)
    return float(values[int(match[0])]) if match.size else float(default)


def _nearest_segment_index(
    point_xz: np.ndarray,
    starts_xz: np.ndarray,
    ends_xz: np.ndarray,
) -> tuple[int, np.ndarray, float]:
    segment = ends_xz - starts_xz
    length_squared = np.einsum("ij,ij->i", segment, segment)
    offset = point_xz[None, :] - starts_xz
    parameter = np.divide(
        np.einsum("ij,ij->i", offset, segment),
        length_squared,
        out=np.zeros(starts_xz.shape[0], dtype=float),
        where=length_squared > 0.0,
    )
    parameter = np.clip(parameter, 0.0, 1.0)
    projected = starts_xz + parameter[:, None] * segment
    distance_squared = np.einsum(
        "ij,ij->i", projected - point_xz[None, :], projected - point_xz[None, :]
    )
    index = int(np.argmin(distance_squared))
    return index, projected[index], float(np.sqrt(distance_squared[index]))


def prepare_fast_binding_target(
    baseline_result: Path,
    output: Path,
    *,
    capture_distance_um: float,
    patch_half_length_um: float,
) -> Path:
    baseline = Path(baseline_result).resolve()
    field_path = baseline / "velocity_and_wall_shear_field.npz"
    trajectory_path = baseline / "microbubble_field_trajectories.npz"
    if not field_path.is_file() or not trajectory_path.is_file():
        raise FileNotFoundError(
            "The baseline result must contain velocity_and_wall_shear_field.npz "
            "and microbubble_field_trajectories.npz."
        )
    if not np.isfinite(capture_distance_um) or capture_distance_um <= 0.0:
        raise ValueError("capture_distance_um must be finite and positive.")
    if not np.isfinite(patch_half_length_um) or patch_half_length_um <= 0.0:
        raise ValueError("patch_half_length_um must be finite and positive.")

    with np.load(trajectory_path, allow_pickle=False) as trajectory:
        required_trajectory = {
            "frame_offsets",
            "record_bubble_id",
            "record_positions_um",
            "record_diameter_um",
            "record_wall_gap_um",
        }
        missing = sorted(required_trajectory.difference(trajectory.files))
        if missing:
            raise KeyError(f"Baseline trajectory is missing arrays: {missing}")
        gaps_um = np.asarray(trajectory["record_wall_gap_um"], dtype=float)
        eligible = np.flatnonzero(gaps_um <= float(capture_distance_um))
        if eligible.size == 0:
            raise ValueError(
                "No baseline microbubble entered the requested capture shell; "
                f"minimum saved wall gap is {float(np.min(gaps_um)):.6g} um."
            )
        record_index = int(eligible[0])
        offsets = np.asarray(trajectory["frame_offsets"], dtype=np.int64)
        frame_index = int(
            np.searchsorted(offsets, record_index, side="right") - 1
        )
        bubble_id = int(trajectory["record_bubble_id"][record_index])
        centre_xz_um = np.asarray(
            trajectory["record_positions_um"][record_index, [0, 2]], dtype=float
        )
        bubble_radius_um = (
            0.5 * float(trajectory["record_diameter_um"][record_index])
        )
        wall_gap_um = float(gaps_um[record_index])
        output_dt_s = _metadata_float(trajectory, "dt_s", 1.0)

    with np.load(field_path, allow_pickle=False) as field:
        required_field = {
            "continuous_geometry_hash_sha256",
            "continuous_wall_start_xz_um",
            "continuous_wall_end_xz_um",
            "continuous_wall_ring_index",
            "continuous_wall_arclength_start_um",
            "continuous_wall_arclength_end_um",
            "continuous_boundary_ring_length_um",
        }
        missing = sorted(required_field.difference(field.files))
        if missing:
            raise KeyError(f"Accepted field is missing arrays: {missing}")
        geometry_hash = str(
            np.asarray(field["continuous_geometry_hash_sha256"]).item()
        )
        starts_xz_um = np.asarray(
            field["continuous_wall_start_xz_um"], dtype=float
        )
        ends_xz_um = np.asarray(field["continuous_wall_end_xz_um"], dtype=float)
        ring_index = np.asarray(
            field["continuous_wall_ring_index"], dtype=np.int64
        )
        arc_start_um = np.asarray(
            field["continuous_wall_arclength_start_um"], dtype=float
        )
        arc_end_um = np.asarray(
            field["continuous_wall_arclength_end_um"], dtype=float
        )
        ring_length_um = float(
            np.asarray(field["continuous_boundary_ring_length_um"]).item()
        )

    wall_index, wall_anchor_xz_um, centre_to_wall_um = _nearest_segment_index(
        centre_xz_um, starts_xz_um, ends_xz_um
    )
    anchor_arclength_um = 0.5 * (
        float(arc_start_um[wall_index]) + float(arc_end_um[wall_index])
    )
    centre_arclength_um = 0.5 * (arc_start_um + arc_end_um)
    raw_distance_um = np.abs(centre_arclength_um - anchor_arclength_um)
    cyclic_distance_um = np.minimum(
        raw_distance_um, ring_length_um - raw_distance_um
    )
    positive = np.ascontiguousarray(
        cyclic_distance_um <= float(patch_half_length_um), dtype=bool
    )
    if not np.any(positive):
        raise RuntimeError("Fast-target selection unexpectedly produced an empty patch.")

    wall_length_um = np.linalg.norm(ends_xz_um - starts_xz_um, axis=1)
    selected_length_um = float(np.sum(wall_length_um[positive]))
    payload = {
        "target_geometry_schema": np.asarray(TARGET_SCHEMA),
        "continuous_geometry_hash_sha256": np.asarray(geometry_hash),
        "selection_semantics": np.asarray(
            "Contiguous v16 wall-arclength patch centred on the earliest saved "
            "baseline observation inside the configured molecular capture shell."
        ),
        "baseline_result_dir": np.asarray(str(baseline)),
        "capture_distance_um": np.asarray(capture_distance_um, dtype=np.float64),
        "patch_half_length_um": np.asarray(
            patch_half_length_um, dtype=np.float64
        ),
        "selected_record_index": np.asarray(record_index, dtype=np.int64),
        "selected_frame_index": np.asarray(frame_index, dtype=np.int64),
        "selected_time_s": np.asarray(
            frame_index * output_dt_s, dtype=np.float64
        ),
        "selected_bubble_id": np.asarray(bubble_id, dtype=np.int64),
        "selected_bubble_centre_xz_um": centre_xz_um,
        "selected_bubble_radius_um": np.asarray(
            bubble_radius_um, dtype=np.float64
        ),
        "selected_wall_gap_um": np.asarray(wall_gap_um, dtype=np.float64),
        "selected_wall_index": np.asarray(wall_index, dtype=np.int64),
        "selected_wall_anchor_xz_um": wall_anchor_xz_um,
        "selected_centre_to_wall_um": np.asarray(
            centre_to_wall_um, dtype=np.float64
        ),
        "selected_wall_arclength_um": np.asarray(
            anchor_arclength_um, dtype=np.float64
        ),
        "selected_target_wall_length_um": np.asarray(
            selected_length_um, dtype=np.float64
        ),
        "wall_ring_index": ring_index,
        "wall_arclength_start_um": arc_start_um,
        "wall_arclength_end_um": arc_end_um,
        "boundary_ring_length_um": np.asarray(
            ring_length_um, dtype=np.float64
        ),
        "wall_start_xz_um": starts_xz_um,
        "wall_end_xz_um": ends_xz_um,
        "wall_center_xz_um": 0.5 * (starts_xz_um + ends_xz_um),
        "wall_length_um": wall_length_um,
        "wall_target_positive": positive,
    }
    destination = Path(output).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, **payload)
    return destination


def main() -> None:
    args = _parse_args()
    output = prepare_fast_binding_target(
        args.baseline_result,
        args.output,
        capture_distance_um=float(args.capture_distance_um),
        patch_half_length_um=float(args.patch_half_length_um),
    )
    with np.load(output, allow_pickle=False) as target:
        print(f"Saved target: {output}")
        print(
            "Earliest capture: "
            f"t={float(target['selected_time_s']):.3f} s, "
            f"bubble={int(target['selected_bubble_id'])}, "
            f"gap={float(target['selected_wall_gap_um']):.3f} um"
        )
        print(
            "Target wall length: "
            f"{float(target['selected_target_wall_length_um']):.1f} um"
        )


if __name__ == "__main__":
    main()
