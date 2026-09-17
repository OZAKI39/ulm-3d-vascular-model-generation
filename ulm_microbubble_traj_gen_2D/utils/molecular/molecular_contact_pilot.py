"""Contact-exposure analysis for dimensionless molecular-binding studies.

This module intentionally does not simulate, fit, infer, or optimise adhesion.
It only measures how long already generated trajectories visit a target-positive
reaction zone.  Intrinsic binding scenarios are defined independently, before
the run, so an observed contact duration can never be reused as a fitted
molecular association-rate parameter.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import yaml

from .molecular_binding_scenarios import (
    MOLECULES_PER_UM2_TO_MOLECULES_PER_M2,
    DaOnScenario,
)


@dataclass(frozen=True)
class MolecularContactEvent:
    """One contiguous sequence of target-positive saved records for one bubble."""

    bubble_id: int
    first_frame: int
    last_frame: int
    first_observed_time_s: float
    last_observed_time_s: float
    positive_record_count: int
    duration_s: float
    mean_reaction_area_um2: float
    mean_abs_tangential_slip_um_s: float
    right_censored: bool


@dataclass(frozen=True)
class BubbleContactSummary:
    """Cumulative target-contact opportunity for one permanent bubble ID."""

    bubble_id: int
    event_count: int
    positive_record_count: int
    cumulative_contact_time_s: float


@dataclass(frozen=True)
class MolecularContactPilotResult:
    """Immutable result of a contact-only trajectory analysis."""

    n_frames: int
    n_records: int
    n_unique_bubbles: int
    n_contacting_bubbles: int
    n_contact_records: int
    output_dt_s: float
    sampling_convention: str
    total_contact_time_s: float
    target_area_time_exposure_um2_s: float
    contacting_bubble_fraction: float
    median_positive_bubble_contact_time_s: float | None
    contact_weighted_mean_reaction_area_um2: float | None
    contact_weighted_mean_abs_tangential_slip_um_s: float | None
    right_censored_event_count: int
    right_censored_contact_time_s: float
    right_censored_contact_time_fraction: float
    numerical_wall_lock_contact_records: int
    numerical_wall_lock_contact_time_s: float
    exposure_summary: str
    events: tuple[MolecularContactEvent, ...]
    bubble_summaries: tuple[BubbleContactSummary, ...]


@dataclass(frozen=True, slots=True)
class _PreparedMolecularContactRecords:
    """Ratio-independent indexing for one immutable ragged trajectory.

    A capture-distance sweep used to rebuild the dense frame index and sort all
    saved records once per ratio.  Those operations depend only on the source
    trajectory, so the configured runner prepares them once and reuses them for
    every ratio.  The public analyzer still prepares its own instance, retaining
    its original standalone validation behaviour.
    """

    offsets: np.ndarray
    ids: np.ndarray
    output_dt_s: float
    n_frames: int
    frame_index: np.ndarray
    order_by_id_and_frame: np.ndarray
    unique_ids: np.ndarray
    group_starts: np.ndarray
    group_stops: np.ndarray


def analyze_molecular_contact_pilot(
    frame_offsets: np.ndarray,
    bubble_id: np.ndarray,
    reaction_area_um2: np.ndarray,
    tangential_slip_um_s: np.ndarray,
    output_dt_s: float,
    numerical_wall_lock: np.ndarray | None = None,
) -> MolecularContactPilotResult:
    """Measure target contact from a ragged saved-trajectory record stream.

    ``frame_offsets`` follows the trajectory NPZ convention: records belonging
    to frame ``k`` occupy ``[frame_offsets[k], frame_offsets[k + 1])``.  Contact
    exists wherever ``reaction_area_um2 > 0``.

    Contact time uses an explicit saved-frame rectangular approximation: every
    target-positive saved record represents one ``output_dt_s`` interval.  This
    makes the estimator deterministic for ragged records and easy to audit, but
    it remains limited by the output cadence.  A convergence study should use a
    smaller saved-frame interval when individual contacts last only one or two
    records.
    """

    prepared = _prepare_molecular_contact_records(
        frame_offsets,
        bubble_id,
        output_dt_s,
    )
    return _analyze_prepared_molecular_contact_pilot(
        prepared,
        reaction_area_um2,
        tangential_slip_um_s,
        numerical_wall_lock=numerical_wall_lock,
    )


def _prepare_molecular_contact_records(
    frame_offsets: np.ndarray,
    bubble_id: np.ndarray,
    output_dt_s: float,
) -> _PreparedMolecularContactRecords:
    """Validate and build the indexing shared by all capture-distance ratios."""

    offsets = _validate_frame_offsets(frame_offsets)
    ids = _validate_bubble_ids(bubble_id)
    dt_s = _finite_positive_float(output_dt_s, "output_dt_s")
    if int(offsets[-1]) != ids.size:
        raise ValueError(
            "frame_offsets[-1] must equal the number of ragged trajectory records."
        )

    n_frames = offsets.size - 1
    records_per_frame = np.diff(offsets)
    frame_index = np.repeat(np.arange(n_frames, dtype=np.int64), records_per_frame)
    if frame_index.size != ids.size:
        raise ValueError("frame_offsets are inconsistent with the record arrays.")

    if ids.size == 0:
        empty = np.empty(0, dtype=np.int64)
        return _PreparedMolecularContactRecords(
            offsets=offsets,
            ids=ids,
            output_dt_s=dt_s,
            n_frames=n_frames,
            frame_index=frame_index,
            order_by_id_and_frame=empty,
            unique_ids=empty,
            group_starts=empty,
            group_stops=empty,
        )

    # Sorting by permanent ID and then frame exposes invalid duplicate records
    # without assuming that the ragged arrays have already been grouped by ID.
    all_order = np.lexsort((frame_index, ids))
    sorted_ids = ids[all_order]
    sorted_frames = frame_index[all_order]
    duplicate = (sorted_ids[1:] == sorted_ids[:-1]) & (
        sorted_frames[1:] == sorted_frames[:-1]
    )
    if np.any(duplicate):
        duplicate_index = int(np.flatnonzero(duplicate)[0] + 1)
        raise ValueError(
            "Each bubble_id may occur at most once per frame; duplicate record "
            f"for bubble_id={int(sorted_ids[duplicate_index])}, "
            f"frame={int(sorted_frames[duplicate_index])}."
        )

    group_starts = np.r_[
        0,
        np.flatnonzero(sorted_ids[1:] != sorted_ids[:-1]) + 1,
    ].astype(np.int64, copy=False)
    group_stops = np.r_[group_starts[1:], sorted_ids.size].astype(
        np.int64,
        copy=False,
    )
    unique_ids = np.asarray(sorted_ids[group_starts], dtype=np.int64)
    return _PreparedMolecularContactRecords(
        offsets=offsets,
        ids=ids,
        output_dt_s=dt_s,
        n_frames=n_frames,
        frame_index=frame_index,
        order_by_id_and_frame=np.asarray(all_order, dtype=np.int64),
        unique_ids=unique_ids,
        group_starts=group_starts,
        group_stops=group_stops,
    )


def _analyze_prepared_molecular_contact_pilot(
    prepared: _PreparedMolecularContactRecords,
    reaction_area_um2: np.ndarray,
    tangential_slip_um_s: np.ndarray,
    *,
    numerical_wall_lock: np.ndarray | None = None,
) -> MolecularContactPilotResult:
    """Analyze one ratio using validated, ratio-independent trajectory indices."""

    ids = prepared.ids
    area = _validate_record_array(
        reaction_area_um2,
        "reaction_area_um2",
        expected_size=ids.size,
        non_negative=True,
    )
    slip = _validate_record_array(
        tangential_slip_um_s,
        "tangential_slip_um_s",
        expected_size=ids.size,
        non_negative=False,
    )
    if numerical_wall_lock is None:
        wall_lock = np.zeros(ids.size, dtype=bool)
    else:
        wall_lock = np.asarray(numerical_wall_lock, dtype=bool)
        if wall_lock.ndim != 1 or wall_lock.size != ids.size:
            raise ValueError(
                "numerical_wall_lock must be one-dimensional and match bubble_id length."
            )
    dt_s = prepared.output_dt_s
    n_frames = prepared.n_frames
    frame_index = prepared.frame_index

    if ids.size == 0:
        return _empty_pilot_result(n_frames=n_frames, output_dt_s=dt_s)

    contact = area > 0.0
    # This is deliberately a saved-frame rectangle rule rather than hidden
    # interpolation.  Its temporal resolution and boundary bias are written to
    # the report so scenario results cannot be mistaken for substep-resolved
    # contact kinetics.
    contact_time_weight_s = contact.astype(np.float64) * dt_s

    unique_ids = prepared.unique_ids
    contact_record_indices = np.flatnonzero(contact)
    events = _build_contact_events(
        ids=ids,
        frame_index=frame_index,
        area_um2=area,
        slip_um_s=slip,
        contact_time_weight_s=contact_time_weight_s,
        output_dt_s=dt_s,
        contact=contact,
        n_frames=n_frames,
        order_by_id_and_frame=prepared.order_by_id_and_frame,
    )
    bubble_summaries = _build_bubble_summaries(
        unique_ids=unique_ids,
        ids=ids,
        contact=contact,
        contact_time_weight_s=contact_time_weight_s,
        events=events,
        order_by_id_and_frame=prepared.order_by_id_and_frame,
        group_starts=prepared.group_starts,
        group_stops=prepared.group_stops,
    )

    positive_contact_times = np.asarray(
        [
            summary.cumulative_contact_time_s
            for summary in bubble_summaries
            if summary.cumulative_contact_time_s > 0.0
        ],
        dtype=np.float64,
    )
    total_contact_time_s = float(np.sum(contact_time_weight_s, dtype=np.float64))
    target_area_time_exposure_um2_s = float(
        np.dot(contact_time_weight_s, area)
    )
    numerical_lock_contact = contact & wall_lock
    numerical_wall_lock_contact_records = int(
        np.count_nonzero(numerical_lock_contact)
    )
    numerical_wall_lock_contact_time_s = float(
        np.sum(contact_time_weight_s[numerical_lock_contact], dtype=np.float64)
    )
    right_censored_events = tuple(event for event in events if event.right_censored)
    right_censored_contact_time_s = float(
        sum(event.duration_s for event in right_censored_events)
    )
    n_contacting_bubbles = sum(
        summary.positive_record_count > 0 for summary in bubble_summaries
    )
    if positive_contact_times.size:
        median_positive_bubble_contact_time_s: float | None = float(
            np.median(positive_contact_times)
        )
        contact_weighted_mean_reaction_area_um2: float | None = float(
            np.dot(contact_time_weight_s, area) / total_contact_time_s
        )
        contact_weighted_mean_abs_tangential_slip_um_s: float | None = float(
            np.dot(contact_time_weight_s, np.abs(slip)) / total_contact_time_s
        )
        if numerical_wall_lock_contact_records > 0:
            message = (
                "Positive target exposure was observed, but some exposed records "
                "overlap a numerical wall-lock diagnostic. Treat the exposure "
                "statistics as numerically suspect until that transport issue is removed."
            )
        else:
            message = (
                "Positive target exposure was observed. These values describe transport "
                "opportunity only and are not used to infer an intrinsic association rate."
            )
    else:
        median_positive_bubble_contact_time_s = None
        contact_weighted_mean_reaction_area_um2 = None
        contact_weighted_mean_abs_tangential_slip_um_s = None
        message = (
            "No target-positive reaction-area exposure was observed. This zero is a "
            "transport result; the predeclared molecular scenario table remains unchanged."
        )

    return MolecularContactPilotResult(
        n_frames=n_frames,
        n_records=ids.size,
        n_unique_bubbles=unique_ids.size,
        n_contacting_bubbles=n_contacting_bubbles,
        n_contact_records=contact_record_indices.size,
        output_dt_s=dt_s,
        sampling_convention="saved_frame_rectangle_each_positive_record_times_output_dt",
        total_contact_time_s=total_contact_time_s,
        target_area_time_exposure_um2_s=target_area_time_exposure_um2_s,
        contacting_bubble_fraction=(
            float(n_contacting_bubbles / unique_ids.size) if unique_ids.size else 0.0
        ),
        median_positive_bubble_contact_time_s=(
            median_positive_bubble_contact_time_s
        ),
        contact_weighted_mean_reaction_area_um2=(
            contact_weighted_mean_reaction_area_um2
        ),
        contact_weighted_mean_abs_tangential_slip_um_s=(
            contact_weighted_mean_abs_tangential_slip_um_s
        ),
        right_censored_event_count=len(right_censored_events),
        right_censored_contact_time_s=right_censored_contact_time_s,
        right_censored_contact_time_fraction=(
            right_censored_contact_time_s / total_contact_time_s
            if total_contact_time_s > 0.0
            else 0.0
        ),
        numerical_wall_lock_contact_records=numerical_wall_lock_contact_records,
        numerical_wall_lock_contact_time_s=numerical_wall_lock_contact_time_s,
        exposure_summary=message,
        events=events,
        bubble_summaries=bubble_summaries,
    )


def contact_pilot_report_mapping(
    pilot: MolecularContactPilotResult,
    scenarios: Sequence[DaOnScenario] = (),
    *,
    study_context: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a YAML-safe, human-readable report without performing I/O."""

    report = {
        "report_kind": "molecular_contact_exposure_and_predeclared_da_on_sweep",
        "interpretation": (
            "sensitivity study only; intrinsic association scenarios were fixed before "
            "transport and were not inferred from observed contact"
        ),
        "density_conversion": {
            "from": "molecule/um^2",
            "to": "molecule/m^2",
            "multiply_by": MOLECULES_PER_UM2_TO_MOLECULES_PER_M2,
        },
        "pilot": {
            "n_frames": pilot.n_frames,
            "n_records": pilot.n_records,
            "n_unique_bubbles": pilot.n_unique_bubbles,
            "n_contacting_bubbles": pilot.n_contacting_bubbles,
            "n_contact_records": pilot.n_contact_records,
            "output_dt_s": pilot.output_dt_s,
            "sampling_convention": pilot.sampling_convention,
            "total_contact_time_s": pilot.total_contact_time_s,
            "target_area_time_exposure_um2_s": (
                pilot.target_area_time_exposure_um2_s
            ),
            "molecular_comparison_ready": (
                pilot.target_area_time_exposure_um2_s > 0.0
            ),
            "contacting_bubble_fraction": pilot.contacting_bubble_fraction,
            "median_positive_bubble_contact_time_s": (
                pilot.median_positive_bubble_contact_time_s
            ),
            "contact_weighted_mean_reaction_area_um2": (
                pilot.contact_weighted_mean_reaction_area_um2
            ),
            "contact_weighted_mean_abs_tangential_slip_um_s": (
                pilot.contact_weighted_mean_abs_tangential_slip_um_s
            ),
            "right_censored_event_count": pilot.right_censored_event_count,
            "right_censored_contact_time_s": pilot.right_censored_contact_time_s,
            "right_censored_contact_time_fraction": (
                pilot.right_censored_contact_time_fraction
            ),
            "numerical_wall_lock_contact_records": (
                pilot.numerical_wall_lock_contact_records
            ),
            "numerical_wall_lock_contact_time_s": (
                pilot.numerical_wall_lock_contact_time_s
            ),
            "exposure_summary": pilot.exposure_summary,
        },
        "per_bubble_contact": [
            {
                "bubble_id": summary.bubble_id,
                "event_count": summary.event_count,
                "positive_record_count": summary.positive_record_count,
                "cumulative_contact_time_s": summary.cumulative_contact_time_s,
            }
            for summary in pilot.bubble_summaries
        ],
        "contact_events": [
            {
                "bubble_id": event.bubble_id,
                "first_frame": event.first_frame,
                "last_frame": event.last_frame,
                "first_observed_time_s": event.first_observed_time_s,
                "last_observed_time_s": event.last_observed_time_s,
                "positive_record_count": event.positive_record_count,
                "duration_s": event.duration_s,
                "mean_reaction_area_um2": event.mean_reaction_area_um2,
                "mean_abs_tangential_slip_um_s": (
                    event.mean_abs_tangential_slip_um_s
                ),
                "right_censored": event.right_censored,
            }
            for event in pilot.events
        ],
        "da_on_scenarios": [
            {
                "scenario_index": scenario.scenario_index,
                "da_on": scenario.da_on,
                "target_density_molecules_per_um2": (
                    scenario.target_density_molecules_per_um2
                ),
                "target_density_molecules_per_m2": (
                    scenario.target_density_molecules_per_m2
                ),
                "ligand_density_molecules_per_um2": (
                    scenario.ligand_density_molecules_per_um2
                ),
                "ligand_density_molecules_per_m2": (
                    scenario.ligand_density_molecules_per_m2
                ),
                "capture_distance_to_rest_length_ratio": (
                    scenario.capture_distance_to_rest_length_ratio
                ),
                "rest_length_um": scenario.rest_length_um,
                "capture_distance_um": scenario.capture_distance_um,
                "da_on_reference_time_s": scenario.da_on_reference_time_s,
                "association_rate_m2_per_molecule_s": (
                    scenario.association_rate_m2_per_molecule_s
                ),
            }
            for scenario in scenarios
        ],
    }
    if study_context is not None:
        report["study_context"] = dict(study_context)
    return report


def render_contact_pilot_yaml(
    pilot: MolecularContactPilotResult,
    scenarios: Sequence[DaOnScenario] = (),
    *,
    study_context: Mapping[str, Any] | None = None,
) -> str:
    """Render a report as YAML without mutating the filesystem."""

    return yaml.safe_dump(
        contact_pilot_report_mapping(
            pilot,
            scenarios,
            study_context=study_context,
        ),
        sort_keys=False,
        allow_unicode=True,
    )


def save_contact_pilot_yaml(
    path: str | Path,
    pilot: MolecularContactPilotResult,
    scenarios: Sequence[DaOnScenario] = (),
    *,
    study_context: Mapping[str, Any] | None = None,
) -> Path:
    """Persist the human-readable YAML report and return its resolved path."""

    output_path = Path(path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        render_contact_pilot_yaml(
            pilot,
            scenarios,
            study_context=study_context,
        ),
        encoding="utf-8",
    )
    return output_path


def _build_contact_events(
    *,
    ids: np.ndarray,
    frame_index: np.ndarray,
    area_um2: np.ndarray,
    slip_um_s: np.ndarray,
    contact_time_weight_s: np.ndarray,
    output_dt_s: float,
    contact: np.ndarray,
    n_frames: int,
    order_by_id_and_frame: np.ndarray,
) -> tuple[MolecularContactEvent, ...]:
    if not np.any(contact):
        return ()
    # ``order_by_id_and_frame`` was already validated for the whole trajectory.
    # Selecting its contact-positive rows is equivalent to sorting the contact
    # subset by the same keys, without repeating an O(N log N) sort per ratio.
    records = order_by_id_and_frame[np.asarray(contact[order_by_id_and_frame], dtype=bool)]
    ordered_ids = ids[records]
    ordered_frames = frame_index[records]
    begins_event = np.ones(records.size, dtype=bool)
    begins_event[1:] = (ordered_ids[1:] != ordered_ids[:-1]) | (
        ordered_frames[1:] != ordered_frames[:-1] + 1
    )
    starts = np.flatnonzero(begins_event)
    stops = np.r_[starts[1:], records.size]
    events: list[MolecularContactEvent] = []
    for start, stop in zip(starts, stops):
        event_records = records[start:stop]
        duration_s = float(np.sum(contact_time_weight_s[event_records], dtype=np.float64))
        if duration_s > 0.0:
            mean_area_um2 = float(
                np.dot(
                    contact_time_weight_s[event_records], area_um2[event_records]
                )
                / duration_s
            )
            mean_abs_slip_um_s = float(
                np.dot(
                    contact_time_weight_s[event_records],
                    np.abs(slip_um_s[event_records]),
                )
                / duration_s
            )
        else:
            mean_area_um2 = float(np.mean(area_um2[event_records]))
            mean_abs_slip_um_s = float(np.mean(np.abs(slip_um_s[event_records])))
        first_frame = int(ordered_frames[start])
        last_frame = int(ordered_frames[stop - 1])
        events.append(
            MolecularContactEvent(
                bubble_id=int(ordered_ids[start]),
                first_frame=first_frame,
                last_frame=last_frame,
                first_observed_time_s=first_frame * output_dt_s,
                last_observed_time_s=last_frame * output_dt_s,
                positive_record_count=int(stop - start),
                duration_s=duration_s,
                mean_reaction_area_um2=mean_area_um2,
                mean_abs_tangential_slip_um_s=mean_abs_slip_um_s,
                right_censored=last_frame == n_frames - 1,
            )
        )
    events.sort(key=lambda event: (event.first_frame, event.bubble_id, event.last_frame))
    return tuple(events)


def _build_bubble_summaries(
    *,
    unique_ids: np.ndarray,
    ids: np.ndarray,
    contact: np.ndarray,
    contact_time_weight_s: np.ndarray,
    events: tuple[MolecularContactEvent, ...],
    order_by_id_and_frame: np.ndarray,
    group_starts: np.ndarray,
    group_stops: np.ndarray,
) -> tuple[BubbleContactSummary, ...]:
    event_count_by_id: dict[int, int] = {}
    for event in events:
        event_count_by_id[event.bubble_id] = event_count_by_id.get(event.bubble_id, 0) + 1
    summaries: list[BubbleContactSummary] = []
    for raw_id, start, stop in zip(unique_ids, group_starts, group_stops):
        bubble = int(raw_id)
        bubble_records = order_by_id_and_frame[start:stop]
        summaries.append(
            BubbleContactSummary(
                bubble_id=bubble,
                event_count=event_count_by_id.get(bubble, 0),
                positive_record_count=int(np.count_nonzero(contact[bubble_records])),
                cumulative_contact_time_s=float(
                    np.sum(contact_time_weight_s[bubble_records], dtype=np.float64)
                ),
            )
        )
    return tuple(summaries)


def _empty_pilot_result(*, n_frames: int, output_dt_s: float) -> MolecularContactPilotResult:
    message = (
        "The trajectory contains no bubble records, so observed target exposure is zero. "
        "The predeclared molecular scenario table remains independent of this outcome."
    )
    return MolecularContactPilotResult(
        n_frames=n_frames,
        n_records=0,
        n_unique_bubbles=0,
        n_contacting_bubbles=0,
        n_contact_records=0,
        output_dt_s=output_dt_s,
        sampling_convention="saved_frame_rectangle_each_positive_record_times_output_dt",
        total_contact_time_s=0.0,
        target_area_time_exposure_um2_s=0.0,
        contacting_bubble_fraction=0.0,
        median_positive_bubble_contact_time_s=None,
        contact_weighted_mean_reaction_area_um2=None,
        contact_weighted_mean_abs_tangential_slip_um_s=None,
        right_censored_event_count=0,
        right_censored_contact_time_s=0.0,
        right_censored_contact_time_fraction=0.0,
        numerical_wall_lock_contact_records=0,
        numerical_wall_lock_contact_time_s=0.0,
        exposure_summary=message,
        events=(),
        bubble_summaries=(),
    )


def _validate_frame_offsets(value: np.ndarray) -> np.ndarray:
    offsets = np.asarray(value)
    if offsets.ndim != 1 or offsets.size < 2:
        raise ValueError("frame_offsets must be a one-dimensional array with at least two entries.")
    if offsets.dtype.kind not in "iu" or offsets.dtype.kind == "b":
        raise ValueError("frame_offsets must contain integers.")
    offsets = np.asarray(offsets, dtype=np.int64)
    if offsets[0] != 0:
        raise ValueError("frame_offsets must start at zero.")
    if np.any(np.diff(offsets) < 0):
        raise ValueError("frame_offsets must be non-decreasing.")
    return offsets


def _validate_bubble_ids(value: np.ndarray) -> np.ndarray:
    ids = np.asarray(value)
    if ids.ndim != 1:
        raise ValueError("bubble_id must be a one-dimensional array.")
    if ids.dtype.kind not in "iu" or ids.dtype.kind == "b":
        raise ValueError("bubble_id must contain permanent integer IDs.")
    ids = np.asarray(ids, dtype=np.int64)
    if np.any(ids < 0):
        raise ValueError("bubble_id must contain non-negative permanent IDs.")
    return ids


def _validate_record_array(
    value: np.ndarray,
    name: str,
    *,
    expected_size: int,
    non_negative: bool,
) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.ndim != 1 or array.size != expected_size:
        raise ValueError(f"{name} must be one-dimensional and match bubble_id length.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values.")
    if non_negative and np.any(array < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    return array


def _finite_positive_float(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be finite and positive.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite and positive.") from exc
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")
    return number
