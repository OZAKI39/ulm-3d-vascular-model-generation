"""Efficient non-adjacent branch collision checks before solid construction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from .config import CFDLumenConfig
from .types import BranchGeometry, CollisionEvent


@dataclass(frozen=True, slots=True)
class _Segment:
    branch_id: int
    segment_index: int
    start: np.ndarray
    end: np.ndarray
    radius_start: float
    radius_end: float
    midpoint: np.ndarray
    half_length: float
    maximum_radius: float


def _closest_points_on_segments(
    first_start: np.ndarray,
    first_end: np.ndarray,
    second_start: np.ndarray,
    second_end: np.ndarray,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Return closest segment parameters and points using a clamped analytic solution."""

    small = 1.0e-15
    u = first_end - first_start
    v = second_end - second_start
    w = first_start - second_start
    a = float(np.dot(u, u))
    b = float(np.dot(u, v))
    c = float(np.dot(v, v))
    d = float(np.dot(u, w))
    e = float(np.dot(v, w))
    denominator = a * c - b * b
    s_denominator = denominator
    t_denominator = denominator
    if denominator < small:
        s_numerator = 0.0
        s_denominator = 1.0
        t_numerator = e
        t_denominator = c
    else:
        s_numerator = b * e - c * d
        t_numerator = a * e - b * d
        if s_numerator < 0.0:
            s_numerator = 0.0
            t_numerator = e
            t_denominator = c
        elif s_numerator > s_denominator:
            s_numerator = s_denominator
            t_numerator = e + b
            t_denominator = c
    if t_numerator < 0.0:
        t_numerator = 0.0
        if -d < 0.0:
            s_numerator = 0.0
        elif -d > a:
            s_numerator = s_denominator
        else:
            s_numerator = -d
            s_denominator = a
    elif t_numerator > t_denominator:
        t_numerator = t_denominator
        if -d + b < 0.0:
            s_numerator = 0.0
        elif -d + b > a:
            s_numerator = s_denominator
        else:
            s_numerator = -d + b
            s_denominator = a
    first_parameter = 0.0 if abs(s_numerator) < small else s_numerator / s_denominator
    second_parameter = 0.0 if abs(t_numerator) < small else t_numerator / t_denominator
    first_point = first_start + first_parameter * u
    second_point = second_start + second_parameter * v
    return first_parameter, second_parameter, first_point, second_point


def _segments(branches: list[BranchGeometry]) -> list[_Segment]:
    output: list[_Segment] = []
    for branch in branches:
        for segment_index, (start, end) in enumerate(zip(branch.points_um[:-1], branch.points_um[1:])):
            length = float(np.linalg.norm(end - start))
            output.append(
                _Segment(
                    branch.branch_id,
                    segment_index,
                    start,
                    end,
                    float(branch.radius_um[segment_index]),
                    float(branch.radius_um[segment_index + 1]),
                    (start + end) * 0.5,
                    length * 0.5,
                    float(max(branch.radius_um[segment_index], branch.radius_um[segment_index + 1])),
                )
            )
    return output


def detect_nonadjacent_collisions(
    branches: list[BranchGeometry],
    config: CFDLumenConfig,
) -> tuple[list[CollisionEvent], dict[str, Any]]:
    if not config.collision_qc.enabled or len(branches) < 2:
        return [], {
            "enabled": config.collision_qc.enabled,
            "candidate_segment_pair_count": 0,
            "candidate_branch_pair_count": 0,
            "hard_collision_count": 0,
            "near_contact_count": 0,
            "status": "PASS",
        }
    segments = _segments(branches)
    midpoints = np.asarray([segment.midpoint for segment in segments], dtype=float)
    maximum_half_length = max(segment.half_length for segment in segments)
    maximum_radius = max(segment.maximum_radius for segment in segments)
    search_radius = (
        2.0 * maximum_half_length
        + 2.0 * maximum_radius
        + config.collision_qc.near_contact_tolerance_um
    )
    raw_pairs = cKDTree(midpoints).query_pairs(search_radius, output_type="ndarray")
    endpoints = {branch.branch_id: set(branch.local_node_ids[:: max(1, len(branch.local_node_ids) - 1)]) for branch in branches}
    minimum_by_branch_pair: dict[tuple[int, int], CollisionEvent] = {}
    exact_candidate_count = 0
    for first_index, second_index in raw_pairs:
        first = segments[int(first_index)]
        second = segments[int(second_index)]
        if first.branch_id == second.branch_id:
            continue
        branch_pair = tuple(sorted((first.branch_id, second.branch_id)))
        if endpoints[first.branch_id] & endpoints[second.branch_id]:
            continue
        midpoint_distance = float(np.linalg.norm(first.midpoint - second.midpoint))
        possible_distance = first.half_length + second.half_length
        possible_distance += first.maximum_radius + second.maximum_radius
        possible_distance += config.collision_qc.near_contact_tolerance_um
        if midpoint_distance > possible_distance:
            continue
        exact_candidate_count += 1
        first_t, second_t, first_point, second_point = _closest_points_on_segments(
            first.start, first.end, second.start, second.end
        )
        distance = float(np.linalg.norm(first_point - second_point))
        first_radius = first.radius_start + first_t * (first.radius_end - first.radius_start)
        second_radius = second.radius_start + second_t * (second.radius_end - second.radius_start)
        clearance = distance - first_radius - second_radius
        if clearance < -config.collision_qc.hard_collision_tolerance_um:
            classification = "HARD_COLLISION"
        elif clearance <= config.collision_qc.near_contact_tolerance_um:
            classification = "NEAR_CONTACT"
        else:
            classification = "CLEAR"
        event = CollisionEvent(
            branch_id_a=first.branch_id,
            branch_id_b=second.branch_id,
            segment_index_a=first.segment_index,
            segment_index_b=second.segment_index,
            distance_um=distance,
            radius_a_um=float(first_radius),
            radius_b_um=float(second_radius),
            clearance_um=float(clearance),
            classification=classification,
            closest_a_um=tuple(map(float, first_point)),
            closest_b_um=tuple(map(float, second_point)),
        )
        previous = minimum_by_branch_pair.get(branch_pair)
        if previous is None or event.clearance_um < previous.clearance_um:
            minimum_by_branch_pair[branch_pair] = event
    events = sorted(
        (event for event in minimum_by_branch_pair.values() if event.classification != "CLEAR"),
        key=lambda event: (event.clearance_um, event.branch_id_a, event.branch_id_b),
    )
    hard_count = sum(event.classification == "HARD_COLLISION" for event in events)
    near_count = sum(event.classification == "NEAR_CONTACT" for event in events)
    report = {
        "enabled": True,
        "candidate_segment_pair_count": exact_candidate_count,
        "candidate_branch_pair_count": len(minimum_by_branch_pair),
        "hard_collision_count": hard_count,
        "near_contact_count": near_count,
        "hard_collision_tolerance_um": config.collision_qc.hard_collision_tolerance_um,
        "near_contact_tolerance_um": config.collision_qc.near_contact_tolerance_um,
        "status": "FAIL" if hard_count else "PASS",
    }
    return events, report
