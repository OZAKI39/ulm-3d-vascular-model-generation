"""Source-constrained C1 centerlines and adaptive v9 discretization."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.spatial import cKDTree

from .config import CFDLumenConfig
from .types import BranchGeometry, GeometryValidationError, PortGeometry
from .unified_polyball import _copy_with_tail


@dataclass(slots=True)
class SmoothCenterlineBuild:
    branches: list[BranchGeometry]
    port_tail_rows: list[dict[str, Any]]
    selected_tangent_angle_deg: float
    selected_sagitta_radius_fraction: float
    sensitivity_rows: list[dict[str, Any]]
    branch_fidelity_rows: list[dict[str, Any]]
    report: dict[str, Any]


def _summary(values: np.ndarray) -> dict[str, float | int | None]:
    data = np.asarray(values, dtype=float)
    if not len(data):
        return {"count": 0, "mean": None, "p95": None, "p99": None, "max": None}
    return {
        "count": int(len(data)),
        "mean": float(np.mean(data)),
        "p95": float(np.percentile(data, 95)),
        "p99": float(np.percentile(data, 99)),
        "max": float(np.max(data)),
    }


def _unit_rows(vectors: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=float)
    norm = np.linalg.norm(vectors, axis=1)
    return np.divide(
        vectors,
        norm[:, None],
        out=np.zeros_like(vectors),
        where=norm[:, None] > 1.0e-14,
    )


def _angles(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    first_unit = _unit_rows(first)
    second_unit = _unit_rows(second)
    cosine = np.einsum("ij,ij->i", first_unit, second_unit)
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def _centripetal_parameter(points: np.ndarray) -> np.ndarray:
    chord = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if np.any(chord <= 1.0e-12):
        raise GeometryValidationError("V9 spline controls contain a zero-length chord")
    return np.concatenate(([0.0], np.cumsum(np.sqrt(chord))))


def _centripetal_tangents(points: np.ndarray, parameter: np.ndarray) -> np.ndarray:
    slopes = np.diff(points, axis=0) / np.diff(parameter)[:, None]
    tangent = np.empty_like(points)
    tangent[0] = slopes[0]
    tangent[-1] = slopes[-1]
    for index in range(1, len(points) - 1):
        before = parameter[index] - parameter[index - 1]
        after = parameter[index + 1] - parameter[index]
        tangent[index] = (after * slopes[index - 1] + before * slopes[index]) / (
            before + after
        )
    return tangent


def _distance_to_segments(points: np.ndarray, polyline: np.ndarray) -> np.ndarray:
    query = np.asarray(points, dtype=float).reshape((-1, 3))
    start = np.asarray(polyline[:-1], dtype=float)
    vector = np.diff(polyline, axis=0)
    squared = np.einsum("ij,ij->i", vector, vector)
    output = np.full(len(query), np.inf, dtype=float)
    for offset in range(0, len(query), 10_000):
        selected = query[offset : offset + 10_000]
        relative = selected[:, None, :] - start[None, :, :]
        t = np.clip(
            np.divide(
                np.einsum("nsi,si->ns", relative, vector),
                squared[None, :],
                out=np.zeros((len(selected), len(vector)), dtype=float),
                where=squared[None, :] > 1.0e-20,
            ),
            0.0,
            1.0,
        )
        closest = start[None, :, :] + t[:, :, None] * vector[None, :, :]
        output[offset : offset + len(selected)] = np.min(
            np.linalg.norm(selected[:, None, :] - closest, axis=2), axis=1
        )
    return output


def _segment_distance(
    first_start: np.ndarray,
    first_end: np.ndarray,
    second_start: np.ndarray,
    second_end: np.ndarray,
) -> float:
    """Shortest distance between two finite 3-D segments."""

    u = first_end - first_start
    v = second_end - second_start
    w = first_start - second_start
    a = float(np.dot(u, u))
    b = float(np.dot(u, v))
    c = float(np.dot(v, v))
    d = float(np.dot(u, w))
    e = float(np.dot(v, w))
    denominator = a * c - b * b
    small = 1.0e-20
    if denominator < small:
        s_numerator, s_denominator = 0.0, 1.0
        t_numerator, t_denominator = e, c
    else:
        s_numerator = b * e - c * d
        t_numerator = a * e - b * d
        s_denominator = t_denominator = denominator
        if s_numerator < 0.0:
            s_numerator = 0.0
            t_numerator, t_denominator = e, c
        elif s_numerator > s_denominator:
            s_numerator = s_denominator
            t_numerator, t_denominator = e + b, c
    if t_numerator < 0.0:
        t_numerator = 0.0
        if -d < 0.0:
            s_numerator = 0.0
        elif -d > a:
            s_numerator = s_denominator
        else:
            s_numerator, s_denominator = -d, a
    elif t_numerator > t_denominator:
        t_numerator = t_denominator
        if -d + b < 0.0:
            s_numerator = 0.0
        elif -d + b > a:
            s_numerator = s_denominator
        else:
            s_numerator, s_denominator = -d + b, a
    s = 0.0 if abs(s_numerator) < small else s_numerator / s_denominator
    t = 0.0 if abs(t_numerator) < small else t_numerator / t_denominator
    return float(np.linalg.norm(w + s * u - t * v))


def _self_intersection_count(points: np.ndarray) -> int:
    starts = np.asarray(points[:-1], dtype=float)
    ends = np.asarray(points[1:], dtype=float)
    lengths = np.linalg.norm(ends - starts, axis=1)
    if len(starts) < 4:
        return 0
    midpoint = 0.5 * (starts + ends)
    candidate_pairs = cKDTree(midpoint).query_pairs(
        r=float(np.max(lengths) + 1.0e-8), output_type="ndarray"
    )
    count = 0
    for first, second in np.asarray(candidate_pairs, dtype=np.int64).reshape((-1, 2)):
        if abs(int(first) - int(second)) <= 2:
            continue
        if _segment_distance(starts[first], ends[first], starts[second], ends[second]) <= 1.0e-8:
            count += 1
    return count


def _raw_control_indices(control: np.ndarray, raw: np.ndarray) -> np.ndarray:
    indices: list[int] = []
    cursor = 0
    for point in raw:
        matches = np.flatnonzero(
            np.all(np.isclose(control[cursor:], point[None, :], rtol=0.0, atol=1.0e-12), axis=1)
        )
        if not len(matches):
            raise GeometryValidationError(
                "V9 source SWC point is absent from spline interpolation controls"
            )
        selected = cursor + int(matches[0])
        indices.append(selected)
        cursor = selected + 1
    return np.asarray(indices, dtype=np.int64)


def _adaptive_parameters(
    curve: CubicHermiteSpline,
    radius: PchipInterpolator,
    knots: np.ndarray,
    *,
    alpha: float,
    max_spacing_um: float,
    theta_max_deg: float,
    eta: float,
) -> np.ndarray:
    output: list[float] = [float(knots[0])]

    def append_interval(first: float, second: float, depth: int) -> None:
        midpoint = 0.5 * (first + second)
        samples = np.asarray((first, midpoint, second), dtype=float)
        point = np.asarray(curve(samples), dtype=float)
        tangent = np.asarray(curve(samples, 1), dtype=float)
        local_radius = np.asarray(radius(samples), dtype=float)
        chord = float(np.linalg.norm(point[2] - point[0]))
        tangent_angle = float(_angles(tangent[[0]], tangent[[2]])[0])
        chord_vector = point[2] - point[0]
        chord_squared = float(np.dot(chord_vector, chord_vector))
        fraction = (
            float(np.clip(np.dot(point[1] - point[0], chord_vector) / chord_squared, 0.0, 1.0))
            if chord_squared > 1.0e-20
            else 0.0
        )
        closest = point[0] + fraction * chord_vector
        sagitta = float(np.linalg.norm(point[1] - closest))
        allowed_spacing = min(max_spacing_um, alpha * float(np.min(local_radius)))
        acceptable = (
            chord <= allowed_spacing * (1.0 + 1.0e-12)
            and tangent_angle <= theta_max_deg * (1.0 + 1.0e-12)
            and sagitta <= eta * float(local_radius[1]) * (1.0 + 1.0e-12)
        )
        if acceptable:
            output.append(float(second))
            return
        if depth >= 28 or second - first <= 1.0e-12:
            raise GeometryValidationError("V9 adaptive spline discretization did not converge")
        append_interval(first, midpoint, depth + 1)
        append_interval(midpoint, second, depth + 1)

    for first, second in zip(knots[:-1], knots[1:]):
        append_interval(float(first), float(second), 0)
    parameter = np.asarray(output, dtype=float)
    for _ in range(12):
        point = np.asarray(curve(parameter), dtype=float)
        chord = np.diff(point, axis=0)
        turns = _angles(chord[:-1], chord[1:]) if len(chord) > 1 else np.empty(0)
        offenders = np.flatnonzero(turns > theta_max_deg * (1.0 + 1.0e-10)) + 1
        if not len(offenders):
            return parameter
        additions: set[float] = set()
        for index in offenders:
            additions.add(float(0.5 * (parameter[index - 1] + parameter[index])))
            additions.add(float(0.5 * (parameter[index] + parameter[index + 1])))
        parameter = np.unique(np.concatenate((parameter, np.asarray(tuple(additions)))))
    raise GeometryValidationError("V9 polyline tangent-angle refinement did not converge")


def _derive_branch(
    source: BranchGeometry,
    control_branch: BranchGeometry,
    config: CFDLumenConfig,
    *,
    theta_max_deg: float,
    eta: float,
) -> tuple[BranchGeometry, dict[str, Any]]:
    control = np.asarray(control_branch.points_um, dtype=float)
    control_radius = np.asarray(control_branch.radius_um, dtype=float)
    parameter = _centripetal_parameter(control)
    tangent = _centripetal_tangents(control, parameter)
    curve = CubicHermiteSpline(parameter, control, tangent, axis=0, extrapolate=False)
    radius = PchipInterpolator(parameter, control_radius, extrapolate=False)
    sample = _adaptive_parameters(
        curve,
        radius,
        parameter,
        alpha=config.v9.spacing_radius_fraction,
        max_spacing_um=config.v9.max_spacing_um,
        theta_max_deg=theta_max_deg,
        eta=eta,
    )
    points = np.asarray(curve(sample), dtype=float)
    radii = np.asarray(radius(sample), dtype=float)
    # Restore every interpolation constraint bit-for-bit in the discretized output.
    for knot_id, knot in enumerate(parameter):
        sample_id = int(np.flatnonzero(sample == knot)[0])
        points[sample_id] = control[knot_id]
        radii[sample_id] = control_radius[knot_id]
    if np.any(~np.isfinite(radii)) or np.any(radii <= 0.0):
        raise GeometryValidationError("V9 PCHIP radius field is non-finite/non-positive")
    arc = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))))
    output = replace(control_branch)
    output.points_um = points
    output.radius_um = radii
    output.arc_length_um = arc

    raw = np.asarray(source.raw_points_um, dtype=float)
    raw_indices = _raw_control_indices(control, raw)
    first_parameter = parameter[raw_indices[0]]
    last_parameter = parameter[raw_indices[-1]]
    source_mask = (sample >= first_parameter) & (sample <= last_parameter)
    source_spline = points[source_mask]
    source_to_spline = _distance_to_segments(raw, source_spline)
    spline_to_source = _distance_to_segments(source_spline, raw)
    combined = np.concatenate((source_to_spline, spline_to_source))
    source_radius_min = float(np.min(source.raw_radius_um))
    raw_length = float(np.linalg.norm(np.diff(raw, axis=0), axis=1).sum())
    spline_length = float(
        np.linalg.norm(np.diff(source_spline, axis=0), axis=1).sum()
    )
    chords = np.diff(points, axis=0)
    turn = _angles(chords[:-1], chords[1:]) if len(chords) > 1 else np.empty(0)
    midpoint_parameter = 0.5 * (sample[:-1] + sample[1:])
    midpoint = np.asarray(curve(midpoint_parameter), dtype=float)
    chord_start = points[:-1]
    chord_vector = np.diff(points, axis=0)
    chord_squared = np.einsum("ij,ij->i", chord_vector, chord_vector)
    chord_fraction = np.divide(
        np.einsum("ij,ij->i", midpoint - chord_start, chord_vector),
        chord_squared,
        out=np.zeros(len(chord_squared), dtype=float),
        where=chord_squared > 1.0e-20,
    )
    chord_fraction = np.clip(chord_fraction, 0.0, 1.0)
    chord_closest = chord_start + chord_fraction[:, None] * chord_vector
    sagitta = np.linalg.norm(midpoint - chord_closest, axis=1)
    midpoint_radius = np.asarray(radius(midpoint_parameter), dtype=float)
    self_intersections = _self_intersection_count(points)
    raw_sample_ids = np.asarray(
        [int(np.flatnonzero(sample == parameter[index])[0]) for index in raw_indices],
        dtype=np.int64,
    )
    endpoint_start_error = float(np.linalg.norm(points[raw_sample_ids[0]] - raw[0]))
    endpoint_end_error = float(np.linalg.norm(points[raw_sample_ids[-1]] - raw[-1]))
    row = {
        "branch_id": int(source.branch_id),
        "method": "CFD_DERIVED_SPLINE_CENTERLINE",
        "interpolator": "centripetal cubic Hermite",
        "radius_interpolator": "PCHIP",
        "raw_source_point_count": int(len(raw)),
        "control_point_count": int(len(control)),
        "discretized_point_count": int(len(points)),
        "source_points_retained_as_exact_constraints": bool(
            np.all(source_to_spline <= 1.0e-12)
        ),
        "source_point_to_spline_distance_um": _summary(source_to_spline),
        "spline_to_source_polyline_distance_um": _summary(spline_to_source),
        "hausdorff_distance_um": float(np.max(combined)),
        "p95_distance_um": float(np.percentile(combined, 95)),
        "hausdorff_distance_min_radius_fraction": float(np.max(combined) / source_radius_min),
        "p95_distance_min_radius_fraction": float(np.percentile(combined, 95) / source_radius_min),
        "source_polyline_length_um": raw_length,
        "spline_centerline_length_um": spline_length,
        "branch_length_change_fraction": float((spline_length - raw_length) / raw_length),
        "junction_position_error_um": max(endpoint_start_error, endpoint_end_error),
        "endpoint_position_error_um": max(endpoint_start_error, endpoint_end_error),
        "continuous_c1_tangent_jump_deg": 0.0,
        "discretized_tangent_angle_deg": _summary(turn),
        "maximum_spacing_um": float(np.max(np.linalg.norm(chords, axis=1))),
        "maximum_sagitta_radius_fraction": float(np.max(sagitta / midpoint_radius)),
        "self_intersection_count": self_intersections,
    }
    row["checks"] = {
        "source_constraints_exact": row["source_points_retained_as_exact_constraints"],
        "junction_position_exact": row["junction_position_error_um"] == 0.0,
        "endpoint_position_exact": row["endpoint_position_error_um"] == 0.0,
        "no_curve_self_intersection": self_intersections == 0,
        "tangent_angle_within_limit": (row["discretized_tangent_angle_deg"]["max"] or 0.0)
        <= theta_max_deg * (1.0 + 1.0e-9),
        "sagitta_within_limit": row["maximum_sagitta_radius_fraction"]
        <= eta * (1.0 + 1.0e-9),
        "source_p95_within_limit": row["p95_distance_min_radius_fraction"]
        <= config.v9.maximum_source_p95_radius_fraction,
        "source_hausdorff_within_limit": row["hausdorff_distance_min_radius_fraction"]
        <= config.v9.maximum_source_hausdorff_radius_fraction,
        "branch_length_change_within_limit": abs(row["branch_length_change_fraction"])
        <= config.v9.maximum_branch_length_change_fraction,
    }
    row["status"] = "PASS" if all(row["checks"].values()) else "FAIL"
    return output, row


def build_smooth_centerline(
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> SmoothCenterlineBuild:
    """Evaluate all v9 adaptive pairs, then select by measured fidelity/turning."""

    raw_controls: list[BranchGeometry] = []
    for branch in branches:
        copied = replace(branch)
        copied.points_um = np.asarray(branch.raw_points_um, dtype=float).copy()
        copied.radius_um = np.asarray(branch.raw_radius_um, dtype=float).copy()
        copied.arc_length_um = np.concatenate(
            ([0.0], np.cumsum(np.linalg.norm(np.diff(copied.points_um, axis=0), axis=1)))
        )
        raw_controls.append(copied)
    controls, tail_rows = _copy_with_tail(raw_controls, ports, config)
    source_by_id = {branch.branch_id: branch for branch in branches}
    candidates: dict[tuple[float, float], tuple[list[BranchGeometry], list[dict[str, Any]]]] = {}
    sensitivity: list[dict[str, Any]] = []
    for theta in config.v9.adaptive_tangent_angle_degrees:
        for eta in config.v9.adaptive_sagitta_radius_fractions:
            derived: list[BranchGeometry] = []
            rows: list[dict[str, Any]] = []
            for control in controls:
                branch, row = _derive_branch(
                    source_by_id[control.branch_id],
                    control,
                    config,
                    theta_max_deg=float(theta),
                    eta=float(eta),
                )
                derived.append(branch)
                rows.append(row)
            aggregate = {
                "theta_max_deg": float(theta),
                "eta_radius_fraction": float(eta),
                "point_count": int(sum(len(branch.points_um) for branch in derived)),
                "maximum_discretized_tangent_angle_deg": max(
                    float(row["discretized_tangent_angle_deg"]["max"] or 0.0)
                    for row in rows
                ),
                "maximum_sagitta_radius_fraction": max(
                    float(row["maximum_sagitta_radius_fraction"]) for row in rows
                ),
                "maximum_source_hausdorff_radius_fraction": max(
                    float(row["hausdorff_distance_min_radius_fraction"]) for row in rows
                ),
                "maximum_source_p95_radius_fraction": max(
                    float(row["p95_distance_min_radius_fraction"]) for row in rows
                ),
                "maximum_absolute_branch_length_change_fraction": max(
                    abs(float(row["branch_length_change_fraction"])) for row in rows
                ),
                "junction_position_error_um": max(
                    float(row["junction_position_error_um"]) for row in rows
                ),
                "endpoint_position_error_um": max(
                    float(row["endpoint_position_error_um"]) for row in rows
                ),
                "self_intersection_count": int(
                    sum(int(row["self_intersection_count"]) for row in rows)
                ),
                "status": "PASS" if all(row["status"] == "PASS" for row in rows) else "FAIL",
            }
            sensitivity.append(aggregate)
            candidates[(float(theta), float(eta))] = (derived, rows)
    eligible = [row for row in sensitivity if row["status"] == "PASS"]
    if not eligible:
        raise GeometryValidationError("V9 no adaptive spline discretization passed fidelity QC")
    selected = min(
        eligible,
        key=lambda row: (
            row["maximum_discretized_tangent_angle_deg"],
            row["maximum_source_hausdorff_radius_fraction"],
            row["maximum_sagitta_radius_fraction"],
            row["point_count"],
        ),
    )
    key = (float(selected["theta_max_deg"]), float(selected["eta_radius_fraction"]))
    selected_branches, fidelity_rows = candidates[key]
    derived_points = np.vstack(
        [np.asarray(branch.points_um, dtype=float) for branch in selected_branches]
    )
    port_errors = [
        float(
            np.min(
                np.linalg.norm(
                    derived_points - np.asarray(port.cap_center_um)[None, :], axis=1
                )
            )
        )
        for port in ports
    ]
    port_position_error = max(port_errors, default=0.0)
    topology_preserved = (
        len(selected_branches) == len(branches)
        and [branch.branch_id for branch in selected_branches]
        == [branch.branch_id for branch in branches]
    )
    report = {
        "method": "CFD_DERIVED_SPLINE_CENTERLINE",
        "source_swc_modified": False,
        "raw_points_um_modified": False,
        "raw_radius_um_modified": False,
        "source_topology_modified": not topology_preserved,
        "source_topology_preserved": topology_preserved,
        "legacy_moving_average_used": False,
        "curve": "centripetal cubic Hermite (C1)",
        "radius": "PCHIP",
        "selection_basis": "measured tangent, source fidelity, sagitta, then point count",
        "selected_theta_max_deg": key[0],
        "selected_eta_radius_fraction": key[1],
        "adaptive_sensitivity": sensitivity,
        "branch_fidelity": fidelity_rows,
        "junction_position_error_um": selected["junction_position_error_um"],
        "endpoint_position_error_um": selected["endpoint_position_error_um"],
        "port_position_error_um": port_position_error,
        "port_position_error_per_port_um": port_errors,
        "port_position_exact": port_position_error == 0.0,
        "continuous_c1_tangent_jump_deg": 0.0,
        "status": (
            "PASS"
            if topology_preserved and port_position_error == 0.0
            else "FAIL"
        ),
    }
    return SmoothCenterlineBuild(
        branches=selected_branches,
        port_tail_rows=tail_rows,
        selected_tangent_angle_deg=key[0],
        selected_sagitta_radius_fraction=key[1],
        sensitivity_rows=sensitivity,
        branch_fidelity_rows=fidelity_rows,
        report=report,
    )
