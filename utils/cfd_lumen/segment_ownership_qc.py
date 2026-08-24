"""Segment-level PolyBall ownership and tangent/defect correlation for v9."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh
from scipy.spatial import cKDTree

from .config import CFDLumenConfig
from .types import BranchGeometry
from .unified_polyball import JunctionBlendSpec, PolyBallLineModel
from .v8_qc import _junction_edge_mask


CROSS_BRANCH_SWITCH = "CROSS_BRANCH_SWITCH"
SAME_BRANCH_ADJACENT_SEGMENT_SWITCH = "SAME_BRANCH_ADJACENT_SEGMENT_SWITCH"
SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH = "SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH"


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


def _correlation(first: np.ndarray, second: np.ndarray) -> float | None:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    finite = np.isfinite(first) & np.isfinite(second)
    if np.count_nonzero(finite) < 2:
        return None
    first = first[finite]
    second = second[finite]
    if np.std(first) <= 0.0 or np.std(second) <= 0.0:
        return None
    return float(np.corrcoef(first, second)[0, 1])


def _polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def _near_fraction(
    query_points: np.ndarray, reference_points: np.ndarray, distance: float
) -> float:
    if not len(query_points) or not len(reference_points):
        return 0.0
    nearest = cKDTree(reference_points).query(query_points, k=1)[0]
    return float(np.mean(nearest <= distance))


def _switch_type(
    first_branch: int,
    second_branch: int,
    first_index: int,
    second_index: int,
) -> str:
    if first_branch != second_branch:
        return CROSS_BRANCH_SWITCH
    if abs(first_index - second_index) == 1:
        return SAME_BRANCH_ADJACENT_SEGMENT_SWITCH
    return SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH


def _write_ownership(
    mesh: trimesh.Trimesh,
    ownership: Any,
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = _polydata(mesh)
    arrays = {
        "winner_segment_id": ownership.winner_segment_id,
        "winner_branch_id": ownership.winner_branch_id,
        "winner_segment_index_in_branch": ownership.winner_segment_index_in_branch,
        "winner_segment_parametric_t": ownership.winner_parametric_t,
        "winner_local_radius_um": ownership.winner_local_radius_um,
        "second_segment_id": ownership.second_segment_id,
        "second_branch_id": ownership.second_branch_id,
        "second_segment_index_in_branch": ownership.second_segment_index_in_branch,
        "second_segment_parametric_t": ownership.second_parametric_t,
        "ownership_margin_um": ownership.ownership_margin_um,
    }
    for name, values in arrays.items():
        data.point_data[name] = np.asarray(values)
    data.save(path)
    return path


def write_switch_edges(
    mesh: trimesh.Trimesh,
    switch_edges: np.ndarray,
    rows: list[dict[str, Any]],
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    edges = np.asarray(switch_edges, dtype=np.int64).reshape((-1, 2))
    if not len(edges):
        pv.PolyData().save(path)
        return path
    points = np.asarray(mesh.vertices, dtype=float)[edges].reshape((-1, 3))
    lines = np.column_stack(
        (
            np.full(len(edges), 2, dtype=np.int64),
            2 * np.arange(len(edges), dtype=np.int64),
            2 * np.arange(len(edges), dtype=np.int64) + 1,
        )
    ).ravel()
    data = pv.PolyData(points, lines=lines)
    type_code = {
        CROSS_BRANCH_SWITCH: 1,
        SAME_BRANCH_ADJACENT_SEGMENT_SWITCH: 2,
        SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH: 3,
    }
    data.cell_data["switch_type"] = np.asarray(
        [type_code[row["switch_type"]] for row in rows], dtype=np.int8
    )
    for key in (
        "winner_segment_0",
        "winner_segment_1",
        "winner_branch_0",
        "winner_branch_1",
        "gradient_jump_angle_deg",
        "surface_dihedral_angle_deg",
        "is_sawtooth_defect",
    ):
        data.cell_data[key] = np.asarray([row[key] for row in rows])
    data.save(path)
    return path


def evaluate_segment_ownership(
    mesh: trimesh.Trimesh,
    model: PolyBallLineModel,
    specs: tuple[JunctionBlendSpec, ...],
    config: CFDLumenConfig,
    *,
    ownership_path: Path | None = None,
    field_model: Any | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, np.ndarray]]:
    """Classify A/B/C switching on a raw or final unified surface."""

    vertices = np.asarray(mesh.vertices, dtype=float)
    ownership = model.evaluate_with_ownership(
        vertices, k=config.v7.k_nearest_segments, gradients=False
    )
    if ownership_path is not None:
        _write_ownership(mesh, ownership, ownership_path)
    adjacency_edges = np.asarray(mesh.face_adjacency_edges, dtype=np.int64)
    dihedral = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    junction = _junction_edge_mask(mesh, specs)
    local_dihedral = dihedral[junction]
    threshold = (
        float(np.percentile(local_dihedral, config.v8.sawtooth_dihedral_percentile))
        if len(local_dihedral)
        else np.inf
    )
    defect = junction & (dihedral >= threshold)
    switch = (
        ownership.winner_segment_id[adjacency_edges[:, 0]]
        != ownership.winner_segment_id[adjacency_edges[:, 1]]
    )
    switch_ids = np.flatnonzero(switch & junction)
    defect_ids = np.flatnonzero(defect)
    edge_points = vertices[adjacency_edges]
    midpoint = edge_points.mean(axis=1)
    median_edge = float(
        np.median(np.linalg.norm(edge_points[:, 1] - edge_points[:, 0], axis=1))
    )
    distance_limit = config.v8.switch_overlap_distance_edge_factors * median_edge

    first_vertex = adjacency_edges[switch_ids, 0]
    second_vertex = adjacency_edges[switch_ids, 1]
    first_segment = ownership.winner_segment_id[first_vertex]
    second_segment = ownership.winner_segment_id[second_vertex]
    first_branch = ownership.winner_branch_id[first_vertex]
    second_branch = ownership.winner_branch_id[second_vertex]
    first_index = ownership.winner_segment_index_in_branch[first_vertex]
    second_index = ownership.winner_segment_index_in_branch[second_vertex]
    switch_types = np.asarray(
        [
            _switch_type(int(a), int(b), int(i), int(j))
            for a, b, i, j in zip(
                first_branch, second_branch, first_index, second_index
            )
        ],
        dtype=object,
    )
    selected_midpoint = midpoint[switch_ids]
    _, first_gradient, first_t, _ = model.evaluate_segments(
        selected_midpoint, first_segment, gradients=True
    )
    _, second_gradient, second_t, _ = model.evaluate_segments(
        selected_midpoint, second_segment, gradients=True
    )
    assert first_gradient is not None and second_gradient is not None
    denominator = np.linalg.norm(first_gradient, axis=1) * np.linalg.norm(
        second_gradient, axis=1
    )
    cosine = np.einsum("ij,ij->i", first_gradient, second_gradient) / np.maximum(
        denominator, 1.0e-15
    )
    gradient_jump = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    # For a competition-aware smooth union the hard model still defines exact
    # ownership, while cross-branch continuity must be measured on the actual
    # field used by FlyingEdges and Newton projection.  Same-branch rows retain
    # the selected-segment analytical gradients so residual segment creases are
    # not hidden by sampling two nearby surface vertices.
    actual_field = field_model or model
    cross = switch_types == CROSS_BRANCH_SWITCH
    if actual_field is not model and np.any(cross):
        _, gradient_0 = actual_field.evaluate(
            vertices[first_vertex[cross]],
            k=config.v7.k_nearest_segments,
            gradients=True,
        )
        _, gradient_1 = actual_field.evaluate(
            vertices[second_vertex[cross]],
            k=config.v7.k_nearest_segments,
            gradients=True,
        )
        assert gradient_0 is not None and gradient_1 is not None
        field_denominator = np.linalg.norm(gradient_0, axis=1) * np.linalg.norm(
            gradient_1, axis=1
        )
        field_cosine = np.einsum("ij,ij->i", gradient_0, gradient_1) / np.maximum(
            field_denominator, 1.0e-15
        )
        gradient_jump[cross] = np.degrees(
            np.arccos(np.clip(field_cosine, -1.0, 1.0))
        )
    switch_dihedral = dihedral[switch_ids]
    rows: list[dict[str, Any]] = []
    for row_id, edge_id in enumerate(switch_ids):
        v0, v1 = adjacency_edges[edge_id]
        rows.append(
            {
                "edge_id": int(edge_id),
                "vertex_0": int(v0),
                "vertex_1": int(v1),
                "switch_type": str(switch_types[row_id]),
                "winner_segment_0": int(first_segment[row_id]),
                "winner_segment_1": int(second_segment[row_id]),
                "winner_branch_0": int(first_branch[row_id]),
                "winner_branch_1": int(second_branch[row_id]),
                "segment_index_in_branch_0": int(first_index[row_id]),
                "segment_index_in_branch_1": int(second_index[row_id]),
                "segment_parametric_t_0": float(first_t[row_id]),
                "segment_parametric_t_1": float(second_t[row_id]),
                "ownership_margin_0_um": float(ownership.ownership_margin_um[v0]),
                "ownership_margin_1_um": float(ownership.ownership_margin_um[v1]),
                "gradient_jump_angle_deg": float(gradient_jump[row_id]),
                "gradient_measurement": (
                    "actual_field_endpoint_gradient"
                    if actual_field is not model
                    and switch_types[row_id] == CROSS_BRANCH_SWITCH
                    else "analytical_selected_segment_midpoint_gradient"
                ),
                "surface_dihedral_angle_deg": float(switch_dihedral[row_id]),
                "is_sawtooth_defect": bool(defect[edge_id]),
                "midpoint_x_um": float(midpoint[edge_id, 0]),
                "midpoint_y_um": float(midpoint[edge_id, 1]),
                "midpoint_z_um": float(midpoint[edge_id, 2]),
            }
        )

    type_reports: dict[str, Any] = {}
    defect_points = midpoint[defect_ids]
    for switch_type in (
        CROSS_BRANCH_SWITCH,
        SAME_BRANCH_ADJACENT_SEGMENT_SWITCH,
        SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH,
    ):
        mask = switch_types == switch_type
        points = selected_midpoint[mask]
        type_reports[switch_type] = {
            "switch_edge_count": int(np.count_nonzero(mask)),
            "gradient_jump_angle_deg": _summary(gradient_jump[mask]),
            "surface_dihedral_angle_deg": _summary(switch_dihedral[mask]),
            "defect_near_switch_fraction": _near_fraction(
                defect_points, points, distance_limit
            ),
            "switch_near_defect_fraction": _near_fraction(
                points, defect_points, distance_limit
            ),
            "switch_edge_directly_classified_as_defect_count": int(
                np.count_nonzero(defect[switch_ids[mask]])
            ),
        }
    cross_fraction = type_reports[CROSS_BRANCH_SWITCH]["defect_near_switch_fraction"]
    adjacent_fraction = type_reports[SAME_BRANCH_ADJACENT_SEGMENT_SWITCH][
        "defect_near_switch_fraction"
    ]
    nonadjacent_fraction = type_reports[SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH][
        "defect_near_switch_fraction"
    ]
    same_points = selected_midpoint[
        np.isin(
            switch_types,
            (
                SAME_BRANCH_ADJACENT_SEGMENT_SWITCH,
                SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH,
            ),
        )
    ]
    same_fraction = _near_fraction(defect_points, same_points, distance_limit)
    visible = 0
    visible_row_mask = np.zeros(len(rows), dtype=bool)
    for row_id, row in enumerate(rows):
        if not row["is_sawtooth_defect"]:
            continue
        limit = (
            config.v9.maximum_cross_branch_gradient_p99_deg
            if row["switch_type"] == CROSS_BRANCH_SWITCH
            else config.v9.maximum_same_branch_gradient_p99_deg
        )
        is_visible = row["gradient_jump_angle_deg"] > limit
        visible += int(is_visible)
        visible_row_mask[row_id] = is_visible
    confirmed = bool(same_fraction >= 0.5 and same_fraction > cross_fraction)
    report = {
        "surface_vertex_count": int(len(vertices)),
        "surface_adjacency_edge_count": int(len(adjacency_edges)),
        "junction_surface_adjacency_edge_count": int(np.count_nonzero(junction)),
        "segment_switch_edge_count": int(len(switch_ids)),
        "sawtooth_defect_edge_count": int(len(defect_ids)),
        "sawtooth_dihedral_threshold_deg": threshold,
        "switch_overlap_distance_um": distance_limit,
        "types": type_reports,
        "defect_near_cross_branch_fraction": cross_fraction,
        "defect_near_same_branch_adjacent_segment_fraction": adjacent_fraction,
        "defect_near_same_branch_nonadjacent_segment_fraction": nonadjacent_fraction,
        "defect_near_any_same_branch_segment_fraction": same_fraction,
        "visible_sawtooth_count": int(visible),
        "gradient_jump_surface_dihedral_correlation": _correlation(
            gradient_jump, switch_dihedral
        ),
        "root_cause_confirmed": confirmed,
        "root_cause_if_confirmed": (
            "PIECEWISE_LINEAR_CENTERLINE_SEGMENT_CREASE" if confirmed else None
        ),
        "cross_branch_gradient_field": (
            "actual_extraction_and_projection_field"
            if actual_field is not model
            else "hard_polyball_field"
        ),
    }
    arrays = {
        "winner_segment_id": ownership.winner_segment_id,
        "winner_branch_id": ownership.winner_branch_id,
        "second_segment_id": ownership.second_segment_id,
        "second_branch_id": ownership.second_branch_id,
        "ownership_margin_um": ownership.ownership_margin_um,
        "switch_edges": adjacency_edges[switch_ids],
        "switch_types": switch_types,
        "defect_edges": adjacency_edges[defect_ids],
        "visible_defect_edges": adjacency_edges[switch_ids[visible_row_mask]],
        "switch_edge_midpoints": selected_midpoint,
        "defect_edge_midpoints": defect_points,
    }
    return report, rows, arrays


def tangent_audit_rows(
    model: PolyBallLineModel,
    branches: list[BranchGeometry],
    specs: tuple[JunctionBlendSpec, ...],
    switch_rows: list[dict[str, Any]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Audit adjacent segment kinks and correlate them with surface defects."""

    branch_by_id = {branch.branch_id: branch for branch in branches}
    rows: list[dict[str, Any]] = []
    kink_by_pair: dict[tuple[int, int], float] = {}
    raw_kink_by_pair: dict[tuple[int, int], float | None] = {}
    for branch_id in model.branch_ids:
        segment_ids = np.flatnonzero(model.segment_branch_id == branch_id)
        order = np.argsort(model.segment_index_in_branch[segment_ids])
        segment_ids = segment_ids[order]
        for first_id, second_id in zip(segment_ids[:-1], segment_ids[1:]):
            first_index = int(model.segment_index_in_branch[first_id])
            second_index = int(model.segment_index_in_branch[second_id])
            if second_index - first_index != 1:
                continue
            first_vector = model.segment_end_local[first_id] - model.segment_start_local[first_id]
            second_vector = model.segment_end_local[second_id] - model.segment_start_local[second_id]
            angle = float(
                np.degrees(
                    np.arccos(
                        np.clip(
                            np.dot(first_vector, second_vector)
                            / max(
                                np.linalg.norm(first_vector) * np.linalg.norm(second_vector),
                                1.0e-15,
                            ),
                            -1.0,
                            1.0,
                        )
                    )
                )
            )
            world = model.local_to_world(model.segment_end_local[first_id][None, :])[0]
            junction_distance = min(
                (
                    float(np.linalg.norm(world - np.asarray(spec.center_world_um)))
                    for spec in specs
                ),
                default=np.inf,
            )
            branch = branch_by_id.get(int(branch_id))
            neighborhood = ""
            raw_kink: float | None = None
            if branch is not None and len(branch.raw_points_um):
                nearest = int(
                    np.argmin(np.linalg.norm(branch.raw_points_um - world[None, :], axis=1))
                )
                lower = max(0, nearest - 1)
                upper = min(len(branch.local_node_ids), nearest + 2)
                neighborhood = ";".join(
                    map(str, branch.local_node_ids[lower:upper])
                )
                if len(branch.raw_points_um) >= 3:
                    raw_index = min(max(nearest, 1), len(branch.raw_points_um) - 2)
                    raw_first = (
                        branch.raw_points_um[raw_index]
                        - branch.raw_points_um[raw_index - 1]
                    )
                    raw_second = (
                        branch.raw_points_um[raw_index + 1]
                        - branch.raw_points_um[raw_index]
                    )
                    raw_kink = float(
                        np.degrees(
                            np.arccos(
                                np.clip(
                                    np.dot(raw_first, raw_second)
                                    / max(
                                        np.linalg.norm(raw_first)
                                        * np.linalg.norm(raw_second),
                                        1.0e-15,
                                    ),
                                    -1.0,
                                    1.0,
                                )
                            )
                        )
                    )
            radius = 0.5 * (
                float(model.radius_end_um[first_id])
                + float(model.radius_start_um[second_id])
            )
            row = {
                "branch_id": int(branch_id),
                "source_node_neighborhood": neighborhood,
                "segment_id_0": int(first_id),
                "segment_id_1": int(second_id),
                "segment_index_in_branch_0": first_index,
                "segment_index_in_branch_1": second_index,
                "theta_deg": angle,
                "source_raw_theta_deg": raw_kink,
                "local_radius_um": radius,
                "distance_to_junction_um": junction_distance,
                "x_um": float(world[0]),
                "y_um": float(world[1]),
                "z_um": float(world[2]),
            }
            rows.append(row)
            pair = tuple(sorted((int(first_id), int(second_id))))
            kink_by_pair[pair] = angle
            raw_kink_by_pair[pair] = raw_kink
    correlation_rows = switch_rows or []
    kink = np.asarray(
        [
            kink_by_pair.get(
                tuple(sorted((int(row["winner_segment_0"]), int(row["winner_segment_1"])))),
                np.nan,
            )
            for row in correlation_rows
        ],
        dtype=float,
    )
    raw_kink = np.asarray(
        [
            raw_kink_by_pair.get(
                tuple(sorted((int(row["winner_segment_0"]), int(row["winner_segment_1"])))),
                np.nan,
            )
            for row in correlation_rows
        ],
        dtype=float,
    )
    for row, angle, raw_angle in zip(correlation_rows, kink, raw_kink):
        row["source_tangent_kink_deg"] = None if not np.isfinite(angle) else float(angle)
        row["source_raw_tangent_kink_deg"] = (
            None if not np.isfinite(raw_angle) else float(raw_angle)
        )
    gradient = np.asarray(
        [float(row["gradient_jump_angle_deg"]) for row in correlation_rows], dtype=float
    )
    dihedral = np.asarray(
        [float(row["surface_dihedral_angle_deg"]) for row in correlation_rows], dtype=float
    )
    defect = np.asarray(
        [float(bool(row["is_sawtooth_defect"])) for row in correlation_rows], dtype=float
    )
    maximum = max(rows, key=lambda row: row["theta_deg"], default=None)
    report = {
        "adjacent_segment_pair_count": len(rows),
        "theta_deg": _summary(np.asarray([row["theta_deg"] for row in rows])),
        "source_raw_theta_deg": _summary(
            np.asarray(
                [
                    row["source_raw_theta_deg"]
                    for row in rows
                    if row["source_raw_theta_deg"] is not None
                ]
            )
        ),
        "maximum_kink": maximum,
        "correlation": {
            "theta_vs_segment_switch_gradient_jump": _correlation(kink, gradient),
            "theta_vs_surface_dihedral": _correlation(kink, dihedral),
            "theta_vs_sawtooth_defect": _correlation(kink, defect),
            "source_raw_theta_vs_segment_switch_gradient_jump": _correlation(
                raw_kink, gradient
            ),
            "source_raw_theta_vs_surface_dihedral": _correlation(raw_kink, dihedral),
            "source_raw_theta_vs_sawtooth_defect": _correlation(raw_kink, defect),
            "matched_switch_count": int(np.count_nonzero(np.isfinite(kink))),
        },
    }
    return rows, report


def competition_support_report(
    mesh: trimesh.Trimesh,
    field_model: Any,
    specs: tuple[JunctionBlendSpec, ...],
) -> dict[str, Any]:
    vertices = np.asarray(mesh.vertices, dtype=float)
    rows: list[dict[str, Any]] = []
    for spec in specs:
        distance = np.linalg.norm(
            vertices - np.asarray(spec.center_world_um)[None, :], axis=1
        )
        selected = np.flatnonzero(distance < spec.blend_length_um)
        active, margin, radius = field_model.competition_mask(vertices[selected], spec)
        rows.append(
            {
                "junction_node_id": spec.junction_node_id,
                "junction_vertex_count": int(len(selected)),
                "competition_active_vertex_count": int(np.count_nonzero(active)),
                "competition_active_fraction": float(np.mean(active)) if len(active) else 0.0,
                "competition_margin_um": _summary(margin),
                "competition_margin_radius_fraction": _summary(margin / radius),
            }
        )
    return {
        "support": "true junction neighborhood AND phi_2-phi_1 < threshold*local_radius",
        "per_junction": rows,
        "total_active_vertex_count": int(
            sum(row["competition_active_vertex_count"] for row in rows)
        ),
    }
