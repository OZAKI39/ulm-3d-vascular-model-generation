"""Field-ownership, stage-localization, and morphology QC for the v8 protocol."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh
from scipy.spatial import cKDTree

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .types import HybridBuildDetails, RadiusFidelitySample
from .unified_polyball import JunctionBlendSpec, PolyBallLineModel


def build_junction_blend_specs(
    roi: ROIRecord,
    details: HybridBuildDetails,
    config: CFDLumenConfig,
) -> tuple[JunctionBlendSpec, ...]:
    """Use the already validated v6 true-junction cores; ports are not involved."""

    specs: list[JunctionBlendSpec] = []
    for node_id, patch in sorted(details.patches.items()):
        incident = tuple(sorted({int(collar.branch_id) for collar in patch.collars}))
        if len(incident) < 2:
            continue
        blend_length = config.v8.blend_length_scale * max(
            float(collar.implicit_extent_um) for collar in patch.collars
        )
        specs.append(
            JunctionBlendSpec(
                junction_node_id=int(node_id),
                center_world_um=np.asarray(
                    roi.local_node_positions_um[node_id], dtype=float
                ),
                radius_um=float(roi.local_node_radius_um[node_id]),
                blend_length_um=float(blend_length),
                incident_branch_ids=incident,
            )
        )
    if not specs:
        raise ValueError("V8 requires at least one true junction with two incident branches")
    return tuple(specs)


def _polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def _summary(values: np.ndarray) -> dict[str, float | int | None]:
    values = np.asarray(values, dtype=float)
    if not len(values):
        return {"count": 0, "mean": None, "p95": None, "p99": None, "max": None}
    return {
        "count": int(len(values)),
        "mean": float(np.mean(values)),
        "p95": float(np.percentile(values, 95)),
        "p99": float(np.percentile(values, 99)),
        "max": float(np.max(values)),
    }


def _junction_edge_mask(
    mesh: trimesh.Trimesh, specs: tuple[JunctionBlendSpec, ...]
) -> np.ndarray:
    edge_points = np.asarray(mesh.vertices, dtype=float)[mesh.face_adjacency_edges]
    midpoint = edge_points.mean(axis=1)
    return np.any(
        np.column_stack(
            [
                np.linalg.norm(
                    midpoint - np.asarray(spec.center_world_um)[None, :], axis=1
                )
                <= spec.blend_length_um
                for spec in specs
            ]
        ),
        axis=1,
    )


def evaluate_ownership_switching(
    mesh: trimesh.Trimesh,
    model: PolyBallLineModel,
    specs: tuple[JunctionBlendSpec, ...],
    config: CFDLumenConfig,
    ownership_path: Path | None = None,
    *,
    include_gradient_details: bool = True,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, np.ndarray]]:
    """Compute branch ownership, switch edges, gradient jumps, and defect overlap."""

    vertices = np.asarray(mesh.vertices, dtype=float)
    branch_phi, _, branch_ids = model.evaluate_branch_fields(
        vertices, gradients=False
    )
    order = np.argsort(branch_phi, axis=1, kind="stable")
    rows = np.arange(len(vertices))
    winner_column = order[:, 0]
    second_column = order[:, 1]
    winner = branch_ids[winner_column]
    second = branch_ids[second_column]
    margin = branch_phi[rows, second_column] - branch_phi[rows, winner_column]

    faces = np.asarray(mesh.faces, dtype=np.int64)
    face_phi = branch_phi[faces].mean(axis=1)
    face_order = np.argsort(face_phi, axis=1, kind="stable")
    face_rows = np.arange(len(faces))
    face_winner = branch_ids[face_order[:, 0]]
    face_second = branch_ids[face_order[:, 1]]
    face_margin = (
        face_phi[face_rows, face_order[:, 1]]
        - face_phi[face_rows, face_order[:, 0]]
    )
    if ownership_path is not None:
        ownership_path.parent.mkdir(parents=True, exist_ok=True)
        polydata = _polydata(mesh)
        polydata.point_data["winner_branch"] = winner.astype(np.int64)
        polydata.point_data["second_best_branch"] = second.astype(np.int64)
        polydata.point_data["ownership_margin"] = margin.astype(float)
        polydata.cell_data["winner_branch"] = face_winner.astype(np.int64)
        polydata.cell_data["second_best_branch"] = face_second.astype(np.int64)
        polydata.cell_data["ownership_margin"] = face_margin.astype(float)
        polydata.save(ownership_path)

    adjacency_edges = np.asarray(mesh.face_adjacency_edges, dtype=np.int64)
    adjacency_angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    switch = winner[adjacency_edges[:, 0]] != winner[adjacency_edges[:, 1]]
    junction = _junction_edge_mask(mesh, specs)
    local_angles = adjacency_angles[junction]
    threshold = (
        float(np.percentile(local_angles, config.v8.sawtooth_dihedral_percentile))
        if len(local_angles)
        else np.inf
    )
    defect = junction & (adjacency_angles >= threshold)
    local_switch = switch & junction
    edge_points = vertices[adjacency_edges]
    midpoint = edge_points.mean(axis=1)
    median_edge = float(np.median(np.linalg.norm(np.diff(edge_points, axis=1)[:, 0], axis=1)))
    distance_limit = config.v8.switch_overlap_distance_edge_factors * median_edge
    defect_ids = np.flatnonzero(defect)
    switch_ids = np.flatnonzero(local_switch)
    if len(defect_ids) and len(switch_ids):
        switch_tree = cKDTree(midpoint[switch_ids])
        defect_distance = switch_tree.query(midpoint[defect_ids], k=1)[0]
        defect_near_switch = defect_distance <= distance_limit
        defect_coverage = float(np.mean(defect_near_switch))
        defect_tree = cKDTree(midpoint[defect_ids])
        switch_distance = defect_tree.query(midpoint[switch_ids], k=1)[0]
        switch_near_defect = switch_distance <= distance_limit
        switch_coverage = float(np.mean(switch_near_defect))
        switch_near_mask = np.zeros(len(adjacency_edges), dtype=bool)
        switch_near_mask[switch_ids[switch_near_defect]] = True
    else:
        defect_coverage = 0.0
        switch_coverage = 0.0
        switch_near_mask = np.zeros(len(adjacency_edges), dtype=bool)

    switch_rows: list[dict[str, Any]] = []
    gradient_array = np.empty(0, dtype=float)
    dihedral_array = np.asarray(adjacency_angles[switch_ids], dtype=float)
    if include_gradient_details and len(switch_ids):
        switch_midpoints = midpoint[switch_ids]
        _, midpoint_gradients, evaluated_ids = model.evaluate_branch_fields(
            switch_midpoints, gradients=True
        )
        assert midpoint_gradients is not None
        column_by_branch = {
            int(branch_id): column for column, branch_id in enumerate(evaluated_ids)
        }
        first_vertices = adjacency_edges[switch_ids, 0]
        second_vertices = adjacency_edges[switch_ids, 1]
        first_branches = winner[first_vertices]
        second_branches = winner[second_vertices]
        first_columns = np.asarray(
            [column_by_branch[int(branch)] for branch in first_branches]
        )
        second_columns = np.asarray(
            [column_by_branch[int(branch)] for branch in second_branches]
        )
        switch_row_ids = np.arange(len(switch_ids))
        first_gradient = midpoint_gradients[switch_row_ids, first_columns]
        second_gradient = midpoint_gradients[switch_row_ids, second_columns]
        denominator = np.linalg.norm(first_gradient, axis=1) * np.linalg.norm(
            second_gradient, axis=1
        )
        cosine = np.einsum("ij,ij->i", first_gradient, second_gradient) / np.maximum(
            denominator, 1.0e-15
        )
        gradient_array = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
        for row_index, edge_id in enumerate(switch_ids):
            first_vertex, second_vertex = adjacency_edges[edge_id]
            switch_rows.append(
                {
                    "edge_id": int(edge_id),
                    "vertex_0": int(first_vertex),
                    "vertex_1": int(second_vertex),
                    "winner_branch_0": int(first_branches[row_index]),
                    "winner_branch_1": int(second_branches[row_index]),
                    "ownership_margin_0_um": float(margin[first_vertex]),
                    "ownership_margin_1_um": float(margin[second_vertex]),
                    "gradient_jump_angle_deg": float(gradient_array[row_index]),
                    "surface_dihedral_angle_deg": float(dihedral_array[row_index]),
                    "is_sawtooth_defect": bool(defect[edge_id]),
                    "switch_within_defect_distance": bool(switch_near_mask[edge_id]),
                    "midpoint_x_um": float(midpoint[edge_id, 0]),
                    "midpoint_y_um": float(midpoint[edge_id, 1]),
                    "midpoint_z_um": float(midpoint[edge_id, 2]),
                }
            )
    correlation = (
        float(np.corrcoef(gradient_array, dihedral_array)[0, 1])
        if len(gradient_array) >= 2
        and np.std(gradient_array) > 0.0
        and np.std(dihedral_array) > 0.0
        else None
    )
    confirmed = bool(
        len(switch_ids) and max(defect_coverage, switch_coverage) >= 0.5
    )
    report = {
        "surface_vertex_count": int(len(vertices)),
        "branch_count": int(len(branch_ids)),
        "branch_ids": branch_ids.astype(int).tolist(),
        "ownership_switch_edge_count": int(len(switch_ids)),
        "junction_surface_adjacency_edge_count": int(np.count_nonzero(junction)),
        "sawtooth_dihedral_percentile": config.v8.sawtooth_dihedral_percentile,
        "sawtooth_dihedral_threshold_deg": threshold,
        "sawtooth_defect_edge_count": int(len(defect_ids)),
        "defect_near_switch_fraction": defect_coverage,
        "switch_near_defect_fraction": switch_coverage,
        "switch_overlap_distance_um": distance_limit,
        "gradient_jump_angle_deg": _summary(gradient_array),
        "surface_dihedral_at_switch_deg": _summary(dihedral_array),
        "gradient_jump_surface_dihedral_correlation": correlation,
        "hard_min_switch_crease_confirmed": confirmed,
        "root_cause_if_confirmed": (
            "HARD_MIN_POLYBALL_SWITCH_CREASE" if confirmed else None
        ),
    }
    arrays = {
        "winner": winner,
        "second": second,
        "margin": margin,
        "switch_edges": adjacency_edges[switch_ids],
        "switch_edge_midpoints": midpoint[switch_ids],
        "defect_edges": adjacency_edges[defect_ids],
        "defect_edge_midpoints": midpoint[defect_ids],
    }
    return report, switch_rows, arrays


def local_normal_metrics(
    mesh: trimesh.Trimesh, specs: tuple[JunctionBlendSpec, ...]
) -> dict[str, Any]:
    mask = _junction_edge_mask(mesh, specs)
    angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))[mask]
    return {"normal_jump_deg": _summary(angles)}


def field_gradient_angles_on_edges(
    mesh: trimesh.Trimesh,
    field_model: Any,
    edges: np.ndarray,
) -> dict[str, float | int | None]:
    edges = np.asarray(edges, dtype=np.int64).reshape((-1, 2))
    if not len(edges):
        return _summary(np.empty(0, dtype=float))
    points = np.asarray(mesh.vertices, dtype=float)[edges].reshape((-1, 3))
    _, gradients = field_model.evaluate(points, gradients=True)
    assert gradients is not None
    paired = gradients.reshape((-1, 2, 3))
    denominator = np.linalg.norm(paired[:, 0], axis=1) * np.linalg.norm(
        paired[:, 1], axis=1
    )
    cosine = np.einsum("ij,ij->i", paired[:, 0], paired[:, 1]) / np.maximum(
        denominator, 1.0e-15
    )
    angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return _summary(angles)


def stage_localization_rows(
    stage_meshes: dict[str, trimesh.Trimesh],
    model: PolyBallLineModel,
    specs: tuple[JunctionBlendSpec, ...],
    config: CFDLumenConfig,
) -> tuple[list[dict[str, Any]], str | None]:
    rows: list[dict[str, Any]] = []
    order = (
        "S0_raw_flying_edges",
        "S1_newton_projected",
        "S2_pyacvd_before_second_projection",
        "S3_final_projected_before_port_clip",
    )
    for stage in order:
        mesh = stage_meshes[stage]
        ownership, _, _ = evaluate_ownership_switching(
            mesh,
            model,
            specs,
            config,
            ownership_path=None,
            include_gradient_details=False,
        )
        edge_mask = _junction_edge_mask(mesh, specs)
        edge_angles = np.degrees(
            np.asarray(mesh.face_adjacency_angles, dtype=float)
        )[edge_mask]
        large_dihedral_count = int(
            np.count_nonzero(edge_angles >= config.v6.silhouette_large_corner_deg)
        )
        rows.append(
            {
                "stage": stage,
                **local_normal_metrics(mesh, specs),
                "ownership_switch_edge_count": ownership[
                    "ownership_switch_edge_count"
                ],
                "sawtooth_defect_edge_count": ownership[
                    "sawtooth_defect_edge_count"
                ],
                "defect_near_switch_fraction": ownership[
                    "defect_near_switch_fraction"
                ],
                "switch_near_defect_fraction": ownership[
                    "switch_near_defect_fraction"
                ],
                "hard_min_switch_crease_present": ownership[
                    "hard_min_switch_crease_confirmed"
                ],
                "large_dihedral_edge_count": large_dihedral_count,
                "sawtooth_present": large_dihedral_count > 0,
            }
        )
    earliest = next(
        (row["stage"] for row in rows if row["sawtooth_present"]),
        None,
    )
    return rows, earliest


def junction_local_volumes(
    model: Any,
    specs: tuple[JunctionBlendSpec, ...],
    config: CFDLumenConfig,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    samples = int(config.v8.junction_volume_samples_across_diameter)
    for spec in specs:
        spacing = 2.0 * spec.blend_length_um / samples
        coordinate = (
            -spec.blend_length_um
            + (np.arange(samples, dtype=float) + 0.5) * spacing
        )
        zz, yy, xx = np.meshgrid(coordinate, coordinate, coordinate, indexing="ij")
        relative = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        core = np.linalg.norm(relative, axis=1) <= spec.blend_length_um
        points = relative[core] + np.asarray(spec.center_world_um)[None, :]
        phi, _ = model.evaluate(points, gradients=False)
        volume = float(np.count_nonzero(phi <= 0.0) * spacing**3)
        rows.append(
            {
                "junction_node_id": spec.junction_node_id,
                "blend_length_um": spec.blend_length_um,
                "sample_spacing_um": spacing,
                "sample_count_in_blend_sphere": int(np.count_nonzero(core)),
                "junction_local_volume_um3": volume,
            }
        )
    return rows


def hydraulic_resistance_comparison(
    baseline: list[RadiusFidelitySample],
    candidate: list[RadiusFidelitySample],
) -> dict[str, Any]:
    """Compare Poiseuille geometry factors; viscosity and 8/pi cancel."""

    baseline_by_key = {
        (sample.branch_id, sample.sample_index): sample for sample in baseline
    }
    candidate_by_key = {
        (sample.branch_id, sample.sample_index): sample for sample in candidate
    }
    rows: list[dict[str, Any]] = []
    for branch_id in sorted({key[0] for key in baseline_by_key}):
        keys = sorted(
            (
                key
                for key in baseline_by_key
                if key[0] == branch_id and key in candidate_by_key
            ),
            key=lambda key: baseline_by_key[key].arc_length_um,
        )
        if len(keys) < 2:
            continue
        arc = np.asarray([baseline_by_key[key].arc_length_um for key in keys])
        hard_radius = np.asarray(
            [baseline_by_key[key].reconstructed_radius_um for key in keys]
        )
        smooth_radius = np.asarray(
            [candidate_by_key[key].reconstructed_radius_um for key in keys]
        )
        hard_factor = float(np.trapz(1.0 / hard_radius**4, arc))
        smooth_factor = float(np.trapz(1.0 / smooth_radius**4, arc))
        relative_error = (smooth_factor - hard_factor) / hard_factor
        rows.append(
            {
                "branch_id": int(branch_id),
                "matched_sample_count": len(keys),
                "v7_resistance_geometry_factor_um-3": hard_factor,
                "candidate_resistance_geometry_factor_um-3": smooth_factor,
                "relative_error": float(relative_error),
                "absolute_relative_error": float(abs(relative_error)),
            }
        )
    absolute = np.asarray([row["absolute_relative_error"] for row in rows])
    return {
        "equation": "R=(8*mu/pi)*integral(ds/r^4); reported comparison cancels 8*mu/pi",
        "per_branch": rows,
        "p95_absolute_relative_error": (
            float(np.percentile(absolute, 95)) if len(absolute) else None
        ),
        "max_absolute_relative_error": float(np.max(absolute)) if len(absolute) else None,
    }


def local_hydraulic_resistance_from_area_profiles(
    baseline: list[dict[str, Any]],
    candidate: list[dict[str, Any]],
) -> dict[str, Any]:
    """Integrate the junction-local Poiseuille factor from measured areas."""

    baseline_by_key = {
        (
            int(row["junction_node_id"]),
            int(row["branch_id"]),
            int(row["sample_index"]),
        ): row
        for row in baseline
    }
    candidate_by_key = {
        (
            int(row["junction_node_id"]),
            int(row["branch_id"]),
            int(row["sample_index"]),
        ): row
        for row in candidate
    }
    branch_keys = sorted({key[:2] for key in baseline_by_key})
    rows: list[dict[str, Any]] = []
    for junction_id, branch_id in branch_keys:
        keys = sorted(
            (
                key
                for key in baseline_by_key
                if key[:2] == (junction_id, branch_id)
                and key in candidate_by_key
                and baseline_by_key[key]["cross_section_area_um2"] is not None
                and candidate_by_key[key]["cross_section_area_um2"] is not None
            ),
            key=lambda key: float(baseline_by_key[key]["distance_from_junction_um"]),
        )
        if len(keys) < 2:
            continue
        distance = np.asarray(
            [baseline_by_key[key]["distance_from_junction_um"] for key in keys],
            dtype=float,
        )
        hard_area = np.asarray(
            [baseline_by_key[key]["cross_section_area_um2"] for key in keys],
            dtype=float,
        )
        smooth_area = np.asarray(
            [candidate_by_key[key]["cross_section_area_um2"] for key in keys],
            dtype=float,
        )
        # 1/r^4 = pi^2/A^2; the common constants cancel in the comparison.
        hard_factor = float(np.trapz(np.pi**2 / hard_area**2, distance))
        smooth_factor = float(np.trapz(np.pi**2 / smooth_area**2, distance))
        relative_error = (smooth_factor - hard_factor) / hard_factor
        rows.append(
            {
                "junction_node_id": junction_id,
                "branch_id": branch_id,
                "matched_sample_count": len(keys),
                "v7_resistance_geometry_factor_um-3": hard_factor,
                "candidate_resistance_geometry_factor_um-3": smooth_factor,
                "relative_error": float(relative_error),
                "absolute_relative_error": float(abs(relative_error)),
            }
        )
    absolute = np.asarray([row["absolute_relative_error"] for row in rows])
    return {
        "equation": "R=(8*mu/pi)*integral(ds/r^4), r=sqrt(A/pi); common constants cancel",
        "scope": "true-junction branch-local area profiles",
        "per_junction_branch": rows,
        "p95_absolute_relative_error": (
            float(np.percentile(absolute, 95)) if len(absolute) else None
        ),
        "max_absolute_relative_error": float(np.max(absolute)) if len(absolute) else None,
    }


def compare_volume_rows(
    baseline: list[dict[str, Any]], candidate: list[dict[str, Any]]
) -> dict[str, Any]:
    baseline_by_id = {int(row["junction_node_id"]): row for row in baseline}
    rows: list[dict[str, Any]] = []
    for row in candidate:
        node_id = int(row["junction_node_id"])
        hard = float(baseline_by_id[node_id]["junction_local_volume_um3"])
        smooth = float(row["junction_local_volume_um3"])
        rows.append(
            {
                "junction_node_id": node_id,
                "v7_volume_um3": hard,
                "candidate_volume_um3": smooth,
                "relative_change": float((smooth - hard) / hard),
            }
        )
    return {
        "per_junction": rows,
        "maximum_relative_increase": max(
            (float(row["relative_change"]) for row in rows), default=None
        ),
        "mean_relative_change": (
            float(np.mean([row["relative_change"] for row in rows])) if rows else None
        ),
    }
