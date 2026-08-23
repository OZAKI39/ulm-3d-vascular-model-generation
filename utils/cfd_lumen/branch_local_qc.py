"""Radius-normalized surface ownership and branch-local cross-section QC."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh
from scipy.spatial import cKDTree
from shapely import Polygon, contains_xy

from .config import CFDLumenConfig
from .export import write_csv, write_json
from .local_implicit_junction import sample_from_junction
from .surface_qc import _orthogonal_basis, _section_polygon
from .types import BranchGeometry, HybridBuildDetails


BRANCH_FIDELITY_REGION = 0
JUNCTION_CORE_REGION = 1


class _BranchDistanceIndex:
    """Spatially shortlist segments, then retain exact segment/radius projection."""

    def __init__(self, branch: BranchGeometry) -> None:
        self.branch = branch
        points = np.asarray(branch.points_um, dtype=float)
        self.radii = np.asarray(branch.radius_um, dtype=float)
        self.starts = points[:-1]
        self.vectors = points[1:] - self.starts
        self.squared_lengths = np.einsum("ij,ij->i", self.vectors, self.vectors)
        self.segment_lengths = np.sqrt(self.squared_lengths)
        self.arc_starts = np.asarray(branch.arc_length_um[:-1], dtype=float)
        self.tree = cKDTree(self.starts + 0.5 * self.vectors)
        self.maximum_half_length = 0.5 * float(self.segment_lengths.max())
        self.maximum_radius = float(self.radii.max())

    def _evaluate_candidates(
        self,
        queries: np.ndarray,
        candidate_indices: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        if candidate_indices.ndim == 1:
            candidate_indices = candidate_indices[:, None]
        starts = self.starts[candidate_indices]
        vectors = self.vectors[candidate_indices]
        relative = queries[:, None, :] - starts
        projection = np.einsum("qki,qki->qk", relative, vectors)
        projection /= self.squared_lengths[candidate_indices]
        projection = np.clip(projection, 0.0, 1.0)
        closest = starts + projection[:, :, None] * vectors
        distance = np.linalg.norm(queries[:, None, :] - closest, axis=2)
        local_radius = self.radii[candidate_indices] + projection * (
            self.radii[candidate_indices + 1] - self.radii[candidate_indices]
        )
        normalized = distance / np.maximum(local_radius, np.finfo(float).eps)
        best_column = np.argmin(normalized, axis=1)
        row = np.arange(len(queries))
        segment = candidate_indices[row, best_column]
        arc = (
            self.arc_starts[segment]
            + projection[row, best_column] * self.segment_lengths[segment]
        )
        return normalized[row, best_column], arc

    def _evaluate_all(
        self,
        queries: np.ndarray,
        *,
        chunk_size: int = 2_048,
    ) -> tuple[np.ndarray, np.ndarray]:
        distance = np.empty(len(queries), dtype=float)
        arc = np.empty(len(queries), dtype=float)
        all_segments = np.arange(len(self.starts), dtype=np.int64)
        for begin in range(0, len(queries), chunk_size):
            end = min(begin + chunk_size, len(queries))
            candidates = np.broadcast_to(
                all_segments[None, :], (end - begin, len(all_segments))
            )
            distance[begin:end], arc[begin:end] = self._evaluate_candidates(
                queries[begin:end], candidates
            )
        return distance, arc

    def query(self, queries: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        segment_count = len(self.starts)
        candidate_count = min(32, segment_count)
        if candidate_count == segment_count:
            return self._evaluate_all(queries)
        midpoint_distance, candidate_indices = self.tree.query(
            queries, k=candidate_count + 1, workers=-1
        )
        best, arc = self._evaluate_candidates(
            queries, np.asarray(candidate_indices)[:, :candidate_count]
        )
        # Every omitted segment is at least next-midpoint distance minus the
        # largest half-segment away. Dividing by the largest source radius gives
        # a conservative lower bound, so unresolved points alone need all
        # segments. The returned result therefore remains exact, not an ANN QC.
        omitted_lower_bound = np.maximum(
            np.asarray(midpoint_distance)[:, candidate_count] - self.maximum_half_length,
            0.0,
        ) / self.maximum_radius
        unresolved = omitted_lower_bound < best
        if np.any(unresolved):
            best[unresolved], arc[unresolved] = self._evaluate_all(queries[unresolved])
        return best, arc


def _all_branch_distances(
    queries: np.ndarray,
    indices: list[_BranchDistanceIndex],
) -> tuple[np.ndarray, np.ndarray]:
    normalized = np.empty((len(queries), len(indices)), dtype=float)
    arcs = np.empty_like(normalized)
    for column, index in enumerate(indices):
        normalized[:, column], arcs[:, column] = index.query(queries)
    return normalized, arcs


def _ownership(
    normalized: np.ndarray,
    branches: list[BranchGeometry],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    owner_column = np.argmin(normalized, axis=1)
    minimum = normalized[np.arange(len(normalized)), owner_column]
    if normalized.shape[1] > 1:
        second = np.partition(normalized, 1, axis=1)[:, 1]
        margin = second - minimum
    else:
        margin = np.full(len(normalized), np.inf, dtype=float)
    branch_ids = np.asarray([branch.branch_id for branch in branches], dtype=np.int32)
    return branch_ids[owner_column], owner_column, margin


def _region_labels(
    owner_column: np.ndarray,
    arcs: np.ndarray,
    branches: list[BranchGeometry],
    details: HybridBuildDetails,
) -> np.ndarray:
    regions = np.full(len(owner_column), BRANCH_FIDELITY_REGION, dtype=np.uint8)
    collars: dict[int, list[Any]] = {}
    for patch in details.patches.values():
        for collar in patch.collars:
            collars.setdefault(collar.branch_id, []).append(collar)
    for column, branch in enumerate(branches):
        selected = owner_column == column
        if not np.any(selected):
            continue
        branch_arc = arcs[selected, column]
        is_core = np.zeros(len(branch_arc), dtype=bool)
        for collar in collars.get(branch.branch_id, []):
            distance = (
                branch_arc
                if collar.endpoint_index == 0
                else float(branch.arc_length_um[-1]) - branch_arc
            )
            is_core |= distance <= collar.implicit_extent_um
        selected_indices = np.flatnonzero(selected)
        regions[selected_indices[is_core]] = JUNCTION_CORE_REGION
    return regions


def _surface_ownership(
    mesh: trimesh.Trimesh,
    branches: list[BranchGeometry],
    indices: list[_BranchDistanceIndex],
    details: HybridBuildDetails,
    margin_threshold: float,
) -> tuple[pv.PolyData, dict[str, Any], np.ndarray]:
    vertices = np.asarray(mesh.vertices, dtype=float)
    centers = np.asarray(mesh.triangles_center, dtype=float)
    vertex_distance, vertex_arcs = _all_branch_distances(vertices, indices)
    face_distance, face_arcs = _all_branch_distances(centers, indices)
    vertex_owner, vertex_column, vertex_margin = _ownership(vertex_distance, branches)
    face_owner, face_column, face_margin = _ownership(face_distance, branches)
    vertex_region = _region_labels(vertex_column, vertex_arcs, branches, details)
    face_region = _region_labels(face_column, face_arcs, branches, details)
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces, dtype=np.int64))
    ).ravel()
    polydata = pv.PolyData(vertices, faces)
    polydata.point_data["owner_branch"] = vertex_owner
    polydata.point_data["ownership_margin"] = vertex_margin
    polydata.point_data["ambiguous_junction_surface"] = (
        vertex_margin < margin_threshold
    ).astype(np.uint8)
    polydata.point_data["qc_region"] = vertex_region
    polydata.cell_data["owner_branch"] = face_owner
    polydata.cell_data["ownership_margin"] = face_margin
    polydata.cell_data["ambiguous_junction_surface"] = (
        face_margin < margin_threshold
    ).astype(np.uint8)
    polydata.cell_data["qc_region"] = face_region
    for column, branch in enumerate(branches):
        name = f"d_norm_branch_{branch.branch_id}"
        polydata.point_data[name] = vertex_distance[:, column]
        polydata.cell_data[name] = face_distance[:, column]
    report = {
        "surface_point_count": int(len(vertices)),
        "surface_triangle_count": int(len(mesh.faces)),
        "point_ambiguous_count": int(np.count_nonzero(vertex_margin < margin_threshold)),
        "triangle_ambiguous_count": int(np.count_nonzero(face_margin < margin_threshold)),
        "point_ambiguous_fraction": float(np.mean(vertex_margin < margin_threshold)),
        "triangle_ambiguous_fraction": float(np.mean(face_margin < margin_threshold)),
        "branch_fidelity_triangle_count": int(
            np.count_nonzero(face_region == BRANCH_FIDELITY_REGION)
        ),
        "junction_core_triangle_count": int(
            np.count_nonzero(face_region == JUNCTION_CORE_REGION)
        ),
    }
    return polydata, report, face_region


def _branch_local_section(
    mesh: trimesh.Trimesh,
    center: np.ndarray,
    tangent: np.ndarray,
    source_radius: float,
    target_branch_id: int,
    branches: list[BranchGeometry],
    indices: list[_BranchDistanceIndex],
    config: CFDLumenConfig,
) -> dict[str, float | int | None]:
    section = _section_polygon(mesh, center, tangent)
    if section is None:
        return {
            "old_total_plane_area_um2": None,
            "new_branch_local_area_um2": None,
            "ownership_confidence": None,
            "classified_area_fraction": None,
            "multi_branch_contamination_area_um2": None,
            "ambiguous_area_um2": None,
            "interior_branch_partition_area_um2": None,
            "surface_owned_boundary_fraction": None,
            "cross_section_grid_point_count": 0,
        }
    old_area, coordinates = section
    polygon = Polygon(coordinates)
    boundary = np.asarray(coordinates, dtype=float)
    boundary_midpoints = 0.5 * (boundary[:-1] + boundary[1:])
    first, second = _orthogonal_basis(tangent)
    boundary_points = (
        center[None, :]
        + boundary_midpoints[:, 0, None] * first[None, :]
        + boundary_midpoints[:, 1, None] * second[None, :]
    )
    boundary_distance, _ = _all_branch_distances(boundary_points, indices)
    boundary_owner, _, boundary_margin = _ownership(boundary_distance, branches)
    edge_lengths = np.linalg.norm(np.diff(boundary, axis=0), axis=1)
    boundary_confident = boundary_margin >= config.branch_local_qc.ownership_margin_threshold
    boundary_target = boundary_owner == target_branch_id
    selected_edges = boundary_target & boundary_confident
    target_length = float(edge_lengths[boundary_target].sum())
    target_confident_length = float(edge_lengths[selected_edges].sum())
    local_area = 0.0
    if np.all(selected_edges):
        local_area = float(old_area)
    elif np.any(selected_edges):
        edge_count = len(selected_edges)
        false_indices = np.flatnonzero(~selected_edges)
        cursor = int((false_indices[0] + 1) % edge_count)
        run: list[np.ndarray] = []
        for offset in range(edge_count):
            edge_index = (cursor + offset) % edge_count
            if selected_edges[edge_index]:
                if not run:
                    run.append(boundary[edge_index])
                run.append(boundary[edge_index + 1])
            elif run:
                if len(run) >= 3:
                    component = Polygon(np.asarray(run, dtype=float))
                    if component.is_valid and not component.is_empty:
                        local_area += float(component.area)
                run = []
        if run and len(run) >= 3:
            component = Polygon(np.asarray(run, dtype=float))
            if component.is_valid and not component.is_empty:
                local_area += float(component.area)
    minimum_x, minimum_y, maximum_x, maximum_y = polygon.bounds
    qc = config.branch_local_qc
    spacing = 2.0 * source_radius / qc.cross_section_cells_across_source_diameter
    columns = max(1, int(np.ceil((maximum_x - minimum_x) / spacing)))
    rows = max(1, int(np.ceil((maximum_y - minimum_y) / spacing)))
    if columns * rows > qc.maximum_cross_section_grid_points:
        scale = np.sqrt(columns * rows / qc.maximum_cross_section_grid_points)
        spacing *= scale
        columns = max(1, int(np.ceil((maximum_x - minimum_x) / spacing)))
        rows = max(1, int(np.ceil((maximum_y - minimum_y) / spacing)))
    x = minimum_x + (np.arange(columns) + 0.5) * (maximum_x - minimum_x) / columns
    y = minimum_y + (np.arange(rows) + 0.5) * (maximum_y - minimum_y) / rows
    grid_x, grid_y = np.meshgrid(x, y, indexing="xy")
    flat_x = grid_x.ravel()
    flat_y = grid_y.ravel()
    inside = contains_xy(polygon, flat_x, flat_y)
    inside_count = int(np.count_nonzero(inside))
    if inside_count == 0:
        return {
            "old_total_plane_area_um2": float(old_area),
            "new_branch_local_area_um2": None,
            "ownership_confidence": None,
            "classified_area_fraction": None,
            "multi_branch_contamination_area_um2": None,
            "ambiguous_area_um2": None,
            "interior_branch_partition_area_um2": None,
            "surface_owned_boundary_fraction": (
                target_confident_length / float(edge_lengths.sum())
                if float(edge_lengths.sum()) > 0
                else None
            ),
            "cross_section_grid_point_count": int(columns * rows),
        }
    plane_points = (
        center[None, :]
        + flat_x[inside, None] * first[None, :]
        + flat_y[inside, None] * second[None, :]
    )
    normalized, _ = _all_branch_distances(plane_points, indices)
    owner, _, margin = _ownership(normalized, branches)
    confident = margin >= qc.ownership_margin_threshold
    target = owner == target_branch_id
    target_confident = int(np.count_nonzero(target & confident))
    other_count = int(np.count_nonzero((~target) & confident))
    ambiguous_count = int(np.count_nonzero(~confident))
    scale_area = float(old_area) / inside_count
    return {
        "old_total_plane_area_um2": float(old_area),
        "new_branch_local_area_um2": local_area,
        "ownership_confidence": (
            target_confident_length / target_length if target_length > 0 else 0.0
        ),
        "classified_area_fraction": float(
            edge_lengths[boundary_confident].sum() / edge_lengths.sum()
        ),
        "multi_branch_contamination_area_um2": scale_area * other_count,
        "ambiguous_area_um2": scale_area * ambiguous_count,
        "interior_branch_partition_area_um2": scale_area * target_confident,
        "surface_owned_boundary_fraction": float(
            target_confident_length / edge_lengths.sum()
        ),
        "cross_section_grid_point_count": int(columns * rows),
    }


def _cross_section_rows(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
    indices: list[_BranchDistanceIndex],
    config: CFDLumenConfig,
    *,
    samples_per_branch: int = 25,
) -> list[dict[str, Any]]:
    branch_by_id = {branch.branch_id: branch for branch in branches}
    rows: list[dict[str, Any]] = []
    for patch in details.patches.values():
        for collar in patch.collars:
            branch = branch_by_id[collar.branch_id]
            diameter = 2.0 * collar.collar_radius_um
            maximum = min(3.0 * diameter, 0.9 * branch.length_um)
            minimum = min(0.1 * diameter, 0.5 * maximum)
            for sample_index, distance in enumerate(
                np.linspace(minimum, maximum, samples_per_branch)
            ):
                center, source_radius, tangent = sample_from_junction(
                    branch, collar.endpoint_index, float(distance)
                )
                section = _branch_local_section(
                    mesh,
                    center,
                    tangent,
                    source_radius,
                    branch.branch_id,
                    branches,
                    indices,
                    config,
                )
                source_area = float(np.pi * source_radius**2)
                old_area = section["old_total_plane_area_um2"]
                local_area = section["new_branch_local_area_um2"]
                region = (
                    "JUNCTION_CORE_REGION"
                    if distance <= collar.implicit_extent_um
                    else "BRANCH_FIDELITY_REGION"
                )
                confidence = section["ownership_confidence"]
                confidence_pass = bool(
                    confidence is not None
                    and confidence >= config.branch_local_qc.minimum_ownership_confidence
                )
                local_ratio = local_area / source_area if local_area is not None else None
                rows.append(
                    {
                        "junction_node_id": collar.junction_node_id,
                        "branch_id": collar.branch_id,
                        "sample_index": sample_index,
                        "distance_from_junction_um": float(distance),
                        "distance_in_diameters": float(distance / diameter),
                        "qc_region": region,
                        "collar_distance_um": collar.collar_distance_um,
                        "implicit_extent_um": collar.implicit_extent_um,
                        "source_radius_um": source_radius,
                        "source_area_um2": source_area,
                        **section,
                        "old_total_plane_area_ratio": (
                            old_area / source_area if old_area is not None else None
                        ),
                        "new_branch_local_area_ratio": local_ratio,
                        "new_branch_local_equivalent_radius_um": (
                            float(np.sqrt(local_area / np.pi))
                            if local_area is not None and local_area >= 0
                            else None
                        ),
                        "branch_local_radius_relative_error": (
                            float(np.sqrt(local_area / np.pi) / source_radius - 1.0)
                            if local_area is not None and local_area >= 0
                            else None
                        ),
                        "ownership_confidence_pass": confidence_pass,
                        "severe_branch_local_abnormality": bool(
                            region == "BRANCH_FIDELITY_REGION"
                            and confidence_pass
                            and local_ratio is not None
                            and abs(local_ratio - 1.0)
                            > config.branch_local_qc.severe_area_relative_error
                        ),
                    }
                )
    return rows


def _junction_summaries(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for junction_id in sorted({int(row["junction_node_id"]) for row in rows}):
        valid = [
            row
            for row in rows
            if row["junction_node_id"] == junction_id
            and row["old_total_plane_area_ratio"] is not None
            and row["new_branch_local_area_ratio"] is not None
        ]
        if not valid:
            continue
        peak = max(valid, key=lambda row: row["old_total_plane_area_ratio"])
        old_ratio = float(peak["old_total_plane_area_ratio"])
        new_ratio = float(peak["new_branch_local_area_ratio"])
        source_area = float(peak["source_area_um2"])
        other_branch_ratio = float(peak["multi_branch_contamination_area_um2"]) / source_area
        ambiguous_ratio = float(peak["ambiguous_area_um2"]) / source_area
        removed = other_branch_ratio + ambiguous_ratio
        summaries.append(
            {
                "junction_node_id": junction_id,
                "old_peak_total_plane_area_ratio": old_ratio,
                "new_branch_local_ratio_at_old_peak": new_ratio,
                "ratio_attributed_to_other_branches": other_branch_ratio,
                "ratio_attributed_to_ambiguous_surface": ambiguous_ratio,
                "ratio_attributed_to_other_branches_or_ambiguity": removed,
                "old_peak_other_branch_fraction": (
                    other_branch_ratio / old_ratio if old_ratio else None
                ),
                "old_peak_ambiguous_fraction": (
                    ambiguous_ratio / old_ratio if old_ratio else None
                ),
                "old_peak_total_contamination_fraction": (
                    removed / old_ratio if old_ratio else None
                ),
                "peak_branch_id": peak["branch_id"],
                "peak_distance_from_junction_um": peak["distance_from_junction_um"],
                "peak_qc_region": peak["qc_region"],
                "peak_ownership_confidence": peak["ownership_confidence"],
            }
        )
    return summaries


def _fidelity_summaries(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for junction_id in sorted({int(row["junction_node_id"]) for row in rows}):
        fidelity = [
            row
            for row in rows
            if row["junction_node_id"] == junction_id
            and row["qc_region"] == "BRANCH_FIDELITY_REGION"
            and row["ownership_confidence_pass"]
            and row["new_branch_local_area_ratio"] is not None
        ]
        summaries.append(
            {
                "junction_node_id": junction_id,
                "evaluated_cross_section_count": len(fidelity),
                "maximum_area_absolute_relative_error": max(
                    (
                        abs(float(row["new_branch_local_area_ratio"]) - 1.0)
                        for row in fidelity
                    ),
                    default=None,
                ),
                "maximum_radius_absolute_relative_error": max(
                    (
                        abs(float(row["branch_local_radius_relative_error"]))
                        for row in fidelity
                        if row["branch_local_radius_relative_error"] is not None
                    ),
                    default=None,
                ),
                "minimum_ownership_confidence": min(
                    (float(row["ownership_confidence"]) for row in fidelity),
                    default=None,
                ),
            }
        )
    return summaries


def evaluate_branch_local_cross_section_qc(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
    config: CFDLumenConfig,
) -> tuple[pv.PolyData, list[dict[str, Any]], dict[str, Any]]:
    """Evaluate v4 branch ownership, region-specific area fidelity, and core safety."""

    indices = [_BranchDistanceIndex(branch) for branch in branches]
    polydata, surface_report, face_region = _surface_ownership(
        mesh,
        branches,
        indices,
        details,
        config.branch_local_qc.ownership_margin_threshold,
    )
    rows = _cross_section_rows(mesh, details, branches, indices, config)
    fidelity = [row for row in rows if row["qc_region"] == "BRANCH_FIDELITY_REGION"]
    confidence_failures = [row for row in fidelity if not row["ownership_confidence_pass"]]
    severe = [row for row in fidelity if row["severe_branch_local_abnormality"]]
    core = [
        row
        for row in rows
        if row["qc_region"] == "JUNCTION_CORE_REGION"
        and row["old_total_plane_area_ratio"] is not None
    ]
    minimum_core_ratio = min(
        (float(row["old_total_plane_area_ratio"]) for row in core), default=None
    )
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    core_adjacency = (
        (face_region[adjacency[:, 0]] == JUNCTION_CORE_REGION)
        & (face_region[adjacency[:, 1]] == JUNCTION_CORE_REGION)
        if len(adjacency)
        else np.empty(0, dtype=bool)
    )
    core_angles = angles[core_adjacency]
    normal_jump_p99 = float(np.percentile(core_angles, 99)) if len(core_angles) else None
    checks = {
        "surface_ownership_complete": bool(
            surface_report["surface_point_count"] == len(mesh.vertices)
            and surface_report["surface_triangle_count"] == len(mesh.faces)
        ),
        "branch_fidelity_ownership_confidence": len(confidence_failures) == 0,
        "no_severe_branch_local_area_abnormality": len(severe) == 0,
        "no_extreme_artificial_throat": bool(
            minimum_core_ratio is None
            or minimum_core_ratio
            >= config.branch_local_qc.minimum_junction_total_area_ratio
        ),
        "junction_surface_smoothness": bool(
            normal_jump_p99 is None
            or normal_jump_p99
            <= config.branch_local_qc.maximum_junction_normal_jump_p99_deg
        ),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "method": "BRANCH_LOCAL_CROSS_SECTION_QC",
        "distance_definition": "d_norm_i = distance_to_branch_i / interpolated_source_radius_i",
        "branch_local_area_definition": (
            "area enclosed by confident owner_branch surface-section arcs, with straight "
            "closure at ownership transitions"
        ),
        "contamination_area_definition": (
            "radius-normalized branch partition of the selected total-plane lumen section"
        ),
        "ambiguous_label": "AMBIGUOUS_JUNCTION_SURFACE",
        "region_codes": {
            str(BRANCH_FIDELITY_REGION): "BRANCH_FIDELITY_REGION",
            str(JUNCTION_CORE_REGION): "JUNCTION_CORE_REGION",
        },
        "ownership_margin_threshold": config.branch_local_qc.ownership_margin_threshold,
        "minimum_ownership_confidence": (
            config.branch_local_qc.minimum_ownership_confidence
        ),
        "severe_area_relative_error_threshold": (
            config.branch_local_qc.severe_area_relative_error
        ),
        "checks": checks,
        **surface_report,
        "cross_section_count": len(rows),
        "branch_fidelity_cross_section_count": len(fidelity),
        "junction_core_cross_section_count": len(core),
        "ownership_confidence_failure_count": len(confidence_failures),
        "severe_branch_local_abnormality_count": len(severe),
        "controlled_local_implicit_required": len(severe) > 0,
        "branch_fidelity_max_area_absolute_relative_error": max(
            (
                abs(float(row["new_branch_local_area_ratio"]) - 1.0)
                for row in fidelity
                if row["new_branch_local_area_ratio"] is not None
                and row["ownership_confidence_pass"]
            ),
            default=None,
        ),
        "branch_fidelity_max_radius_absolute_relative_error": max(
            (
                abs(float(row["branch_local_radius_relative_error"]))
                for row in fidelity
                if row["branch_local_radius_relative_error"] is not None
                and row["ownership_confidence_pass"]
            ),
            default=None,
        ),
        "junction_core_minimum_total_plane_area_ratio": minimum_core_ratio,
        "junction_core_normal_jump_p99_deg": normal_jump_p99,
        "junction_area_ratio_is_descriptive_not_hard_upper_bound": True,
        "junction_summaries": _junction_summaries(rows),
        "branch_fidelity_summaries": _fidelity_summaries(rows),
    }
    return polydata, rows, report


def write_branch_local_qc_artifacts(
    root: Path,
    polydata: pv.PolyData,
    rows: list[dict[str, Any]],
    report: dict[str, Any],
) -> list[Path]:
    root.mkdir(parents=True, exist_ok=True)
    surface_path = root / "surface_branch_ownership.vtp"
    rows_path = root / "branch_local_cross_sections.csv"
    report_path = root / "branch_local_cross_section_qc.json"
    polydata.save(surface_path)
    write_csv(rows_path, rows)
    write_json(report_path, report)
    return [surface_path, rows_path, report_path]
