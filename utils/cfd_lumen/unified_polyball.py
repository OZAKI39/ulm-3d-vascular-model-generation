"""One-field, one-iso-surface vascular reconstruction for the v7 protocol."""

from __future__ import annotations

import importlib.util
import gc
import os
import shutil
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh
import vtk
from scipy.spatial import cKDTree
from scipy.optimize import linear_sum_assignment
from skimage.measure import marching_cubes
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from .config import CFDLumenConfig
from .mesh_defects import triangle_quality
from .types import (
    BranchGeometry,
    GeometryValidationError,
    HybridBuildDetails,
    PatchResult,
    PortGeometry,
)


@dataclass(slots=True)
class PolyBallOwnership:
    """Exact first/second segment ownership without changing ``evaluate``."""

    phi_um: np.ndarray
    gradient_world: np.ndarray | None
    winner_segment_id: np.ndarray
    winner_branch_id: np.ndarray
    winner_segment_index_in_branch: np.ndarray
    winner_parametric_t: np.ndarray
    winner_local_radius_um: np.ndarray
    second_segment_id: np.ndarray
    second_branch_id: np.ndarray
    second_segment_index_in_branch: np.ndarray
    second_parametric_t: np.ndarray
    ownership_margin_um: np.ndarray


@dataclass(slots=True)
class PolyBallLineModel:
    """Piecewise-linear centerline with continuously interpolated radius."""

    segment_start_local: np.ndarray
    segment_end_local: np.ndarray
    radius_start_um: np.ndarray
    radius_end_um: np.ndarray
    segment_branch_id: np.ndarray
    segment_index_in_branch: np.ndarray
    origin_world_um: np.ndarray
    local_axes: np.ndarray
    tree: cKDTree
    provider: str

    @property
    def segment_count(self) -> int:
        return int(len(self.segment_start_local))

    def world_to_local(self, points: np.ndarray) -> np.ndarray:
        return (np.asarray(points, dtype=float) - self.origin_world_um) @ self.local_axes.T

    def local_to_world(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(points, dtype=float) @ self.local_axes + self.origin_world_um

    def evaluate(
        self,
        points_world_um: np.ndarray,
        *,
        k: int = 32,
        gradients: bool = False,
        chunk_size: int = 50_000,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """Evaluate exact variable-radius segment envelopes near arbitrary points."""

        query_world = np.asarray(points_world_um, dtype=float).reshape((-1, 3))
        values = np.empty(len(query_world), dtype=float)
        gradient_world = np.empty_like(query_world) if gradients else None
        candidate_count = min(max(1, int(k)), self.segment_count)
        for start in range(0, len(query_world), chunk_size):
            stop = min(start + chunk_size, len(query_world))
            query = self.world_to_local(query_world[start:stop])
            _, candidates = self.tree.query(query, k=candidate_count, workers=-1)
            if candidate_count == 1:
                candidates = np.asarray(candidates, dtype=np.int64)[:, None]
            phi, grad_local = _candidate_polyball_values(
                query,
                np.asarray(candidates, dtype=np.int64),
                self.segment_start_local,
                self.segment_end_local,
                self.radius_start_um,
                self.radius_end_um,
                gradients=gradients,
            )
            values[start:stop] = phi
            if gradients and gradient_world is not None and grad_local is not None:
                gradient_world[start:stop] = grad_local @ self.local_axes
        return values, gradient_world

    def evaluate_with_ownership(
        self,
        points_world_um: np.ndarray,
        *,
        k: int = 32,
        gradients: bool = False,
        chunk_size: int = 50_000,
    ) -> PolyBallOwnership:
        """Return exact segment winners, runners-up, margins, and segment ``t``."""

        query_world = np.asarray(points_world_um, dtype=float).reshape((-1, 3))
        count = len(query_world)
        values = np.empty(count, dtype=float)
        gradient_world = np.empty_like(query_world) if gradients else None
        winner_segment = np.empty(count, dtype=np.int64)
        winner_t = np.empty(count, dtype=float)
        winner_radius = np.empty(count, dtype=float)
        second_segment = np.full(count, -1, dtype=np.int64)
        second_t = np.full(count, np.nan, dtype=float)
        margin = np.full(count, np.inf, dtype=float)
        candidate_count = min(max(1, int(k)), self.segment_count)
        for start in range(0, count, chunk_size):
            stop = min(start + chunk_size, count)
            query = self.world_to_local(query_world[start:stop])
            _, candidates = self.tree.query(query, k=candidate_count, workers=-1)
            if candidate_count == 1:
                candidates = np.asarray(candidates, dtype=np.int64)[:, None]
            else:
                candidates = np.asarray(candidates, dtype=np.int64)
            candidate_phi, parametric_t, radial, distance = (
                _candidate_polyball_components(
                    query,
                    candidates,
                    self.segment_start_local,
                    self.segment_end_local,
                    self.radius_start_um,
                    self.radius_end_um,
                )
            )
            order = np.argsort(candidate_phi, axis=1, kind="stable")
            rows = np.arange(len(query))
            winner_column = order[:, 0]
            winner_ids = candidates[rows, winner_column]
            values[start:stop] = candidate_phi[rows, winner_column]
            winner_segment[start:stop] = winner_ids
            winner_t[start:stop] = parametric_t[rows, winner_column]
            winner_radius[start:stop] = (
                self.radius_start_um[winner_ids]
                + winner_t[start:stop]
                * (self.radius_end_um[winner_ids] - self.radius_start_um[winner_ids])
            )
            if gradients and gradient_world is not None:
                chosen_radial = radial[rows, winner_column]
                chosen_distance = distance[rows, winner_column]
                gradient_local = np.divide(
                    chosen_radial,
                    chosen_distance[:, None],
                    out=np.zeros_like(chosen_radial),
                    where=chosen_distance[:, None] > 1.0e-14,
                )
                gradient_world[start:stop] = gradient_local @ self.local_axes
            if candidate_count > 1:
                second_column = order[:, 1]
                second_ids = candidates[rows, second_column]
                second_segment[start:stop] = second_ids
                second_t[start:stop] = parametric_t[rows, second_column]
                margin[start:stop] = (
                    candidate_phi[rows, second_column]
                    - candidate_phi[rows, winner_column]
                )
        valid_second = second_segment >= 0
        second_branch = np.full(count, -1, dtype=np.int64)
        second_index = np.full(count, -1, dtype=np.int64)
        second_branch[valid_second] = self.segment_branch_id[
            second_segment[valid_second]
        ]
        second_index[valid_second] = self.segment_index_in_branch[
            second_segment[valid_second]
        ]
        return PolyBallOwnership(
            phi_um=values,
            gradient_world=gradient_world,
            winner_segment_id=winner_segment,
            winner_branch_id=self.segment_branch_id[winner_segment],
            winner_segment_index_in_branch=self.segment_index_in_branch[
                winner_segment
            ],
            winner_parametric_t=winner_t,
            winner_local_radius_um=winner_radius,
            second_segment_id=second_segment,
            second_branch_id=second_branch,
            second_segment_index_in_branch=second_index,
            second_parametric_t=second_t,
            ownership_margin_um=margin,
        )

    def evaluate_segments(
        self,
        points_world_um: np.ndarray,
        segment_ids: np.ndarray,
        *,
        gradients: bool = False,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray]:
        """Evaluate one explicitly selected source segment per query point."""

        points = np.asarray(points_world_um, dtype=float).reshape((-1, 3))
        ids = np.asarray(segment_ids, dtype=np.int64).reshape(-1)
        if len(points) != len(ids) or np.any((ids < 0) | (ids >= self.segment_count)):
            raise ValueError("segment_ids must contain one valid segment per point")
        query = self.world_to_local(points)
        candidates = ids[:, None]
        candidate_phi, parametric_t, radial, distance = (
            _candidate_polyball_components(
                query,
                candidates,
                self.segment_start_local,
                self.segment_end_local,
                self.radius_start_um,
                self.radius_end_um,
            )
        )
        phi = candidate_phi[:, 0]
        t = parametric_t[:, 0]
        local_radius = self.radius_start_um[ids] + t * (
            self.radius_end_um[ids] - self.radius_start_um[ids]
        )
        gradient_world = None
        if gradients:
            gradient_local = np.divide(
                radial[:, 0],
                distance[:, 0, None],
                out=np.zeros_like(radial[:, 0]),
                where=distance[:, 0, None] > 1.0e-14,
            )
            gradient_world = gradient_local @ self.local_axes
        return phi, gradient_world, t, local_radius

    @property
    def branch_ids(self) -> np.ndarray:
        return np.unique(self.segment_branch_id)

    def evaluate_branch_fields(
        self,
        points_world_um: np.ndarray,
        *,
        branch_ids: np.ndarray | None = None,
        gradients: bool = False,
        chunk_size: int = 20_000,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
        """Evaluate the exact minimum field separately for every SWC branch."""

        query_world = np.asarray(points_world_um, dtype=float).reshape((-1, 3))
        selected_branches = (
            self.branch_ids
            if branch_ids is None
            else np.asarray(branch_ids, dtype=np.int64).reshape(-1)
        )
        values = np.full((len(query_world), len(selected_branches)), np.inf, dtype=float)
        gradient_world = (
            np.zeros((len(query_world), len(selected_branches), 3), dtype=float)
            if gradients
            else None
        )
        for column, branch_id in enumerate(selected_branches):
            segment_ids = np.flatnonzero(self.segment_branch_id == branch_id)
            if not len(segment_ids):
                continue
            for start in range(0, len(query_world), chunk_size):
                stop = min(start + chunk_size, len(query_world))
                query = self.world_to_local(query_world[start:stop])
                candidates = np.broadcast_to(
                    segment_ids[None, :], (len(query), len(segment_ids))
                )
                phi, gradient_local = _candidate_polyball_values(
                    query,
                    candidates,
                    self.segment_start_local,
                    self.segment_end_local,
                    self.radius_start_um,
                    self.radius_end_um,
                    gradients=gradients,
                )
                values[start:stop, column] = phi
                if (
                    gradients
                    and gradient_world is not None
                    and gradient_local is not None
                ):
                    gradient_world[start:stop, column] = (
                        gradient_local @ self.local_axes
                    )
        return values, gradient_world, selected_branches


@dataclass(frozen=True, slots=True)
class JunctionBlendSpec:
    junction_node_id: int
    center_world_um: np.ndarray
    radius_um: float
    blend_length_um: float
    incident_branch_ids: tuple[int, ...]

    def compact_weight_and_gradient(
        self, points_world_um: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Quintic compact weight: one at the core and zero with C2 boundary."""

        points = np.asarray(points_world_um, dtype=float).reshape((-1, 3))
        relative = points - np.asarray(self.center_world_um, dtype=float)[None, :]
        distance = np.linalg.norm(relative, axis=1)
        coordinate = np.clip(distance / self.blend_length_um, 0.0, 1.0)
        smootherstep = coordinate**3 * (
            coordinate * (coordinate * 6.0 - 15.0) + 10.0
        )
        weight = 1.0 - smootherstep
        derivative = -30.0 * coordinate**2 * (coordinate - 1.0) ** 2
        gradient = np.zeros_like(relative)
        active = (distance > 1.0e-14) & (distance < self.blend_length_um)
        gradient[active] = (
            derivative[active] / self.blend_length_um / distance[active]
        )[:, None] * relative[active]
        return weight, gradient


@dataclass(slots=True)
class SmoothJunctionPolyBallModel:
    """Hard PolyBall everywhere except compact true-junction smooth unions."""

    hard_model: PolyBallLineModel
    junctions: tuple[JunctionBlendSpec, ...]
    k_radius_ratio: float
    competition_threshold_radius_fraction: float | None = None

    @property
    def radius_start_um(self) -> np.ndarray:
        return self.hard_model.radius_start_um

    @property
    def radius_end_um(self) -> np.ndarray:
        return self.hard_model.radius_end_um

    def evaluate(
        self,
        points_world_um: np.ndarray,
        *,
        k: int = 32,
        gradients: bool = False,
        chunk_size: int = 50_000,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        points = np.asarray(points_world_um, dtype=float).reshape((-1, 3))
        values, output_gradient = self.hard_model.evaluate(
            points, k=k, gradients=gradients, chunk_size=chunk_size
        )
        if not self.junctions or not len(points):
            return values, output_gradient

        normalized = np.column_stack(
            [
                np.linalg.norm(
                    points - np.asarray(spec.center_world_um, dtype=float)[None, :],
                    axis=1,
                )
                / spec.blend_length_um
                for spec in self.junctions
            ]
        )
        nearest = np.argmin(normalized, axis=1)
        inside = np.min(normalized, axis=1) < 1.0
        all_branch_ids = self.hard_model.branch_ids
        for junction_index, spec in enumerate(self.junctions):
            selected = np.flatnonzero(inside & (nearest == junction_index))
            if not len(selected):
                continue
            local_points = points[selected]
            branch_phi, branch_gradient, branch_ids = (
                self.hard_model.evaluate_branch_fields(
                    local_points,
                    branch_ids=all_branch_ids,
                    gradients=gradients,
                    chunk_size=min(chunk_size, 20_000),
                )
            )
            incident_columns = np.flatnonzero(
                np.isin(branch_ids, np.asarray(spec.incident_branch_ids, dtype=np.int64))
            )
            if len(incident_columns) < 2:
                continue
            weight, weight_gradient = spec.compact_weight_and_gradient(local_points)
            incident_values = branch_phi[:, incident_columns]
            incident_order = np.argsort(incident_values, axis=1, kind="stable")
            incident_rows = np.arange(len(local_points))
            competition_margin = (
                incident_values[incident_rows, incident_order[:, 1]]
                - incident_values[incident_rows, incident_order[:, 0]]
            )
            if self.competition_threshold_radius_fraction is None:
                competition_active = weight > 0.0
            else:
                ownership = self.hard_model.evaluate_with_ownership(
                    local_points,
                    k=k,
                    gradients=False,
                    chunk_size=min(chunk_size, 20_000),
                )
                competition_active = (
                    (weight > 0.0)
                    & (
                        competition_margin
                        < self.competition_threshold_radius_fraction
                        * ownership.winner_local_radius_um
                    )
                )
            spatial_k = self.k_radius_ratio * spec.radius_um * weight
            spatial_k_gradient = (
                self.k_radius_ratio * spec.radius_um * weight_gradient
            )
            incident_gradient = (
                branch_gradient[:, incident_columns]
                if branch_gradient is not None
                else None
            )
            smooth_value, smooth_gradient = _stable_smooth_min_reduce(
                branch_phi[:, incident_columns],
                incident_gradient,
                spatial_k,
                spatial_k_gradient,
            )
            other_columns = np.flatnonzero(
                ~np.isin(branch_ids, np.asarray(spec.incident_branch_ids, dtype=np.int64))
            )
            if len(other_columns):
                other_local = np.argmin(branch_phi[:, other_columns], axis=1)
                rows = np.arange(len(selected))
                other_value = branch_phi[rows, other_columns[other_local]]
                use_smooth = smooth_value <= other_value
                combined = np.where(use_smooth, smooth_value, other_value)
                if gradients and output_gradient is not None:
                    assert branch_gradient is not None and smooth_gradient is not None
                    other_gradient = branch_gradient[
                        rows, other_columns[other_local]
                    ]
                    combined_gradient = np.where(
                        use_smooth[:, None], smooth_gradient, other_gradient
                    )
            else:
                combined = smooth_value
                combined_gradient = smooth_gradient
            selected_values = values[selected].copy()
            selected_values[competition_active] = combined[competition_active]
            values[selected] = selected_values
            if gradients and output_gradient is not None:
                assert combined_gradient is not None
                selected_gradient = output_gradient[selected].copy()
                selected_gradient[competition_active] = combined_gradient[
                    competition_active
                ]
                output_gradient[selected] = selected_gradient
        return values, output_gradient

    def competition_mask(
        self, points_world_um: np.ndarray, spec: JunctionBlendSpec
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return v9 competition support, margin, and local radius."""

        points = np.asarray(points_world_um, dtype=float).reshape((-1, 3))
        branch_phi, _, branch_ids = self.hard_model.evaluate_branch_fields(
            points, gradients=False
        )
        incident_columns = np.flatnonzero(
            np.isin(
                branch_ids,
                np.asarray(spec.incident_branch_ids, dtype=np.int64),
            )
        )
        if len(incident_columns) < 2:
            return (
                np.zeros(len(points), dtype=bool),
                np.full(len(points), np.inf),
                np.full(len(points), np.nan),
            )
        values = branch_phi[:, incident_columns]
        order = np.argsort(values, axis=1, kind="stable")
        rows = np.arange(len(points))
        margin = values[rows, order[:, 1]] - values[rows, order[:, 0]]
        ownership = self.hard_model.evaluate_with_ownership(
            points, gradients=False
        )
        weight, _ = spec.compact_weight_and_gradient(points)
        if self.competition_threshold_radius_fraction is None:
            active = weight > 0.0
        else:
            active = (weight > 0.0) & (
                margin
                < self.competition_threshold_radius_fraction
                * ownership.winner_local_radius_um
            )
        return active, margin, ownership.winner_local_radius_um


@dataclass(slots=True)
class UnifiedPolyBallBuild:
    mesh: trimesh.Trimesh
    wall_mesh_before_clip: trimesh.Trimesh
    patch: PatchResult
    face_boundary_id: np.ndarray
    model: PolyBallLineModel
    field_model: PolyBallLineModel | SmoothJunctionPolyBallModel
    constructed_branches: list[BranchGeometry]
    stage_meshes: dict[str, trimesh.Trimesh]
    metadata: dict[str, Any]


@dataclass(slots=True)
class PreparedPolyBallRaster:
    model: PolyBallLineModel
    constructed_branches: list[BranchGeometry]
    port_tail_rows: list[dict[str, Any]]
    minimum_local_um: np.ndarray
    dimensions_xyz: tuple[int, int, int]
    spacing_um: float
    field: np.ndarray
    backing_path: Path | None
    field_report: dict[str, Any]


def _stable_smooth_min_reduce(
    values: np.ndarray,
    gradients: np.ndarray | None,
    spatial_k: np.ndarray,
    spatial_k_gradient: np.ndarray,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Deterministic pairwise reduction using the requested polynomial smin."""

    current = np.asarray(values[:, 0], dtype=float).copy()
    current_gradient = (
        np.asarray(gradients[:, 0], dtype=float).copy()
        if gradients is not None
        else None
    )
    positive = spatial_k > 1.0e-14
    for column in range(1, values.shape[1]):
        following = np.asarray(values[:, column], dtype=float)
        h = np.zeros(len(current), dtype=float)
        h[positive] = np.clip(
            0.5
            + 0.5
            * (following[positive] - current[positive])
            / spatial_k[positive],
            0.0,
            1.0,
        )
        hard_first = current <= following
        h[~positive] = hard_first[~positive].astype(float)
        current = (
            (1.0 - h) * following
            + h * current
            - spatial_k * h * (1.0 - h)
        )
        if current_gradient is not None and gradients is not None:
            following_gradient = np.asarray(gradients[:, column], dtype=float)
            current_gradient = (
                (1.0 - h)[:, None] * following_gradient
                + h[:, None] * current_gradient
                - (h * (1.0 - h))[:, None] * spatial_k_gradient
            )
    return current, current_gradient


def _candidate_polyball_components(
    query: np.ndarray,
    candidates: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    radius_start: np.ndarray,
    radius_end: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return exact per-candidate PolyBall values and minimizing parameters."""

    p0 = starts[candidates]
    direction = ends[candidates] - p0
    delta_radius = radius_end[candidates] - radius_start[candidates]
    relative = query[:, None, :] - p0
    squared_length = np.einsum("nki,nki->nk", direction, direction)
    projection = np.einsum("nki,nki->nk", relative, direction)
    t = np.clip(
        np.divide(
            projection,
            squared_length,
            out=np.zeros_like(projection),
            where=squared_length > 1.0e-20,
        ),
        0.0,
        1.0,
    )
    # The objective is convex for a physically valid tapered segment
    # (|dr| < segment length). Newton converges from the axial projection;
    # clamping also handles endpoint minima exactly.
    for _ in range(6):
        radial = relative - t[:, :, None] * direction
        distance = np.linalg.norm(radial, axis=2)
        numerator = t * squared_length - projection
        derivative = np.divide(
            numerator,
            distance,
            out=np.zeros_like(numerator),
            where=distance > 1.0e-14,
        ) - delta_radius
        second = np.divide(
            squared_length * distance**2 - numerator**2,
            distance**3,
            out=np.zeros_like(distance),
            where=distance > 1.0e-14,
        )
        update = np.divide(
            derivative,
            second,
            out=np.zeros_like(derivative),
            where=second > 1.0e-12,
        )
        t = np.clip(t - update, 0.0, 1.0)
    closest = p0 + t[:, :, None] * direction
    radial = query[:, None, :] - closest
    distance = np.linalg.norm(radial, axis=2)
    local_radius = radius_start[candidates] + t * delta_radius
    candidate_phi = distance - local_radius
    return candidate_phi, t, radial, distance


def _candidate_polyball_values(
    query: np.ndarray,
    candidates: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    radius_start: np.ndarray,
    radius_end: np.ndarray,
    *,
    gradients: bool,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Minimum of exact ``distance-to-line(t) - linear-radius(t)`` values."""

    candidate_phi, _, radial, distance = _candidate_polyball_components(
        query,
        candidates,
        starts,
        ends,
        radius_start,
        radius_end,
    )
    owner = np.argmin(candidate_phi, axis=1)
    rows = np.arange(len(query))
    phi = candidate_phi[rows, owner]
    if not gradients:
        return phi, None
    chosen_radial = radial[rows, owner]
    chosen_distance = distance[rows, owner]
    gradient = np.divide(
        chosen_radial,
        chosen_distance[:, None],
        out=np.zeros_like(chosen_radial),
        where=chosen_distance[:, None] > 1.0e-14,
    )
    return phi, gradient


def _single_segment_phi(
    query: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    radius_start: float,
    radius_end: float,
) -> np.ndarray:
    candidates = np.zeros((len(query), 1), dtype=np.int64)
    values, _ = _candidate_polyball_values(
        query,
        candidates,
        np.asarray(start, dtype=float)[None, :],
        np.asarray(end, dtype=float)[None, :],
        np.asarray((radius_start,), dtype=float),
        np.asarray((radius_end,), dtype=float),
        gradients=False,
    )
    return values


def _copy_with_tail(
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> tuple[list[BranchGeometry], list[dict[str, Any]]]:
    """Append an implicit tail beyond each final CFD plane without changing source data."""

    from dataclasses import replace

    output = [replace(branch) for branch in branches]
    for copied, source in zip(output, branches):
        copied.points_um = np.asarray(source.points_um, dtype=float).copy()
        copied.radius_um = np.asarray(source.radius_um, dtype=float).copy()
        copied.arc_length_um = np.asarray(source.arc_length_um, dtype=float).copy()
    endpoint_rows: list[tuple[float, int, int, PortGeometry]] = []
    for port in ports:
        for branch_index, branch in enumerate(output):
            for endpoint_index in (0, -1):
                endpoint_rows.append(
                    (
                        float(
                            np.linalg.norm(
                                branch.points_um[endpoint_index] - port.cap_center_um
                            )
                        ),
                        branch_index,
                        endpoint_index,
                        port,
                    )
                )
    assignments: list[dict[str, Any]] = []
    used: set[tuple[int, int]] = set()
    for port in ports:
        candidates = [row for row in endpoint_rows if row[3].port_id == port.port_id]
        distance, branch_index, endpoint_index, _ = min(candidates, key=lambda row: row[0])
        tolerance = max(0.10 * port.radius_um, config.v6.port_extension_min_spacing_um)
        if distance > tolerance:
            original_candidates: list[tuple[float, int, int]] = []
            for candidate_index, candidate_branch in enumerate(output):
                for candidate_endpoint in (0, -1):
                    original_candidates.append(
                        (
                            float(
                                np.linalg.norm(
                                    candidate_branch.points_um[candidate_endpoint]
                                    - port.original_position_um
                                )
                            ),
                            candidate_index,
                            candidate_endpoint,
                        )
                    )
            original_distance, branch_index, endpoint_index = min(
                original_candidates, key=lambda row: row[0]
            )
            if original_distance > tolerance:
                raise GeometryValidationError(
                    f"V7_PORT_TAIL_MATCH_FAILED: port {port.port_id}, cap distance "
                    f"{distance:.6g} um, cut distance {original_distance:.6g} um"
                )
            branch = output[branch_index]
            if endpoint_index == -1:
                branch.points_um = np.vstack((branch.points_um, port.cap_center_um))
                branch.radius_um = np.concatenate((branch.radius_um, (port.radius_um,)))
            else:
                branch.points_um = np.vstack((port.cap_center_um, branch.points_um))
                branch.radius_um = np.concatenate(((port.radius_um,), branch.radius_um))
            distance = 0.0
        key = (branch_index, endpoint_index)
        if key in used:
            raise GeometryValidationError(
                f"V7_PORT_TAIL_MATCH_FAILED: endpoint reused for port {port.port_id}"
            )
        used.add(key)
        branch = output[branch_index]
        tail_length = config.v7.port_tail_diameters * 2.0 * port.radius_um
        cap = np.asarray(port.cap_center_um, dtype=float)
        tail = cap + tail_length * np.asarray(port.outward_tangent, dtype=float)
        if endpoint_index == -1:
            branch.points_um = np.vstack((branch.points_um, tail))
            branch.radius_um = np.concatenate((branch.radius_um, (port.radius_um,)))
        else:
            branch.points_um = np.vstack((tail, branch.points_um))
            branch.radius_um = np.concatenate(((port.radius_um,), branch.radius_um))
        branch.arc_length_um = np.concatenate(
            (
                (0.0,),
                np.cumsum(np.linalg.norm(np.diff(branch.points_um, axis=0), axis=1)),
            )
        )
        assignments.append(
            {
                "port_id": port.port_id,
                "cut_port_id": port.cut_port_id,
                "branch_id": branch.branch_id,
                "matched_endpoint_index": 0 if endpoint_index == 0 else -1,
                "endpoint_to_cap_distance_um": distance,
                "port_tail_length_um": tail_length,
                "cap_center_um": cap.tolist(),
                "tail_end_um": tail.tolist(),
                "tail_passes_final_plane": True,
            }
        )
    return output, assignments


def build_polyball_model(
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    v6_details: HybridBuildDetails | None = None,
    constructed_override: list[BranchGeometry] | None = None,
    port_tail_rows_override: list[dict[str, Any]] | None = None,
) -> tuple[PolyBallLineModel, list[BranchGeometry], list[dict[str, Any]]]:
    if constructed_override is not None:
        constructed = constructed_override
        tail_rows = list(port_tail_rows_override or ())
    else:
        base = (
            v6_details.constructed_branches
            if v6_details is not None and v6_details.constructed_branches
            else branches
        )
        constructed, tail_rows = _copy_with_tail(base, ports, config)
    all_points = np.vstack([branch.points_um for branch in constructed])
    origin = np.mean(all_points, axis=0)
    if config.v7.oriented_grid and len(all_points) >= 3:
        _, _, axes = np.linalg.svd(all_points - origin, full_matrices=False)
        if np.linalg.det(axes) < 0.0:
            axes[-1] *= -1.0
    else:
        axes = np.eye(3)
    starts: list[np.ndarray] = []
    ends: list[np.ndarray] = []
    radius_start: list[float] = []
    radius_end: list[float] = []
    branch_ids: list[int] = []
    segment_indices: list[int] = []
    for branch in constructed:
        local = (np.asarray(branch.points_um) - origin) @ axes.T
        for index in range(len(local) - 1):
            if np.linalg.norm(local[index + 1] - local[index]) <= 1.0e-12:
                continue
            starts.append(local[index])
            ends.append(local[index + 1])
            radius_start.append(float(branch.radius_um[index]))
            radius_end.append(float(branch.radius_um[index + 1]))
            branch_ids.append(branch.branch_id)
            segment_indices.append(index)
    start_array = np.asarray(starts, dtype=float)
    end_array = np.asarray(ends, dtype=float)
    if not len(start_array):
        raise GeometryValidationError("V7 unified PolyBallLine has no valid segments")
    vmtk_available = bool(importlib.util.find_spec("vmtk"))
    requested = config.v7.polyball_provider
    if requested == "vmtk" and not vmtk_available:
        raise GeometryValidationError("V7 requested VMTK, but vmtk is unavailable")
    # The compatible VMTK module is preferred when present. The current Windows
    # environment cannot solve it against VTK 9.6, so the mathematically exact
    # line-envelope evaluator is the selected formal fallback.
    provider = "vtkvmtkPolyBallLine" if vmtk_available and requested != "ckdtree" else "ckdtree_exact_polyballline"
    midpoint = 0.5 * (start_array + end_array)
    return (
        PolyBallLineModel(
            segment_start_local=start_array,
            segment_end_local=end_array,
            radius_start_um=np.asarray(radius_start, dtype=float),
            radius_end_um=np.asarray(radius_end, dtype=float),
            segment_branch_id=np.asarray(branch_ids, dtype=np.int64),
            segment_index_in_branch=np.asarray(segment_indices, dtype=np.int64),
            origin_world_um=origin,
            local_axes=axes,
            tree=cKDTree(midpoint),
            provider=provider,
        ),
        constructed,
        tail_rows,
    )


def _grid_specification(
    model: PolyBallLineModel,
    config: CFDLumenConfig,
    cells_across_min_diameter: int,
) -> tuple[np.ndarray, tuple[int, int, int], float]:
    minimum_radius = float(
        min(model.radius_start_um.min(), model.radius_end_um.min())
    )
    spacing = 2.0 * minimum_radius / float(cells_across_min_diameter)
    maximum_radius = float(
        max(model.radius_start_um.max(), model.radius_end_um.max())
    )
    points = np.vstack((model.segment_start_local, model.segment_end_local))
    padding = maximum_radius + config.v7.field_padding_cells * spacing
    minimum = points.min(axis=0) - padding
    maximum = points.max(axis=0) + padding
    dimensions = tuple(
        int(np.ceil((maximum[axis] - minimum[axis]) / spacing)) + 1
        for axis in range(3)
    )
    cells = int(np.prod(dimensions, dtype=np.int64))
    if cells > config.v7.max_grid_cells:
        raise GeometryValidationError(
            f"V7_GRID_LIMIT_EXCEEDED: {cells} > {config.v7.max_grid_cells} cells "
            f"at {cells_across_min_diameter} cells/min-diameter"
        )
    return minimum, dimensions, spacing


def _allocate_field(
    dimensions_xyz: tuple[int, int, int],
    config: CFDLumenConfig,
    fill_value: float,
) -> tuple[np.ndarray, Path | None]:
    shape = (dimensions_xyz[2], dimensions_xyz[1], dimensions_xyz[0])
    cells = int(np.prod(dimensions_xyz, dtype=np.int64))
    if cells < config.v7.memory_map_threshold_cells:
        return np.full(shape, fill_value, dtype=np.float32), None
    backing_root = (
        Path(config.v7.memory_map_directory).resolve()
        if config.v7.memory_map_directory
        else Path(tempfile.gettempdir()).resolve()
    )
    backing_root.mkdir(parents=True, exist_ok=True)
    required_bytes = cells * np.dtype(np.float32).itemsize
    free_bytes = shutil.disk_usage(backing_root).free
    if free_bytes < int(1.25 * required_bytes):
        raise GeometryValidationError(
            f"V7_MEMMAP_SPACE_INSUFFICIENT: need {required_bytes} bytes for field, "
            f"free {free_bytes} bytes in {backing_root}"
        )
    descriptor, name = tempfile.mkstemp(
        prefix="ulm_v7_polyball_",
        suffix=".float32",
        dir=backing_root,
    )
    os.close(descriptor)
    path = Path(name)
    field = np.memmap(path, mode="w+", dtype=np.float32, shape=shape)
    field.fill(np.float32(fill_value))
    return field, path


def _rasterize_sparse_polyball(
    model: PolyBallLineModel,
    minimum: np.ndarray,
    dimensions_xyz: tuple[int, int, int],
    spacing: float,
    config: CFDLumenConfig,
) -> tuple[np.ndarray, Path | None, dict[str, Any]]:
    started = time.perf_counter()
    maximum_radius = float(
        max(model.radius_start_um.max(), model.radius_end_um.max())
    )
    fill_value = maximum_radius + (config.v7.field_padding_cells + 1) * spacing
    field, backing_path = _allocate_field(dimensions_xyz, config, fill_value)
    nx, ny, nz = dimensions_xyz
    evaluations = 0
    for segment_id in range(model.segment_count):
        start = model.segment_start_local[segment_id]
        end = model.segment_end_local[segment_id]
        radius0 = float(model.radius_start_um[segment_id])
        radius1 = float(model.radius_end_um[segment_id])
        band = max(radius0, radius1) + 2.0 * spacing
        lower = np.floor((np.minimum(start, end) - band - minimum) / spacing).astype(int)
        upper = np.ceil((np.maximum(start, end) + band - minimum) / spacing).astype(int) + 1
        lower = np.maximum(lower, 0)
        upper = np.minimum(upper, np.asarray((nx, ny, nz)))
        x_ids = np.arange(lower[0], upper[0])
        y_ids = np.arange(lower[1], upper[1])
        z_ids = np.arange(lower[2], upper[2])
        # Bound transient allocations while retaining vectorized exact segment evaluation.
        plane_points = max(1, len(x_ids) * len(y_ids))
        z_chunk = max(1, min(len(z_ids), 2_000_000 // plane_points))
        for offset in range(0, len(z_ids), z_chunk):
            selected_z = z_ids[offset : offset + z_chunk]
            zz, yy, xx = np.meshgrid(
                minimum[2] + selected_z * spacing,
                minimum[1] + y_ids * spacing,
                minimum[0] + x_ids * spacing,
                indexing="ij",
            )
            query = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            phi = _single_segment_phi(query, start, end, radius0, radius1).astype(
                np.float32
            )
            view = field[
                selected_z[0] : selected_z[-1] + 1,
                lower[1] : upper[1],
                lower[0] : upper[0],
            ]
            view[...] = np.minimum(view, phi.reshape(view.shape))
            evaluations += len(query)
    if isinstance(field, np.memmap):
        field.flush()
    if not float(np.min(field)) < 0.0 < float(np.max(field)):
        raise GeometryValidationError("V7 unified field does not bracket phi=0")
    return field, backing_path, {
        "sparse_segment_voxel_evaluations": evaluations,
        "field_runtime_s": time.perf_counter() - started,
        "field_storage": "numpy.memmap" if backing_path is not None else "RAM_float32",
        "field_backing_bytes": int(field.nbytes),
    }


def _rasterize_smooth_junction_polyball(
    model: SmoothJunctionPolyBallModel,
    minimum: np.ndarray,
    dimensions_xyz: tuple[int, int, int],
    spacing: float,
    config: CFDLumenConfig,
) -> tuple[np.ndarray, Path | None, dict[str, Any]]:
    """Rasterize hard PolyBall, replacing only compact junction-core voxels."""

    field, backing_path, report = _rasterize_sparse_polyball(
        model.hard_model, minimum, dimensions_xyz, spacing, config
    )
    _, smooth_report = _apply_smooth_junction_polyball(
        field, model, minimum, dimensions_xyz, spacing, config
    )
    report.update(smooth_report)
    return field, backing_path, report


def _apply_smooth_junction_polyball(
    field: np.ndarray,
    model: SmoothJunctionPolyBallModel,
    minimum: np.ndarray,
    dimensions_xyz: tuple[int, int, int],
    spacing: float,
    config: CFDLumenConfig,
) -> tuple[list[tuple[tuple[slice, slice, slice], np.ndarray]], dict[str, Any]]:
    """In-place compact-core replacement with exact backups for shared rasters."""

    started = time.perf_counter()
    nx, ny, nz = dimensions_xyz
    evaluations = 0
    backups: list[tuple[tuple[slice, slice, slice], np.ndarray]] = []
    for spec in model.junctions:
        center_local = model.hard_model.world_to_local(spec.center_world_um[None, :])[0]
        lower = np.floor(
            (center_local - spec.blend_length_um - minimum) / spacing
        ).astype(int)
        upper = (
            np.ceil((center_local + spec.blend_length_um - minimum) / spacing)
            .astype(int)
            + 1
        )
        lower = np.maximum(lower, 0)
        upper = np.minimum(upper, np.asarray((nx, ny, nz)))
        x_ids = np.arange(lower[0], upper[0])
        y_ids = np.arange(lower[1], upper[1])
        z_ids = np.arange(lower[2], upper[2])
        plane_points = max(1, len(x_ids) * len(y_ids))
        z_chunk = max(1, min(len(z_ids), 1_000_000 // plane_points))
        for offset in range(0, len(z_ids), z_chunk):
            selected_z = z_ids[offset : offset + z_chunk]
            zz, yy, xx = np.meshgrid(
                minimum[2] + selected_z * spacing,
                minimum[1] + y_ids * spacing,
                minimum[0] + x_ids * spacing,
                indexing="ij",
            )
            query_local = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            relative = query_local - center_local[None, :]
            inside = np.linalg.norm(relative, axis=1) < spec.blend_length_um
            if not np.any(inside):
                continue
            query_world = model.hard_model.local_to_world(query_local[inside])
            phi, _ = model.evaluate(
                query_world,
                k=config.v7.k_nearest_segments,
                gradients=False,
            )
            selection = (
                slice(int(selected_z[0]), int(selected_z[-1]) + 1),
                slice(int(lower[1]), int(upper[1])),
                slice(int(lower[0]), int(upper[0])),
            )
            view = field[selection]
            backups.append((selection, np.asarray(view).copy()))
            local_values = np.asarray(view, dtype=np.float32).reshape(-1).copy()
            local_values[inside] = phi.astype(np.float32)
            view[...] = local_values.reshape(view.shape)
            evaluations += int(np.count_nonzero(inside))
    if isinstance(field, np.memmap):
        field.flush()
    return backups, {
        "smooth_junction_voxel_evaluations": evaluations,
        "smooth_junction_runtime_s": time.perf_counter() - started,
        "smooth_union_equation": (
            "h=clip(0.5+0.5*(b-a)/k,0,1); "
            "smin=(1-h)*b+h*a-k*h*(1-h)"
        ),
        "smooth_union_reduction": "deterministic pairwise branch-id order",
        "compact_support": "quintic smootherstep; k=0 outside true junction cores",
        "competition_aware_support": (
            model.competition_threshold_radius_fraction is not None
        ),
        "competition_threshold_radius_fraction": (
            model.competition_threshold_radius_fraction
        ),
    }


def _restore_raster_backups(
    field: np.ndarray,
    backups: list[tuple[tuple[slice, slice, slice], np.ndarray]],
) -> None:
    for selection, values in reversed(backups):
        field[selection] = values
    if isinstance(field, np.memmap):
        field.flush()


def prepare_polyball_raster(
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    v6_details: HybridBuildDetails | None = None,
    cells_across_min_diameter: int | None = None,
    constructed_override: list[BranchGeometry] | None = None,
    port_tail_rows_override: list[dict[str, Any]] | None = None,
) -> PreparedPolyBallRaster:
    """Prepare one immutable hard-field raster for the three v8 k/r variants."""

    cells = int(cells_across_min_diameter or config.v7.cells_across_min_diameter)
    model, constructed, tail_rows = build_polyball_model(
        branches,
        ports,
        config,
        v6_details=v6_details,
        constructed_override=constructed_override,
        port_tail_rows_override=port_tail_rows_override,
    )
    minimum, dimensions, spacing = _grid_specification(model, config, cells)
    field, backing_path, field_report = _rasterize_sparse_polyball(
        model, minimum, dimensions, spacing, config
    )
    return PreparedPolyBallRaster(
        model=model,
        constructed_branches=constructed,
        port_tail_rows=tail_rows,
        minimum_local_um=minimum,
        dimensions_xyz=dimensions,
        spacing_um=spacing,
        field=field,
        backing_path=backing_path,
        field_report=field_report,
    )


def release_prepared_polyball_raster(prepared: PreparedPolyBallRaster) -> None:
    """Close and remove the formal temporary field created for v8 sensitivity."""

    memory_map = None
    if isinstance(prepared.field, np.memmap):
        prepared.field.flush()
        memory_map = prepared.field._mmap
    if memory_map is not None:
        memory_map.close()
    gc.collect()
    if prepared.backing_path is not None:
        try:
            prepared.backing_path.unlink()
        except PermissionError:
            pass


def _vtk_faces(polydata: vtk.vtkPolyData) -> np.ndarray:
    polygons = polydata.GetPolys()
    if hasattr(polygons, "GetOffsetsArray"):
        offsets = np.asarray(vtk_to_numpy(polygons.GetOffsetsArray()), dtype=np.int64)
        connectivity = np.asarray(
            vtk_to_numpy(polygons.GetConnectivityArray()), dtype=np.int64
        )
        if not np.all(np.diff(offsets) == 3):
            raise GeometryValidationError("V7 FlyingEdges output contains non-triangles")
        return connectivity.reshape((-1, 3))
    raw = vtk_to_numpy(polygons.GetData())
    return np.asarray(raw, dtype=np.int64).reshape((-1, 4))[:, 1:]


def _extract_flying_edges(
    field: np.ndarray,
    minimum: np.ndarray,
    spacing: float,
) -> trimesh.Trimesh:
    nz, ny, nx = field.shape
    image = vtk.vtkImageData()
    image.SetDimensions(nx, ny, nz)
    image.SetOrigin(*map(float, minimum))
    image.SetSpacing(spacing, spacing, spacing)
    scalars = numpy_to_vtk(
        np.asarray(field).reshape(-1, order="C"),
        deep=False,
        array_type=vtk.VTK_FLOAT,
    )
    scalars.SetName("polyball_phi_um")
    image.GetPointData().SetScalars(scalars)
    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputData(image)
    contour.SetValue(0, 0.0)
    contour.ComputeNormalsOff()
    contour.ComputeGradientsOff()
    contour.Update()
    output = contour.GetOutput()
    vertices = np.asarray(vtk_to_numpy(output.GetPoints().GetData()), dtype=float)
    faces = _vtk_faces(output)
    return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)


def _extract_marching_cubes(
    field: np.ndarray,
    minimum: np.ndarray,
    spacing: float,
) -> trimesh.Trimesh:
    vertices_zyx, faces, _, _ = marching_cubes(
        field, level=0.0, spacing=(spacing, spacing, spacing)
    )
    vertices = vertices_zyx[:, ::-1] + minimum
    return trimesh.Trimesh(vertices=vertices, faces=faces, process=False)


def _clean_mesh(mesh: trimesh.Trimesh) -> trimesh.Trimesh:
    clean = trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float),
        faces=np.asarray(mesh.faces, dtype=np.int64),
        process=True,
        validate=True,
    )
    clean.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(clean, multibody=True)
    if clean.volume < 0.0:
        clean.invert()
    return clean


def project_vertices_to_polyball(
    mesh: trimesh.Trimesh,
    model: PolyBallLineModel | SmoothJunctionPolyBallModel,
    config: CFDLumenConfig,
) -> tuple[trimesh.Trimesh, dict[str, Any]]:
    vertices = np.asarray(mesh.vertices, dtype=float).copy()
    before, _ = model.evaluate(
        vertices, k=config.v7.k_nearest_segments, gradients=False
    )
    iteration_rows: list[dict[str, Any]] = []
    for iteration in range(config.v7.newton_iterations):
        phi, gradient = model.evaluate(
            vertices, k=config.v7.k_nearest_segments, gradients=True
        )
        assert gradient is not None
        squared_norm = np.einsum("ij,ij->i", gradient, gradient)
        active = (np.abs(phi) > config.v7.projection_tolerance_um) & (
            squared_norm > 1.0e-16
        )
        vertices[active] -= (
            phi[active] / squared_norm[active]
        )[:, None] * gradient[active]
        iteration_rows.append(
            {
                "iteration": iteration + 1,
                "active_vertex_count": int(np.count_nonzero(active)),
                "phi_abs_p95_um": float(np.percentile(np.abs(phi), 95)),
                "phi_abs_max_um": float(np.max(np.abs(phi))),
            }
        )
        if not np.any(active):
            break
    after, _ = model.evaluate(
        vertices, k=config.v7.k_nearest_segments, gradients=False
    )
    projected = trimesh.Trimesh(
        vertices=vertices, faces=np.asarray(mesh.faces, dtype=np.int64), process=False
    )
    trimesh.repair.fix_normals(projected, multibody=True)
    minimum_radius = float(
        min(model.radius_start_um.min(), model.radius_end_um.min())
    )
    report = {
        "equation": "x_new = x - phi(x) * grad(phi) / ||grad(phi)||^2",
        "configured_iterations": config.v7.newton_iterations,
        "tolerance_um": config.v7.projection_tolerance_um,
        "pre_projection_phi_abs_mean_um": float(np.mean(np.abs(before))),
        "pre_projection_phi_abs_p95_um": float(np.percentile(np.abs(before), 95)),
        "pre_projection_phi_abs_max_um": float(np.max(np.abs(before))),
        "post_projection_phi_abs_mean_um": float(np.mean(np.abs(after))),
        "post_projection_phi_abs_p95_um": float(np.percentile(np.abs(after), 95)),
        "post_projection_phi_abs_max_um": float(np.max(np.abs(after))),
        "pre_projection_radius_error_p95": float(
            np.percentile(np.abs(before) / minimum_radius, 95)
        ),
        "post_projection_radius_error_p95": float(
            np.percentile(np.abs(after) / minimum_radius, 95)
        ),
        "iterations": iteration_rows,
    }
    return projected, report


def _to_pyvista(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def isotropic_remesh(
    mesh: trimesh.Trimesh,
    spacing: float,
    config: CFDLumenConfig,
) -> tuple[trimesh.Trimesh, dict[str, Any]]:
    before_quality = triangle_quality(mesh)["summary"]
    if config.v7.remesh_backend == "none":
        return mesh.copy(), {
            "backend": "none",
            "pre_triangle_quality": before_quality,
            "post_triangle_quality": before_quality,
        }
    try:
        import pyacvd
    except ImportError as exc:
        raise GeometryValidationError(
            "V7_REMesher_UNAVAILABLE: pyacvd is required for formal v7"
        ) from exc
    target_edge = config.v7.remesh_target_edge_spacing_factor * spacing
    estimated = int(mesh.area / max(0.8660254 * target_edge**2, 1.0e-12))
    target_clusters = min(
        max(20, len(mesh.vertices) // 2),
        max(min(config.v7.remesh_minimum_clusters, len(mesh.vertices) // 2), estimated),
    )
    attempt_rows: list[dict[str, Any]] = []
    output: trimesh.Trimesh | None = None
    selected_clusters: int | None = None
    attempts = [
        (int(target_clusters), False),
        (max(20, int(round(0.80 * target_clusters))), False),
        (max(20, int(round(0.65 * target_clusters))), False),
        (int(target_clusters), True),
    ]
    for attempt_id, (cluster_count, subdivide_input) in enumerate(attempts, start=1):
        cluster = pyacvd.Clustering(_to_pyvista(mesh))
        if subdivide_input:
            cluster.subdivide(1)
        cluster.cluster(cluster_count, maxiter=200, iso_try=50)
        remeshed = cluster.create_mesh(moveclus=True, flipnorm=True, clean=True).triangulate()
        faces = np.asarray(remeshed.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
        candidate = _clean_mesh(
            trimesh.Trimesh(
                vertices=np.asarray(remeshed.points, dtype=float),
                faces=faces,
                process=False,
            )
        )
        _, edge_counts = np.unique(
            np.asarray(candidate.edges_sorted), axis=0, return_counts=True
        )
        topology = {
            "attempt": attempt_id,
            "cluster_count": cluster_count,
            "subdivide_input_once": subdivide_input,
            "watertight": bool(candidate.is_watertight),
            "component_count": int(
                len(candidate.split(only_watertight=False))
            ),
            "boundary_edge_count": int(np.count_nonzero(edge_counts == 1)),
            "nonmanifold_edge_count": int(np.count_nonzero(edge_counts > 2)),
            "triangle_count": int(len(candidate.faces)),
        }
        attempt_rows.append(topology)
        if (
            topology["watertight"]
            and topology["component_count"] == 1
            and topology["boundary_edge_count"] == 0
            and topology["nonmanifold_edge_count"] == 0
        ):
            output = candidate
            selected_clusters = cluster_count
            break
    if output is None:
        raise GeometryValidationError(
            "V7_REMESH_TOPOLOGY_FAILED after deterministic PyACVD attempts: "
            + repr(attempt_rows)
        )
    return output, {
        "backend": "pyacvd_isotropic_clustering",
        "geometry_smoothing": False,
        "target_edge_um": target_edge,
        "target_cluster_count": int(target_clusters),
        "selected_cluster_count": selected_clusters,
        "deterministic_topology_attempts": attempt_rows,
        "pre_vertex_count": int(len(mesh.vertices)),
        "pre_triangle_count": int(len(mesh.faces)),
        "post_vertex_count": int(len(output.vertices)),
        "post_triangle_count": int(len(output.faces)),
        "pre_triangle_quality": before_quality,
        "post_triangle_quality": triangle_quality(output)["summary"],
    }


def _boundary_loops(mesh: trimesh.Trimesh) -> list[np.ndarray]:
    inverse = np.asarray(mesh.edges_unique_inverse, dtype=np.int64)
    counts = np.bincount(inverse, minlength=len(mesh.edges_unique))
    boundary = np.asarray(mesh.edges_unique[counts == 1], dtype=np.int64)
    adjacency: dict[int, list[int]] = {}
    for first, second in boundary:
        adjacency.setdefault(int(first), []).append(int(second))
        adjacency.setdefault(int(second), []).append(int(first))
    loops: list[np.ndarray] = []
    unused = {tuple(sorted(map(int, edge))) for edge in boundary}
    while unused:
        seed = next(iter(unused))
        ordered = [seed[0], seed[1]]
        unused.remove(seed)
        previous, current = seed
        while current != ordered[0]:
            choices = [
                node
                for node in adjacency[current]
                if tuple(sorted((current, node))) in unused and node != previous
            ]
            if not choices:
                break
            following = choices[0]
            unused.remove(tuple(sorted((current, following))))
            if following != ordered[0]:
                ordered.append(following)
            previous, current = current, following
        if current != ordered[0] or len(ordered) < 3:
            raise GeometryValidationError("V7_PORT_CLIP produced an open/non-simple boundary")
        loops.append(np.asarray(ordered, dtype=np.int64))
    return loops


def _clip_local_port_tail(
    mesh: trimesh.Trimesh,
    port: PortGeometry,
) -> trimesh.Trimesh:
    """Clip only the outlet tail, not the complete mesh by an infinite halfspace."""

    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    tangent = np.asarray(port.outward_tangent, dtype=float)
    relative = vertices - np.asarray(port.cap_center_um, dtype=float)[None, :]
    signed = relative @ tangent
    radial = np.linalg.norm(relative - signed[:, None] * tangent[None, :], axis=1)
    output_vertices = vertices.tolist()
    output_faces: list[tuple[int, int, int]] = []
    intersection_by_edge: dict[tuple[int, int], int] = {}

    def intersection(first: int, second: int) -> int:
        key = tuple(sorted((int(first), int(second))))
        existing = intersection_by_edge.get(key)
        if existing is not None:
            return existing
        first_signed = float(signed[first])
        second_signed = float(signed[second])
        fraction = first_signed / (first_signed - second_signed)
        point = vertices[first] + fraction * (vertices[second] - vertices[first])
        # Remove roundoff in the normal coordinate so the port is exactly planar.
        point -= float(np.dot(point - port.cap_center_um, tangent)) * tangent
        vertex_id = len(output_vertices)
        output_vertices.append(point.tolist())
        intersection_by_edge[key] = vertex_id
        return vertex_id

    for face in faces:
        face_signed = signed[face]
        face_radial = radial[face]
        local_tail = (
            float(np.min(face_radial)) <= 2.0 * port.radius_um
            and float(np.max(face_signed)) >= -2.0 * port.radius_um
            and float(np.min(face_signed)) <= 4.0 * port.radius_um
        )
        if not local_tail or np.all(face_signed <= 0.0):
            output_faces.append(tuple(map(int, face)))
            continue
        if np.all(face_signed > 0.0):
            continue
        polygon: list[int] = []
        previous = int(face[-1])
        previous_inside = bool(signed[previous] <= 0.0)
        for current_raw in face:
            current = int(current_raw)
            current_inside = bool(signed[current] <= 0.0)
            if current_inside:
                if not previous_inside:
                    polygon.append(intersection(previous, current))
                polygon.append(current)
            elif previous_inside:
                polygon.append(intersection(previous, current))
            previous = current
            previous_inside = current_inside
        if len(polygon) < 3:
            continue
        for index in range(1, len(polygon) - 1):
            output_faces.append((polygon[0], polygon[index], polygon[index + 1]))
    clipped = trimesh.Trimesh(
        vertices=np.asarray(output_vertices, dtype=float),
        faces=np.asarray(output_faces, dtype=np.int64),
        process=False,
    )
    clipped.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(clipped, multibody=True)
    return clipped


def clip_and_cap_ports(
    wall: trimesh.Trimesh,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> tuple[trimesh.Trimesh, np.ndarray, list[dict[str, Any]]]:
    clipped = wall.copy()
    clip_rows: list[dict[str, Any]] = []
    for port in ports:
        candidate = _clip_local_port_tail(clipped, port)
        if not len(candidate.faces):
            raise GeometryValidationError(f"V7_PORT_CLIP_FAILED: port {port.port_id}")
        clipped = candidate
        clip_rows.append(
            {
                "port_id": port.port_id,
                "cut_port_id": port.cut_port_id,
                "plane_center_um": np.asarray(port.cap_center_um).tolist(),
                "plane_outward_normal": np.asarray(port.outward_tangent).tolist(),
                "exact_plane_clip": True,
                "clip_scope": "local implicit port tail only",
                "infinite_global_halfspace_clip": False,
            }
        )
    loops = _boundary_loops(clipped)
    if len(loops) != len(ports):
        raise GeometryValidationError(
            f"V7_PORT_CLIP produced {len(loops)} boundary loops for {len(ports)} ports"
        )
    assignment_cost = np.empty((len(ports), len(loops)), dtype=float)
    assignment_metrics: dict[tuple[int, int], tuple[float, float, float]] = {}
    for port_index, port in enumerate(ports):
        for loop_id, loop in enumerate(loops):
            points = np.asarray(clipped.vertices)[loop]
            relative = points - np.asarray(port.cap_center_um)[None, :]
            axial = relative @ np.asarray(port.outward_tangent)
            planar = relative - axial[:, None] * np.asarray(port.outward_tangent)[None, :]
            radial = np.linalg.norm(planar, axis=1)
            axial_error = float(np.mean(np.abs(axial)))
            radial_p95 = float(np.percentile(radial, 95))
            center_offset = float(np.linalg.norm(np.mean(planar, axis=0)))
            # Exact-plane residual dominates. Radial excess disqualifies a
            # remote loop, while the unweighted vertex centroid has only a
            # small tie-breaking weight because ACVD sampling is nonuniform.
            assignment_cost[port_index, loop_id] = (
                axial_error / max(port.radius_um, 1.0e-12)
                + 10.0 * max(0.0, radial_p95 / port.radius_um - 1.25)
                + 0.05 * center_offset / max(port.radius_um, 1.0e-12)
            )
            assignment_metrics[(port_index, loop_id)] = (
                axial_error,
                radial_p95,
                center_offset,
            )
    assigned_ports, assigned_loops = linear_sum_assignment(assignment_cost)
    loop_by_port = {
        int(port_index): int(loop_id)
        for port_index, loop_id in zip(assigned_ports, assigned_loops)
    }
    vertices = np.asarray(clipped.vertices, dtype=float).tolist()
    faces = np.asarray(clipped.faces, dtype=np.int64).tolist()
    boundary_ids = [0] * len(faces)
    for port_index, port in enumerate(ports):
        loop_id = loop_by_port[port_index]
        axial_error, radial_p95, center_offset = assignment_metrics[
            (port_index, loop_id)
        ]
        plane_tolerance = max(config.ports.plane_tolerance_um, 1.0e-6)
        if axial_error > plane_tolerance or radial_p95 > 1.5 * port.radius_um:
            raise GeometryValidationError(
                f"V7_PORT_CAP_LOOP_MISMATCH: port {port.port_id}, axial "
                f"{axial_error:.6g} um, radial P95 {radial_p95:.6g} um, "
                f"center offset {center_offset:.6g} um"
            )
        loop = loops[loop_id]
        center_id = len(vertices)
        vertices.append(np.asarray(port.cap_center_um, dtype=float).tolist())
        for index, first in enumerate(loop):
            second = int(loop[(index + 1) % len(loop)])
            first = int(first)
            normal = np.cross(
                np.asarray(vertices[second]) - np.asarray(vertices[first]),
                np.asarray(vertices[center_id]) - np.asarray(vertices[first]),
            )
            if float(np.dot(normal, port.outward_tangent)) < 0.0:
                first, second = second, first
            faces.append((first, second, center_id))
            boundary_ids.append(port.port_id + 1)
        clip_rows[port.port_id].update(
            {
                "boundary_loop_vertex_count": int(len(loop)),
                "boundary_loop_plane_residual_um": axial_error,
                "boundary_loop_radial_p95_um": radial_p95,
                "boundary_loop_vertex_mean_center_offset_um": center_offset,
                "boundary_loop_global_assignment_cost": float(
                    assignment_cost[port_index, loop_id]
                ),
                "flat_cap_triangle_count": int(len(loop)),
                "boundary_id": port.port_id + 1,
            }
        )
    capped = trimesh.Trimesh(
        vertices=np.asarray(vertices, dtype=float),
        faces=np.asarray(faces, dtype=np.int64),
        process=False,
    )
    trimesh.repair.fix_normals(capped, multibody=True)
    if capped.volume < 0.0:
        capped.invert()
        # Inversion reverses winding only; face order and IDs are unchanged.
    if not capped.is_watertight:
        raise GeometryValidationError("V7 flat-port capping did not produce a watertight mesh")
    return capped, np.asarray(boundary_ids, dtype=np.int32), clip_rows


def _patch_from_boundary_ids(
    mesh: trimesh.Trimesh,
    boundary_ids: np.ndarray,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> PatchResult:
    from .surface_qc import identify_port_patches

    measured = identify_port_patches(mesh, ports, config)
    patch_id = np.asarray(boundary_ids, dtype=np.int32)
    patch_type = (patch_id > 0).astype(np.int8)
    port_id = np.where(patch_id > 0, patch_id - 1, -1).astype(np.int32)
    return PatchResult(
        patch_id=patch_id,
        patch_type=patch_type,
        port_id=port_id,
        port_rows=measured.port_rows,
        detected_port_count=measured.detected_port_count,
        all_ports_pass=measured.all_ports_pass,
    )


def build_unified_polyball_surface(
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    v6_details: HybridBuildDetails | None = None,
    cells_across_min_diameter: int | None = None,
    compare_extractors: bool | None = None,
    remesh: bool = True,
    junction_specs: tuple[JunctionBlendSpec, ...] = (),
    smooth_k_radius_ratio: float | None = None,
    prepared_raster: PreparedPolyBallRaster | None = None,
    constructed_override: list[BranchGeometry] | None = None,
    port_tail_rows_override: list[dict[str, Any]] | None = None,
    competition_threshold_radius_fraction: float | None = None,
) -> UnifiedPolyBallBuild:
    """Build one complete implicit wall, then remesh/project and cut final ports."""

    started = time.perf_counter()
    cells = int(cells_across_min_diameter or config.v7.cells_across_min_diameter)
    if prepared_raster is None:
        model, constructed, tail_rows = build_polyball_model(
            branches,
            ports,
            config,
            v6_details=v6_details,
            constructed_override=constructed_override,
            port_tail_rows_override=port_tail_rows_override,
        )
        minimum, dimensions, spacing = _grid_specification(model, config, cells)
    else:
        model = prepared_raster.model
        constructed = prepared_raster.constructed_branches
        tail_rows = prepared_raster.port_tail_rows
        minimum = prepared_raster.minimum_local_um
        dimensions = prepared_raster.dimensions_xyz
        spacing = prepared_raster.spacing_um
        expected_cells = int(
            2.0
            * min(model.radius_start_um.min(), model.radius_end_um.min())
            / spacing
        )
        if abs(expected_cells - cells) > 1:
            raise GeometryValidationError(
                "V8 prepared hard raster resolution does not match requested cells"
            )
    field_model: PolyBallLineModel | SmoothJunctionPolyBallModel = model
    if smooth_k_radius_ratio is not None:
        if smooth_k_radius_ratio <= 0.0 or not junction_specs:
            raise GeometryValidationError(
                "V8 smooth union requires positive k/r and true junction specifications"
            )
        field_model = SmoothJunctionPolyBallModel(
            hard_model=model,
            junctions=tuple(junction_specs),
            k_radius_ratio=float(smooth_k_radius_ratio),
            competition_threshold_radius_fraction=(
                float(competition_threshold_radius_fraction)
                if competition_threshold_radius_fraction is not None
                else None
            ),
        )
    field: np.ndarray | None = None
    backing_path: Path | None = None
    raster_backups: list[tuple[tuple[slice, slice, slice], np.ndarray]] = []
    try:
        if prepared_raster is not None:
            field = prepared_raster.field
            backing_path = prepared_raster.backing_path
            field_report = dict(prepared_raster.field_report)
            if isinstance(field_model, SmoothJunctionPolyBallModel):
                raster_backups, smooth_report = _apply_smooth_junction_polyball(
                    field, field_model, minimum, dimensions, spacing, config
                )
                field_report.update(smooth_report)
        elif isinstance(field_model, SmoothJunctionPolyBallModel):
            field, backing_path, field_report = _rasterize_smooth_junction_polyball(
                field_model, minimum, dimensions, spacing, config
            )
        else:
            field, backing_path, field_report = _rasterize_sparse_polyball(
                model, minimum, dimensions, spacing, config
            )
        extraction_rows: list[dict[str, Any]] = []
        methods = [config.v7.extraction_backend]
        compare = config.v7.compare_marching_cubes if compare_extractors is None else compare_extractors
        if compare:
            alternate = (
                "marching_cubes"
                if config.v7.extraction_backend == "flying_edges"
                else "flying_edges"
            )
            methods.append(alternate)
        candidates: dict[str, trimesh.Trimesh] = {}
        for method in methods:
            extraction_started = time.perf_counter()
            local_mesh = (
                _extract_flying_edges(field, minimum, spacing)
                if method == "flying_edges"
                else _extract_marching_cubes(field, minimum, spacing)
            )
            local_mesh = _clean_mesh(local_mesh)
            world_mesh = local_mesh.copy()
            world_mesh.vertices = model.local_to_world(local_mesh.vertices)
            trimesh.repair.fix_normals(world_mesh, multibody=True)
            candidates[method] = world_mesh
            extraction_rows.append(
                {
                    "backend": method,
                    "runtime_s": time.perf_counter() - extraction_started,
                    "vertex_count": int(len(world_mesh.vertices)),
                    "triangle_count": int(len(world_mesh.faces)),
                    "watertight": bool(world_mesh.is_watertight),
                    "component_count": int(
                        len(world_mesh.split(only_watertight=False))
                    ),
                    "triangle_quality": triangle_quality(world_mesh)["summary"],
                }
            )
        raw = candidates[config.v7.extraction_backend]
    finally:
        if prepared_raster is not None:
            if field is not None and raster_backups:
                _restore_raster_backups(field, raster_backups)
            field = None
        else:
            memory_map = None
            if isinstance(field, np.memmap):
                field.flush()
                memory_map = field._mmap
            field = None
            if backing_path is not None:
                gc.collect()
                if memory_map is not None:
                    memory_map.close()
                try:
                    backing_path.unlink()
                except PermissionError:
                    pass
    if not raw.is_watertight or len(raw.split(only_watertight=False)) != 1:
        raise GeometryValidationError(
            "V7_UNIFIED_ISOSURFACE_TOPOLOGY_FAILED before port clipping"
        )
    projected_raw, first_projection = project_vertices_to_polyball(
        raw, field_model, config
    )
    remesh_started = time.perf_counter()
    if remesh:
        remeshed, remesh_report = isotropic_remesh(projected_raw, spacing, config)
    else:
        remeshed, remesh_report = projected_raw, {
            "backend": "disabled_for_convergence",
            "pre_triangle_quality": triangle_quality(projected_raw)["summary"],
            "post_triangle_quality": triangle_quality(projected_raw)["summary"],
        }
    remesh_report["runtime_s"] = time.perf_counter() - remesh_started
    projected_remesh, second_projection = project_vertices_to_polyball(
        remeshed, field_model, config
    )
    capped, boundary_ids, clip_rows = clip_and_cap_ports(
        projected_remesh, ports, config
    )
    patch = _patch_from_boundary_ids(capped, boundary_ids, ports, config)
    is_smooth = isinstance(field_model, SmoothJunctionPolyBallModel)
    is_v9 = bool(
        is_smooth
        and field_model.competition_threshold_radius_fraction is not None
    )
    is_v8 = is_smooth
    metadata = {
        "protocol": "v8 local smooth-union diagnosis" if is_v8 else "(new) 子图建模修改v7",
        "backend": (
            "c1_spline_polyball_competition_union"
            if is_v9
            else "unified_polyball_smooth_junction"
            if is_smooth
            else "unified_polyball"
        ),
        "polyball_provider": model.provider,
        "vmtk_available": bool(importlib.util.find_spec("vmtk")),
        "vmtk_environment_decision": (
            "available"
            if importlib.util.find_spec("vmtk")
            else "conda-forge win-64 VMTK requires obsolete VTK 9.1/9.2; preserve VTK 9.6 environment and use exact cKDTree fallback"
        ),
        "source_swc_modified": False,
        "source_radius_modified": False,
        "single_continuous_implicit_field": True,
        "one_iso_surface_extraction": True,
        "vtk_tube_filter_surface_count": 0,
        "surface_stitch_count": 0,
        "surface_boolean_count": 0,
        "collar_surface_count": 0,
        "hybrid_interface_edge_count": 0,
        "projection_field": (
            "Phi_v9_competition_aware"
            if is_v9
            else "Phi_v8"
            if is_smooth
            else "Phi_v7_hard_min"
        ),
        "same_field_for_flying_edges_and_both_newton_projections": True,
        "smooth_k_radius_ratio": smooth_k_radius_ratio,
        "competition_threshold_radius_fraction": competition_threshold_radius_fraction,
        "competition_aware_support": is_v9,
        "hard_field_reused_for_k_sensitivity": prepared_raster is not None,
        "junction_blends": [
            {
                "junction_node_id": spec.junction_node_id,
                "center_world_um": np.asarray(spec.center_world_um).tolist(),
                "junction_radius_um": spec.radius_um,
                "blend_length_um": spec.blend_length_um,
                "incident_branch_ids": list(spec.incident_branch_ids),
            }
            for spec in junction_specs
        ],
        "port_tail_rows": tail_rows,
        "port_clip_rows": clip_rows,
        "segment_count": model.segment_count,
        "cells_across_min_diameter": cells,
        "grid_spacing_um": spacing,
        "grid_dimensions_xyz": list(dimensions),
        "grid_cell_count": int(np.prod(dimensions, dtype=np.int64)),
        "grid_coordinate_system": "PCA-oriented physical Cartesian grid" if config.v7.oriented_grid else "world Cartesian grid",
        "grid_origin_local_um": minimum.tolist(),
        "local_axes_rows": model.local_axes.tolist(),
        **field_report,
        "extractor_comparison": extraction_rows,
        "selected_extractor": config.v7.extraction_backend,
        "projection_before_remesh": first_projection,
        "isotropic_remesh": remesh_report,
        "projection_after_remesh": second_projection,
        "boundary_ids": {"WALL": 0, **{f"PORT_{port.port_id}": port.port_id + 1 for port in ports}},
        "runtime_s": time.perf_counter() - started,
    }
    if is_v9:
        metadata["protocol"] = "v9 C1 spline plus competition-aware junction union"
    return UnifiedPolyBallBuild(
        mesh=capped,
        wall_mesh_before_clip=projected_remesh,
        patch=patch,
        face_boundary_id=boundary_ids,
        model=model,
        field_model=field_model,
        constructed_branches=constructed,
        stage_meshes={
            "S0_raw_flying_edges": raw.copy(),
            "S1_newton_projected": projected_raw.copy(),
            "S2_pyacvd_before_second_projection": remeshed.copy(),
            "S3_final_projected_before_port_clip": projected_remesh.copy(),
        },
        metadata=metadata,
    )
