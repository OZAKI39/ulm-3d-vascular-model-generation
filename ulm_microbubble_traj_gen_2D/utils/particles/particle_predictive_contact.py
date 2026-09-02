"""Revised-v15 predictive rigid-wall transaction for particle transport.

The wall geometry in this module is exclusively the finite authoritative
solid-face segment set.  Contact is predicted before a particle is moved and
the single-wall low-Re mobility constraint is solved in closed form.  A trial
whose complete centre chord cannot be certified is rejected so the caller can
retry two genuine half-time intervals.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from .particle_geometry_tolerance import (
    PARTICLE_GEOMETRY_ROUNDOFF_FACTOR,
    particle_geometry_roundoff_tolerance_um,
)

_STATUS_OK = np.int8(0)
_STATUS_START_PENETRATION = np.int8(1)
_STATUS_INVALID_MOBILITY = np.int8(2)
_STATUS_ENDPOINT_PENETRATION = np.int8(3)
_STATUS_SWEPT_PATH_PENETRATION = np.int8(4)


@dataclass(frozen=True, slots=True)
class PredictiveContactStep:
    """Result of one speculative v15 particle-position transaction."""

    accepted_position_grid: np.ndarray
    angle_increment_rad: np.ndarray
    constrained_generalized_velocity: np.ndarray
    reaction_force_pn: np.ndarray
    endpoint_gap_um: np.ndarray
    minimum_path_gap_um: np.ndarray
    predicted_normal_gap_um: np.ndarray
    complementarity_residual_pn_um: np.ndarray
    residual_correction_um: np.ndarray
    contact_normal_xz: np.ndarray
    active_contact: np.ndarray
    outlet_fraction: np.ndarray
    outlet_position_grid: np.ndarray
    outlet_label: np.ndarray
    need_time_refinement: np.ndarray
    failure_codes: np.ndarray
    failure_position_grid: np.ndarray
    failure_event_fraction: np.ndarray

    @property
    def earliest_outlet_fraction(self) -> float | None:
        finite = self.outlet_fraction[np.isfinite(self.outlet_fraction)]
        return float(np.min(finite)) if finite.size else None


@dataclass(frozen=True, slots=True)
class PredictiveCandidateWall:
    solid_face_index: int
    center_xz_um: tuple[float, float]
    length_um: float
    inward_normal_xz: tuple[float, float]
    exact_distance_um: float
    exact_gap_um: float


@dataclass(frozen=True, slots=True)
class PredictiveWallContactFailure:
    local_lane: int
    event_fraction: float
    position_grid: tuple[float, float]
    position_xz_um: tuple[float, float]
    gap_um: float
    bubble_radius_um: float
    candidate_walls: tuple[PredictiveCandidateWall, ...]
    candidate_wall_source: str
    multiple_wall_contact: bool
    reason: str


class PredictiveWallContactGeometryError(RuntimeError):
    """A v15 contact cannot be represented by the single-wall model."""

    def __init__(
        self,
        message: str,
        failures: tuple[PredictiveWallContactFailure, ...],
    ) -> None:
        self.failures = failures
        lanes = [failure.local_lane for failure in failures]
        super().__init__(f"{message} local_lanes={lanes}.")


@dataclass(frozen=True, slots=True)
class PredictiveBoundaryExitFailure:
    """One trial that left the continuous lumen through an unclaimed opening."""

    local_lane: int
    start_position_grid: tuple[float, float]
    end_position_grid: tuple[float, float]
    start_position_xz_um: tuple[float, float]
    end_position_xz_um: tuple[float, float]
    bubble_radius_um: float
    nearest_open_section_index: int
    nearest_open_section_kind: int
    nearest_open_section_label: int
    section_signed_start_um: float
    section_signed_end_um: float
    section_lateral_end_um: float
    section_half_width_um: float


class PredictiveBoundaryLifecycleError(RuntimeError):
    """A trial crossed an open boundary that the lifecycle did not claim."""

    def __init__(
        self,
        failures: tuple[PredictiveBoundaryExitFailure, ...],
    ) -> None:
        self.failures = failures
        lanes = [failure.local_lane for failure in failures]
        super().__init__(
            "A v16 particle path ended outside the continuous lumen without "
            f"crossing an authoritative outlet; local_lanes={lanes}."
        )


def solve_predictive_contact_step(
    start_positions_grid: np.ndarray,
    alive: np.ndarray,
    dt_s: float,
    free_generalized_velocity: np.ndarray,
    generalized_mobility: np.ndarray,
    bubble_radii_um: np.ndarray,
    *,
    grid_spacing_um: float,
    boundary_geometry: ContinuousVesselGeometry,
    geometry_tolerance_um: float,
    use_numba: bool,
    inputs_prevalidated: bool = False,
    start_wall_distance_um: np.ndarray | None = None,
    start_wall_normal_xz: np.ndarray | None = None,
) -> PredictiveContactStep:
    """Predict, constrain, and certify one complete v15 batch trial.

    For one active wall with generalized normal ``a=(nx,nz,0)``, the reaction
    is the closed-form unilateral solution

    ``lambda=max(0, -(g/dt + a.T@U0)/(a.T@M@a))``.

    Exact endpoint and swept-disc checks are performed after this local planar
    solve.  A nonlocal/curved-wall inconsistency requests physical-time
    refinement; it is never hidden by shortening a displacement while still
    consuming the original interval.
    """

    positions = np.ascontiguousarray(start_positions_grid, dtype=np.float64)
    live = np.ascontiguousarray(alive, dtype=bool)
    free_velocity = np.ascontiguousarray(
        free_generalized_velocity, dtype=np.float64
    )
    mobility = np.ascontiguousarray(generalized_mobility, dtype=np.float64)
    radii = np.ascontiguousarray(bubble_radii_um, dtype=np.float64)
    count = int(positions.shape[0])
    spacing = float(grid_spacing_um)
    dt = float(dt_s)
    tolerance = float(geometry_tolerance_um)
    if not inputs_prevalidated:
        _validate_inputs(
            positions,
            live,
            dt,
            free_velocity,
            mobility,
            radii,
            spacing,
            tolerance,
            boundary_geometry,
        )
    cached_start_distance = None
    cached_start_normal = None
    if start_wall_distance_um is not None or start_wall_normal_xz is not None:
        if start_wall_distance_um is None or start_wall_normal_xz is None:
            raise ValueError(
                "Cached wall distance and wall normal must be supplied together."
            )
        cached_start_distance = np.asarray(
            start_wall_distance_um, dtype=np.float64
        ).reshape(-1)
        cached_start_normal = np.asarray(
            start_wall_normal_xz, dtype=np.float64
        ).reshape(-1, 2)
        if (
            cached_start_distance.shape != (count,)
            or cached_start_normal.shape != (count, 2)
        ):
            raise ValueError(
                "Cached wall geometry must match the predictive-contact batch."
            )

    accepted_position = positions.copy()
    angle_increment = np.zeros(count, dtype=np.float64)
    constrained_velocity = np.zeros((count, 3), dtype=np.float64)
    reaction = np.zeros(count, dtype=np.float64)
    endpoint_gap = np.full(count, np.nan, dtype=np.float64)
    minimum_path_gap = np.full(count, np.nan, dtype=np.float64)
    predicted_gap = np.full(count, np.nan, dtype=np.float64)
    complementarity = np.zeros(count, dtype=np.float64)
    residual_correction = np.zeros(count, dtype=np.float64)
    contact_normal = np.zeros((count, 2), dtype=np.float64)
    active_contact = np.zeros(count, dtype=bool)
    outlet_fraction = np.full(count, np.nan, dtype=np.float64)
    outlet_position = np.full((count, 2), np.nan, dtype=np.float64)
    outlet_label = np.full(count, -1, dtype=np.int32)
    need_refinement = np.zeros(count, dtype=bool)
    failure_codes = np.zeros(count, dtype=np.int8)
    failure_position = np.full((count, 2), np.nan, dtype=np.float64)
    failure_fraction = np.full(count, np.nan, dtype=np.float64)

    world_start = np.asarray(
        boundary_geometry.grid_to_world_xz(positions), dtype=np.float64
    )
    free_world_end = world_start + dt * free_velocity[:, :2]
    live_lanes = np.flatnonzero(live)
    start_distance = np.full(count, np.nan, dtype=np.float64)
    start_normal = np.zeros((count, 2), dtype=np.float64)
    free_endpoint_gap_all = np.full(count, np.nan, dtype=np.float64)
    free_path_gap = np.full(count, np.nan, dtype=np.float64)
    free_path_fraction = np.full(count, np.nan, dtype=np.float64)
    free_path_multiple = np.zeros(count, dtype=bool)
    if live_lanes.size:
        if cached_start_distance is None:
            start_batch = (
                boundary_geometry.exact_solid_wall_state_xz_um_accelerated(
                    world_start[live_lanes]
                )
                if use_numba
                else boundary_geometry.exact_solid_wall_state_xz_um(
                    world_start[live_lanes]
                )
            )
            start_distance[live_lanes] = np.asarray(
                start_batch.distance_um, dtype=np.float64
            ).reshape(-1)
            start_normal[live_lanes] = np.asarray(
                start_batch.inward_normal_xz, dtype=np.float64
            ).reshape(-1, 2)
        else:
            start_distance[live_lanes] = cached_start_distance[live_lanes]
            start_normal[live_lanes] = cached_start_normal[live_lanes]
        if use_numba:
            fused = (
                boundary_geometry.inspect_swept_solid_wall_paths_with_end_distance_xz_um_precomputed(
                world_start[live_lanes],
                free_world_end[live_lanes],
                radii[live_lanes],
                start_distance[live_lanes],
                tolerance_um=0.0,
            )
            )
            end_distance = np.asarray(fused[0], dtype=np.float64).reshape(-1)
            gap = np.asarray(fused[1], dtype=np.float64).reshape(-1)
            fraction = np.asarray(fused[2], dtype=np.float64).reshape(-1)
            multiple = np.asarray(fused[4], dtype=bool).reshape(-1)
            free_endpoint_gap_all[live_lanes] = end_distance - radii[live_lanes]
            free_path_gap[live_lanes] = gap
            free_path_fraction[live_lanes] = fraction
            free_path_multiple[live_lanes] = multiple
        else:
            end_batch = boundary_geometry.exact_solid_wall_state_xz_um(
                free_world_end[live_lanes]
            )
            end_distance = np.asarray(
                end_batch.distance_um, dtype=np.float64
            ).reshape(-1)
            free_endpoint_gap_all[live_lanes] = end_distance - radii[live_lanes]
            for lane_value in live_lanes:
                lane = int(lane_value)
                strict = _strict_geometry_epsilon(
                    world_start[lane], free_world_end[lane], radii[lane], spacing
                )
                inspected = boundary_geometry.inspect_swept_solid_wall_path_xz_um(
                    world_start[lane],
                    free_world_end[lane],
                    float(radii[lane]),
                    tolerance_um=strict,
                )
                free_path_gap[lane] = float(inspected.minimum_gap_um)
                if inspected.first_contact_fraction is not None:
                    free_path_fraction[lane] = float(
                        inspected.first_contact_fraction
                    )
                free_path_multiple[lane] = bool(
                    inspected.multiple_wall_contact
                )
    free_outlet_fraction, free_outlet_position, free_outlet_label = (
        _first_outlet_crossing_data_grid(
            boundary_geometry,
            positions[live_lanes],
            boundary_geometry.world_xz_to_grid(free_world_end[live_lanes]),
            use_numba=bool(use_numba),
        )
    )
    fraction_tolerance = _fraction_tolerance(dt)
    strict_epsilon_all = np.full(count, np.nan, dtype=np.float64)
    if live_lanes.size:
        coordinate_scale = np.maximum.reduce(
            (
                np.max(np.abs(world_start), axis=1),
                np.max(np.abs(free_world_end), axis=1),
                np.abs(radii),
                np.full(count, abs(spacing), dtype=np.float64),
                np.ones(count, dtype=np.float64),
            )
        )
        strict_epsilon_all = (
            PARTICLE_GEOMETRY_ROUNDOFF_FACTOR
            * np.finfo(np.float64).eps
            * coordinate_scale
        )
    start_gap_all = start_distance - radii
    invalid_start = live & (start_gap_all < -strict_epsilon_all)
    if np.any(invalid_start):
        lane = int(np.flatnonzero(invalid_start)[0])
        raise PredictiveWallContactGeometryError(
            "A v15 trial started from an already penetrating accepted state;",
            (
                _wall_contact_failure(
                    local_lane=lane,
                    event_fraction=0.0,
                    position_xz_um=world_start[lane],
                    bubble_radius_um=float(radii[lane]),
                    boundary_geometry=boundary_geometry,
                    geometry_tolerance_um=tolerance,
                    reason="accepted_start_has_negative_exact_gap",
                    known_gap_um=float(start_gap_all[lane]),
                ),
            ),
        )

    if (
        live_lanes.size
        and not np.any(np.isfinite(free_outlet_fraction))
        and np.all(
            free_path_gap[live_lanes]
            >= -strict_epsilon_all[live_lanes]
        )
        and np.all(free_endpoint_gap_all[live_lanes] >= 0.0)
    ):
        lanes = live_lanes
        accepted_position[lanes] = boundary_geometry.world_xz_to_grid(
            free_world_end[lanes]
        )
        angle_increment[lanes] = dt * free_velocity[lanes, 2]
        constrained_velocity[lanes] = free_velocity[lanes]
        endpoint_gap[lanes] = free_endpoint_gap_all[lanes]
        minimum_path_gap[lanes] = np.maximum(free_path_gap[lanes], 0.0)
        predicted_gap[lanes] = start_gap_all[lanes] + dt * np.sum(
            start_normal[lanes] * free_velocity[lanes, :2], axis=1
        )
        contact_normal[lanes] = start_normal[lanes]
        # The complete centre chord has already been certified against every
        # solid wall and every directed outlet.  For continuous-v16 geometry,
        # an explicit reverse-inlet guard therefore proves that the endpoint
        # remains in the same closed lumen.  This removes one GEOS/Shapely
        # point-in-polygon dispatch from every ordinary internal particle step.
        # Legacy/test geometries retain the reference membership query.
        inlet_query = (
            getattr(
                boundary_geometry,
                "directed_inlet_crossing_mask_grid_accelerated",
                None,
            )
            if use_numba
            else None
        )
        if callable(inlet_query):
            reverse_inlet = np.asarray(
                inlet_query(positions[lanes], accepted_position[lanes]),
                dtype=bool,
            ).reshape(-1)
            if reverse_inlet.shape != (lanes.size,):
                raise RuntimeError(
                    "The accelerated inlet-crossing guard returned an invalid shape."
                )
            inside = ~reverse_inlet
        else:
            inside = np.asarray(
                boundary_geometry.contains_grid(
                    accepted_position[lanes], tolerance_um=0.0
                ),
                dtype=bool,
            )
        if np.all(inside):
            return PredictiveContactStep(
                accepted_position_grid=accepted_position,
                angle_increment_rad=angle_increment,
                constrained_generalized_velocity=constrained_velocity,
                reaction_force_pn=reaction,
                endpoint_gap_um=endpoint_gap,
                minimum_path_gap_um=minimum_path_gap,
                predicted_normal_gap_um=predicted_gap,
                complementarity_residual_pn_um=complementarity,
                residual_correction_um=residual_correction,
                contact_normal_xz=contact_normal,
                active_contact=active_contact,
                outlet_fraction=outlet_fraction,
                outlet_position_grid=outlet_position,
                outlet_label=outlet_label,
                need_time_refinement=need_refinement,
                failure_codes=failure_codes,
                failure_position_grid=failure_position,
                failure_event_fraction=failure_fraction,
            )

    completed = np.zeros(count, dtype=bool)
    outlet_fraction_all = np.full(count, np.nan, dtype=np.float64)
    outlet_position_all = np.full((count, 2), np.nan, dtype=np.float64)
    outlet_label_all = np.full(count, -1, dtype=np.int32)
    outlet_fraction_all[live_lanes] = free_outlet_fraction
    outlet_position_all[live_lanes] = free_outlet_position
    outlet_label_all[live_lanes] = free_outlet_label
    wall_fraction_all = np.where(
        np.isfinite(free_path_fraction), free_path_fraction, np.inf
    )
    finite_outlet = np.isfinite(outlet_fraction_all)
    outlet_before_wall = (
        outlet_fraction_all
        <= wall_fraction_all + fraction_tolerance
    )
    # A finite-radius particle may touch a solid wall before its centre reaches
    # the adjoining outlet plane and then slide out without ever penetrating
    # that wall.  At the cap corner, the wall-contact and outlet fractions can
    # differ by roundoff even though the complete swept-disc chord is certified
    # nonpenetrating.  Treat that chord as a legal outlet prefix as well.
    #
    # The swept-path condition is essential: a centre-line outlet crossing is
    # not allowed to hide a genuine wall penetration that occurred first.
    nonpenetrating_outlet_chord = (
        finite_outlet
        & (free_path_gap >= -strict_epsilon_all)
    )
    direct_outlet = (
        live
        & finite_outlet
        & (outlet_before_wall | nonpenetrating_outlet_chord)
    )
    if np.any(direct_outlet):
        lanes = np.flatnonzero(direct_outlet)
        fractions = np.clip(outlet_fraction_all[lanes], 0.0, 1.0)
        accepted_position[lanes] = outlet_position_all[lanes]
        outlet_fraction[lanes] = fractions
        outlet_position[lanes] = outlet_position_all[lanes]
        outlet_label[lanes] = outlet_label_all[lanes]
        constrained_velocity[lanes] = free_velocity[lanes]
        angle_increment[lanes] = fractions * dt * free_velocity[lanes, 2]
        contact_normal[lanes] = start_normal[lanes]
        predicted_gap[lanes] = start_gap_all[lanes] + fractions * dt * np.sum(
            start_normal[lanes] * free_velocity[lanes, :2], axis=1
        )
        outlet_world = np.asarray(
            boundary_geometry.grid_to_world_xz(outlet_position_all[lanes]),
            dtype=np.float64,
        )
        outlet_state = (
            boundary_geometry.exact_solid_wall_state_xz_um_accelerated(
                outlet_world
            )
            if use_numba
            else boundary_geometry.exact_solid_wall_state_xz_um(outlet_world)
        )
        outlet_distance = np.asarray(
            outlet_state.distance_um, dtype=np.float64
        ).reshape(-1)
        endpoint_gap[lanes] = outlet_distance - radii[lanes]
        # Audit only the physically consumed prefix; the free proposal beyond
        # an outlet is not part of the accepted trajectory.
        prefix_gap = np.empty(lanes.size, dtype=np.float64)
        if use_numba:
            prefix_gap[:] = (
                boundary_geometry.inspect_swept_solid_wall_paths_with_end_distance_xz_um_precomputed(
                    world_start[lanes],
                    outlet_world,
                    radii[lanes],
                    start_distance[lanes],
                    tolerance_um=0.0,
                )[1]
            )
        else:
            for local, lane in enumerate(lanes):
                prefix_gap[local] = (
                    boundary_geometry.inspect_swept_solid_wall_path_xz_um(
                        world_start[lane],
                        outlet_world[local],
                        float(radii[lane]),
                        tolerance_um=float(strict_epsilon_all[lane]),
                    ).minimum_gap_um
                )
        minimum_path_gap[lanes] = np.maximum(prefix_gap, 0.0)
        completed[lanes] = True

    safe_free = (
        live
        & ~completed
        & (free_path_gap >= -strict_epsilon_all)
        & (free_endpoint_gap_all >= 0.0)
    )
    if np.any(safe_free):
        lanes = np.flatnonzero(safe_free)
        accepted_position[lanes] = boundary_geometry.world_xz_to_grid(
            free_world_end[lanes]
        )
        angle_increment[lanes] = dt * free_velocity[lanes, 2]
        constrained_velocity[lanes] = free_velocity[lanes]
        endpoint_gap[lanes] = free_endpoint_gap_all[lanes]
        minimum_path_gap[lanes] = np.maximum(free_path_gap[lanes], 0.0)
        predicted_gap[lanes] = start_gap_all[lanes] + dt * np.sum(
            start_normal[lanes] * free_velocity[lanes, :2], axis=1
        )
        contact_normal[lanes] = start_normal[lanes]
        completed[lanes] = True

    for compact_lane, lane_value in enumerate(live_lanes):
        lane = int(lane_value)
        if completed[lane]:
            continue
        radius = float(radii[lane])
        start = world_start[lane]
        free_end = free_world_end[lane]
        start_gap = float(start_distance[lane]) - radius
        strict_epsilon = _strict_geometry_epsilon(
            start, free_end, radius, spacing
        )
        if start_gap < -strict_epsilon:
            failure_codes[lane] = _STATUS_START_PENETRATION
            raise PredictiveWallContactGeometryError(
                "A v15 trial started from an already penetrating accepted state;",
                (
                    _wall_contact_failure(
                        local_lane=lane,
                        event_fraction=0.0,
                        position_xz_um=start,
                        bubble_radius_um=radius,
                        boundary_geometry=boundary_geometry,
                        geometry_tolerance_um=tolerance,
                        reason="accepted_start_has_negative_exact_gap",
                        known_gap_um=start_gap,
                    ),
                ),
            )
        crossing_fraction = float(free_outlet_fraction[compact_lane])
        wall_fraction = (
            float(free_path_fraction[lane])
            if np.isfinite(free_path_fraction[lane])
            else float("inf")
        )
        if (
            np.isfinite(crossing_fraction)
            and crossing_fraction <= wall_fraction + fraction_tolerance
        ):
            fraction = float(np.clip(crossing_fraction, 0.0, 1.0))
            accepted_position[lane] = free_outlet_position[compact_lane]
            outlet_fraction[lane] = fraction
            outlet_position[lane] = free_outlet_position[compact_lane]
            outlet_label[lane] = int(free_outlet_label[compact_lane])
            constrained_velocity[lane] = free_velocity[lane]
            angle_increment[lane] = fraction * dt * free_velocity[lane, 2]
            endpoint_gap[lane] = float(
                boundary_geometry.exact_true_gap_at_xz_um(
                    boundary_geometry.grid_to_world_xz(
                        free_outlet_position[compact_lane]
                    ),
                    radius,
                )
            )
            minimum_path_gap[lane] = float(free_path_gap[lane])
            predicted_gap[lane] = start_gap + fraction * dt * float(
                start_normal[lane] @ free_velocity[lane, :2]
            )
            contact_normal[lane] = start_normal[lane]
            continue

        free_endpoint_gap = float(free_endpoint_gap_all[lane])
        if (
            free_path_gap[lane] >= -strict_epsilon
            and free_endpoint_gap >= -strict_epsilon
        ):
            candidate_end = free_end.copy()
            correction = 0.0
            if free_endpoint_gap < 0.0:
                candidate_end, correction, free_endpoint_gap = _project_small_residual(
                    candidate_end,
                    radius,
                    tolerance,
                    strict_epsilon,
                    boundary_geometry,
                )
            accepted = (
                (
                    candidate_end,
                    max(free_endpoint_gap, 0.0),
                    max(float(free_path_gap[lane]), 0.0),
                    0.0,
                )
                if correction == 0.0
                else _certify_candidate_path(
                    start,
                    candidate_end,
                    radius,
                    tolerance,
                    strict_epsilon,
                    boundary_geometry,
                )
            )
            if accepted is not None:
                certified_end, certified_gap, certified_path, extra = accepted
                _store_accepted_lane(
                    lane,
                    certified_end,
                    free_velocity[lane],
                    0.0,
                    start_gap + dt * float(
                        start_normal[lane] @ free_velocity[lane, :2]
                    ),
                    certified_gap,
                    certified_path,
                    correction + extra,
                    start_normal[lane],
                    positions,
                    dt,
                    spacing,
                    boundary_geometry,
                    accepted_position,
                    angle_increment,
                    constrained_velocity,
                    reaction,
                    endpoint_gap,
                    minimum_path_gap,
                    predicted_gap,
                    complementarity,
                    residual_correction,
                    contact_normal,
                    active_contact,
                )
                _store_constrained_outlet_if_any(
                    lane,
                    positions[lane],
                    accepted_position[lane],
                    dt,
                    angle_increment,
                    outlet_fraction,
                    outlet_position,
                    outlet_label,
                    boundary_geometry,
                    bool(use_numba),
                )
                continue

        if free_path_multiple[lane]:
            event_position = start + (
                0.0 if not np.isfinite(wall_fraction) else wall_fraction
            ) * (free_end - start)
            raise PredictiveWallContactGeometryError(
                "A particle simultaneously contacts distinct true solid walls, "
                "which is outside the v15 single-wall mobility model;",
                (
                    _wall_contact_failure(
                        local_lane=lane,
                        event_fraction=wall_fraction if np.isfinite(wall_fraction) else 0.0,
                        position_xz_um=event_position,
                        bubble_radius_um=radius,
                        boundary_geometry=boundary_geometry,
                        geometry_tolerance_um=tolerance,
                        reason="simultaneous_distinct_solid_wall_contact",
                        known_multiple_wall_contact=True,
                    ),
                ),
            )

        contact_position = start + (
            0.0 if not np.isfinite(wall_fraction) else wall_fraction
        ) * (free_end - start)
        contact_state = boundary_geometry.exact_solid_wall_state_xz_um(
            contact_position
        )
        normal = np.asarray(contact_state.inward_normal_xz, dtype=np.float64).reshape(2)
        normal_length = float(np.linalg.norm(normal))
        if not np.isfinite(normal_length) or normal_length <= 0.0:
            raise PredictiveWallContactGeometryError(
                "A v15 rigid-wall contact does not have a unique finite wall normal;",
                (
                    _wall_contact_failure(
                        local_lane=lane,
                        event_fraction=wall_fraction if np.isfinite(wall_fraction) else 0.0,
                        position_xz_um=contact_position,
                        bubble_radius_um=radius,
                        boundary_geometry=boundary_geometry,
                        geometry_tolerance_um=tolerance,
                        reason="non_finite_or_zero_exact_contact_normal",
                    ),
                ),
            )
        normal /= normal_length
        generalized_normal = np.asarray([normal[0], normal[1], 0.0])
        denominator = float(
            generalized_normal @ mobility[lane] @ generalized_normal
        )
        denominator_epsilon = (
            512.0
            * np.finfo(np.float64).eps
            * max(float(np.linalg.norm(mobility[lane], ord=2)), 1.0)
        )
        if not np.isfinite(denominator) or denominator <= denominator_epsilon:
            need_refinement[lane] = True
            failure_codes[lane] = _STATUS_INVALID_MOBILITY
            failure_position[lane] = boundary_geometry.world_xz_to_grid(
                contact_position
            )
            failure_fraction[lane] = (
                wall_fraction if np.isfinite(wall_fraction) else 0.0
            )
            continue

        free_normal_velocity = float(generalized_normal @ free_velocity[lane])
        predicted = start_gap + dt * free_normal_velocity
        multiplier = max(
            0.0,
            -(start_gap / dt + free_normal_velocity) / denominator,
        )
        corrected_velocity = free_velocity[lane] + (
            mobility[lane] @ generalized_normal
        ) * multiplier
        candidate_end = start + dt * corrected_velocity[:2]
        candidate_gap = float(
            boundary_geometry.exact_true_gap_at_xz_um(candidate_end, radius)
        )
        correction = 0.0
        if candidate_gap < 0.0 and candidate_gap >= -tolerance:
            candidate_end, correction, candidate_gap = _project_small_residual(
                candidate_end,
                radius,
                tolerance,
                strict_epsilon,
                boundary_geometry,
            )
        if candidate_gap < 0.0:
            need_refinement[lane] = True
            failure_codes[lane] = _STATUS_ENDPOINT_PENETRATION
            failure_position[lane] = boundary_geometry.world_xz_to_grid(
                candidate_end
            )
            failure_fraction[lane] = 1.0
            continue
        certified = _certify_candidate_path(
            start,
            candidate_end,
            radius,
            tolerance,
            strict_epsilon,
            boundary_geometry,
        )
        if certified is None:
            need_refinement[lane] = True
            failure_codes[lane] = _STATUS_SWEPT_PATH_PENETRATION
            failure_position[lane] = boundary_geometry.world_xz_to_grid(
                candidate_end
            )
            failure_fraction[lane] = 1.0
            continue
        certified_end, certified_gap, certified_path, extra = certified
        complementarity_value = abs(
            multiplier
            * (start_gap + dt * float(generalized_normal @ corrected_velocity))
        )
        _store_accepted_lane(
            lane,
            certified_end,
            corrected_velocity,
            multiplier,
            predicted,
            certified_gap,
            certified_path,
            correction + extra,
            normal,
            positions,
            dt,
            spacing,
            boundary_geometry,
            accepted_position,
            angle_increment,
            constrained_velocity,
            reaction,
            endpoint_gap,
            minimum_path_gap,
            predicted_gap,
            complementarity,
            residual_correction,
            contact_normal,
            active_contact,
            complementarity_value=complementarity_value,
        )
        _store_constrained_outlet_if_any(
            lane,
            positions[lane],
            accepted_position[lane],
            dt,
            angle_increment,
            outlet_fraction,
            outlet_position,
            outlet_label,
            boundary_geometry,
            bool(use_numba),
        )

    continuing = live & ~need_refinement & ~np.isfinite(outlet_fraction)
    if np.any(continuing):
        inside = np.asarray(
            boundary_geometry.contains_grid(
                accepted_position[continuing], tolerance_um=0.0
            ),
            dtype=bool,
        )
        if not np.all(inside):
            bad = np.flatnonzero(continuing)[~inside]
            raise PredictiveBoundaryLifecycleError(
                _inspect_unclaimed_boundary_exits(
                    bad,
                    positions,
                    accepted_position,
                    world_start,
                    radii,
                    boundary_geometry,
                )
            )

    return PredictiveContactStep(
        accepted_position_grid=accepted_position,
        angle_increment_rad=angle_increment,
        constrained_generalized_velocity=constrained_velocity,
        reaction_force_pn=reaction,
        endpoint_gap_um=endpoint_gap,
        minimum_path_gap_um=minimum_path_gap,
        predicted_normal_gap_um=predicted_gap,
        complementarity_residual_pn_um=complementarity,
        residual_correction_um=residual_correction,
        contact_normal_xz=contact_normal,
        active_contact=active_contact,
        outlet_fraction=outlet_fraction,
        outlet_position_grid=outlet_position,
        outlet_label=outlet_label,
        need_time_refinement=need_refinement,
        failure_codes=failure_codes,
        failure_position_grid=failure_position,
        failure_event_fraction=failure_fraction,
    )


def _inspect_unclaimed_boundary_exits(
    lanes: np.ndarray,
    start_positions_grid: np.ndarray,
    end_positions_grid: np.ndarray,
    start_positions_xz_um: np.ndarray,
    radii_um: np.ndarray,
    boundary_geometry: ContinuousVesselGeometry,
) -> tuple[PredictiveBoundaryExitFailure, ...]:
    """Describe the nearest anatomical opening without changing the trial."""

    section_points = np.asarray(
        boundary_geometry.open_section_point_xz_um, dtype=np.float64
    )
    section_normals = np.asarray(
        boundary_geometry.open_section_outward_normal_xz, dtype=np.float64
    )
    section_tangents = np.asarray(
        boundary_geometry.open_section_tangent_xz, dtype=np.float64
    )
    section_widths = np.asarray(
        boundary_geometry.open_section_half_width_um, dtype=np.float64
    )
    section_kinds = np.asarray(boundary_geometry.open_section_kind, dtype=np.int8)
    section_labels = np.asarray(
        boundary_geometry.open_section_label, dtype=np.int32
    )
    failures: list[PredictiveBoundaryExitFailure] = []
    for raw_lane in np.asarray(lanes, dtype=np.int64):
        lane = int(raw_lane)
        start = np.asarray(start_positions_xz_um[lane], dtype=np.float64)
        end = np.asarray(
            boundary_geometry.grid_to_world_xz(end_positions_grid[lane]),
            dtype=np.float64,
        )
        relative_end = end[None, :] - section_points
        signed_start = np.einsum(
            "ij,ij->i", section_points - start[None, :], section_normals
        )
        signed_end = np.einsum(
            "ij,ij->i", section_points - end[None, :], section_normals
        )
        lateral_end = np.abs(
            np.einsum("ij,ij->i", relative_end, section_tangents)
        )
        segment_distance = np.hypot(
            np.abs(signed_end), np.maximum(lateral_end - section_widths, 0.0)
        )
        section = int(np.argmin(segment_distance))
        failures.append(
            PredictiveBoundaryExitFailure(
                local_lane=lane,
                start_position_grid=tuple(
                    float(value) for value in start_positions_grid[lane]
                ),
                end_position_grid=tuple(
                    float(value) for value in end_positions_grid[lane]
                ),
                start_position_xz_um=tuple(float(value) for value in start),
                end_position_xz_um=tuple(float(value) for value in end),
                bubble_radius_um=float(radii_um[lane]),
                nearest_open_section_index=section,
                nearest_open_section_kind=int(section_kinds[section]),
                nearest_open_section_label=int(section_labels[section]),
                section_signed_start_um=float(signed_start[section]),
                section_signed_end_um=float(signed_end[section]),
                section_lateral_end_um=float(lateral_end[section]),
                section_half_width_um=float(section_widths[section]),
            )
        )
    return tuple(failures)


def _project_small_residual(
    position_xz_um: np.ndarray,
    radius_um: float,
    tolerance_um: float,
    strict_epsilon_um: float,
    boundary_geometry: ContinuousVesselGeometry,
) -> tuple[np.ndarray, float, float]:
    """Project a roundoff-scale negative endpoint to a strictly feasible state."""

    original = np.asarray(position_xz_um, dtype=np.float64).reshape(2).copy()
    corrected = original.copy()
    gap = float(boundary_geometry.exact_true_gap_at_xz_um(corrected, radius_um))
    if gap >= 0.0:
        return corrected, 0.0, gap
    if gap < -float(tolerance_um):
        return corrected, 0.0, gap

    # A single projection can expose an adjacent face at a tessellation vertex.
    # Re-querying the authoritative nearest wall makes the correction converge
    # at corners instead of leaving a small negative state for the next RHS.
    for _ in range(8):
        state = boundary_geometry.exact_solid_wall_state_xz_um(corrected)
        gradient = np.asarray(state.inward_normal_xz, dtype=np.float64).reshape(2)
        gradient_norm_squared = float(gradient @ gradient)
        if not np.isfinite(gradient_norm_squared) or gradient_norm_squared <= 0.0:
            break
        distance = -gap + float(strict_epsilon_um)
        corrected += distance * gradient / gradient_norm_squared
        gap = float(
            boundary_geometry.exact_true_gap_at_xz_um(corrected, radius_um)
        )
        if gap >= 0.0:
            break
        if float(np.linalg.norm(corrected - original)) > (
            float(tolerance_um) + float(strict_epsilon_um)
        ):
            break
    return corrected, float(np.linalg.norm(corrected - original)), gap


def _certify_candidate_path(
    start_xz_um: np.ndarray,
    end_xz_um: np.ndarray,
    radius_um: float,
    tolerance_um: float,
    strict_epsilon_um: float,
    boundary_geometry: ContinuousVesselGeometry,
) -> tuple[np.ndarray, float, float, float] | None:
    """Return a nonpenetrating endpoint/chord or reject the full trial."""

    end = np.asarray(end_xz_um, dtype=np.float64).reshape(2).copy()
    gap = float(boundary_geometry.exact_true_gap_at_xz_um(end, radius_um))
    correction = 0.0
    if gap < 0.0:
        end, correction, gap = _project_small_residual(
            end,
            radius_um,
            tolerance_um,
            strict_epsilon_um,
            boundary_geometry,
        )
    if gap < 0.0 or correction > tolerance_um + strict_epsilon_um:
        return None
    path = boundary_geometry.inspect_swept_solid_wall_path_xz_um(
        start_xz_um,
        end,
        radius_um,
        tolerance_um=strict_epsilon_um,
    )
    if path.minimum_gap_um < -strict_epsilon_um:
        return None
    if path.multiple_wall_contact:
        return None
    return end, max(gap, 0.0), max(float(path.minimum_gap_um), 0.0), correction


def _store_accepted_lane(
    lane: int,
    end_xz_um: np.ndarray,
    velocity: np.ndarray,
    multiplier: float,
    predicted_gap_value: float,
    endpoint_gap_value: float,
    path_gap_value: float,
    correction_um: float,
    normal_xz: np.ndarray,
    start_positions_grid: np.ndarray,
    dt_s: float,
    spacing_um: float,
    boundary_geometry: ContinuousVesselGeometry,
    accepted_position: np.ndarray,
    angle_increment: np.ndarray,
    constrained_velocity: np.ndarray,
    reaction: np.ndarray,
    endpoint_gap: np.ndarray,
    minimum_path_gap: np.ndarray,
    predicted_gap: np.ndarray,
    complementarity: np.ndarray,
    residual_correction: np.ndarray,
    contact_normal: np.ndarray,
    active_contact: np.ndarray,
    *,
    complementarity_value: float = 0.0,
) -> None:
    del start_positions_grid, spacing_um
    accepted_position[lane] = boundary_geometry.world_xz_to_grid(end_xz_um)
    angle_increment[lane] = float(dt_s) * float(velocity[2])
    constrained_velocity[lane] = velocity
    reaction[lane] = float(multiplier)
    endpoint_gap[lane] = float(endpoint_gap_value)
    minimum_path_gap[lane] = float(path_gap_value)
    predicted_gap[lane] = float(predicted_gap_value)
    complementarity[lane] = float(complementarity_value)
    residual_correction[lane] = float(correction_um)
    contact_normal[lane] = normal_xz
    active_contact[lane] = bool(multiplier > 0.0)


def _store_constrained_outlet_if_any(
    lane: int,
    start_grid: np.ndarray,
    end_grid: np.ndarray,
    dt_s: float,
    angle_increment: np.ndarray,
    outlet_fraction: np.ndarray,
    outlet_position: np.ndarray,
    outlet_label: np.ndarray,
    boundary_geometry: ContinuousVesselGeometry,
    use_numba: bool,
) -> None:
    fractions, positions, labels = _first_outlet_crossing_data_grid(
        boundary_geometry,
        np.asarray(start_grid, dtype=np.float64).reshape(1, 2),
        np.asarray(end_grid, dtype=np.float64).reshape(1, 2),
        use_numba=use_numba,
    )
    if not np.isfinite(fractions[0]):
        return
    fraction = float(np.clip(fractions[0], 0.0, 1.0))
    outlet_fraction[lane] = fraction
    outlet_position[lane] = positions[0]
    outlet_label[lane] = int(labels[0])
    angle_increment[lane] *= fraction


def _first_outlet_crossing_data_grid(
    boundary_geometry: ContinuousVesselGeometry,
    start_grid: np.ndarray,
    end_grid: np.ndarray,
    *,
    use_numba: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    starts = np.asarray(start_grid, dtype=np.float64)
    ends = np.asarray(end_grid, dtype=np.float64)
    if starts.shape != ends.shape or starts.ndim != 2 or starts.shape[1:] != (2,):
        raise ValueError("Outlet path batches must have matching shape (N, 2).")
    count = int(starts.shape[0])
    if use_numba:
        fractions, indices, positions, _, _ = (
            boundary_geometry.first_outlet_crossing_arrays_grid_accelerated(
                starts, ends
            )
        )
        fractions = np.asarray(fractions, dtype=np.float64)
        indices = np.asarray(indices, dtype=np.int32)
        positions = np.asarray(positions, dtype=np.float64)
        labels = np.full(count, -1, dtype=np.int32)
        valid = indices >= 0
        if np.any(valid):
            labels[valid] = np.asarray(
                boundary_geometry.outlet_section_label, dtype=np.int32
            )[indices[valid]]
        return fractions, positions, labels
    crossings = tuple(
        boundary_geometry.first_outlet_crossings_grid(starts, ends)
    )
    fractions = np.full(count, np.nan, dtype=np.float64)
    positions = np.full((count, 2), np.nan, dtype=np.float64)
    labels = np.full(count, -1, dtype=np.int32)
    for lane, crossing in enumerate(crossings):
        if crossing is None:
            continue
        fractions[lane] = float(crossing.fraction)
        positions[lane] = boundary_geometry.world_xz_to_grid(
            crossing.position_xz_um
        )
        labels[lane] = int(crossing.label)
    return fractions, positions, labels


def _wall_contact_failure(
    *,
    local_lane: int,
    event_fraction: float,
    position_xz_um: np.ndarray,
    bubble_radius_um: float,
    boundary_geometry: ContinuousVesselGeometry,
    geometry_tolerance_um: float,
    reason: str,
    known_gap_um: float | None = None,
    known_multiple_wall_contact: bool | None = None,
) -> PredictiveWallContactFailure:
    world = np.asarray(position_xz_um, dtype=np.float64).reshape(2)
    grid = np.asarray(boundary_geometry.world_xz_to_grid(world), dtype=np.float64)
    gap = (
        float(boundary_geometry.exact_true_gap_at_xz_um(world, bubble_radius_um))
        if known_gap_um is None
        else float(known_gap_um)
    )
    walls = _candidate_walls_at_position(
        boundary_geometry, world, bubble_radius_um
    )
    if known_multiple_wall_contact is None:
        multiple = bool(
            np.asarray(
                boundary_geometry.multiple_wall_contact_mask_xz_um(
                    world,
                    bubble_radius_um,
                    tolerance_um=geometry_tolerance_um,
                )
            )
        )
    else:
        multiple = bool(known_multiple_wall_contact)
    return PredictiveWallContactFailure(
        local_lane=int(local_lane),
        event_fraction=float(np.clip(event_fraction, 0.0, 1.0)),
        position_grid=(float(grid[0]), float(grid[1])),
        position_xz_um=(float(world[0]), float(world[1])),
        gap_um=gap,
        bubble_radius_um=float(bubble_radius_um),
        candidate_walls=walls,
        candidate_wall_source="nearest_exact_distance_ties",
        multiple_wall_contact=multiple,
        reason=str(reason),
    )


def inspect_predictive_wall_contact_failure(
    *,
    local_lane: int,
    event_fraction: float,
    position_xz_um: np.ndarray,
    bubble_radius_um: float,
    boundary_geometry: ContinuousVesselGeometry,
    reason: str,
) -> PredictiveWallContactFailure:
    """Build the same exact-face diagnostic for a rejected v15 refinement."""

    world = np.asarray(position_xz_um, dtype=np.float64).reshape(2)
    scale = max(
        float(np.max(np.abs(world))), abs(float(bubble_radius_um)), 1.0
    )
    strict_tolerance = 128.0 * np.finfo(np.float64).eps * scale
    return _wall_contact_failure(
        local_lane=local_lane,
        event_fraction=event_fraction,
        position_xz_um=world,
        bubble_radius_um=bubble_radius_um,
        boundary_geometry=boundary_geometry,
        geometry_tolerance_um=strict_tolerance,
        reason=reason,
    )


def _candidate_walls_at_position(
    boundary_geometry: ContinuousVesselGeometry,
    position_xz_um: np.ndarray,
    radius_um: float,
) -> tuple[PredictiveCandidateWall, ...]:
    projected, primary, tied = boundary_geometry.nearest_solid_face_projection_xz_um(
        np.asarray(position_xz_um, dtype=np.float64).reshape(1, 2)
    )
    primary_index = int(np.asarray(primary).reshape(-1)[0])
    indices = (
        np.asarray(tied[0], dtype=np.int64).reshape(-1)
        if tied
        else np.asarray([primary_index], dtype=np.int64)
    )
    return tuple(
        _candidate_wall_record(
            boundary_geometry, int(index), position_xz_um, radius_um
        )
        for index in np.unique(indices)
    )


def _candidate_wall_record(
    geometry: ContinuousVesselGeometry,
    index: int,
    position_xz_um: np.ndarray,
    radius_um: float,
) -> PredictiveCandidateWall:
    center = np.asarray(geometry.solid_face_center_xz_um[index], dtype=np.float64)
    length = float(geometry.solid_face_length_um[index])
    point = np.asarray(position_xz_um, dtype=np.float64)
    start = np.asarray(geometry.solid_face_start_xz_um[index], dtype=np.float64)
    end = np.asarray(geometry.solid_face_end_xz_um[index], dtype=np.float64)
    delta = end - start
    fraction = float(np.clip(((point - start) @ delta) / (delta @ delta), 0.0, 1.0))
    projection = start + fraction * delta
    distance = float(np.linalg.norm(np.asarray(position_xz_um) - projection))
    inward = np.asarray(
        geometry.solid_face_inward_normal_xz[index], dtype=np.float64
    )
    return PredictiveCandidateWall(
        solid_face_index=index,
        center_xz_um=(float(center[0]), float(center[1])),
        length_um=length,
        inward_normal_xz=(float(inward[0]), float(inward[1])),
        exact_distance_um=distance,
        exact_gap_um=distance - float(radius_um),
    )


def _strict_geometry_epsilon(
    start_xz_um: np.ndarray,
    end_xz_um: np.ndarray,
    radius_um: float,
    spacing_um: float,
) -> float:
    return particle_geometry_roundoff_tolerance_um(
        start_xz_um,
        end_xz_um,
        radius_um,
        spacing_um,
    )


def _fraction_tolerance(dt_s: float) -> float:
    return 128.0 * np.finfo(np.float64).eps * max(abs(float(dt_s)), 1.0)


def _validate_inputs(
    positions: np.ndarray,
    alive: np.ndarray,
    dt_s: float,
    velocity: np.ndarray,
    mobility: np.ndarray,
    radii: np.ndarray,
    spacing_um: float,
    tolerance_um: float,
    boundary_geometry: ContinuousVesselGeometry,
) -> None:
    count = int(positions.shape[0]) if positions.ndim == 2 else -1
    if positions.shape != (count, 2) or alive.shape != (count,):
        raise ValueError("Particle positions/alive flags have inconsistent shapes.")
    if velocity.shape != (count, 3) or mobility.shape != (count, 3, 3):
        raise ValueError("Generalized velocity/mobility arrays have inconsistent shapes.")
    if radii.shape != (count,):
        raise ValueError("Bubble radii have inconsistent particle shapes.")
    if (
        not np.isfinite(dt_s)
        or dt_s <= 0.0
        or not np.isfinite(spacing_um)
        or spacing_um <= 0.0
        or not np.isfinite(tolerance_um)
        or tolerance_um <= 0.0
    ):
        raise ValueError("Time, spacing, and geometry tolerance must be positive.")
    if (
        np.any(~np.isfinite(positions))
        or np.any(~np.isfinite(velocity))
        or np.any(~np.isfinite(mobility))
        or np.any(~np.isfinite(radii))
        or np.any(radii <= 0.0)
    ):
        raise ValueError("The v15 predictive transaction received invalid inputs.")
