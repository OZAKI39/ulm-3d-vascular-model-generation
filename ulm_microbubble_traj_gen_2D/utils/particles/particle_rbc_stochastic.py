"""Zero-time continuous-geometry sweep for RBC stochastic displacement."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from .particle_geometry_tolerance import (
    particle_geometry_roundoff_tolerance_um,
)
from .particle_topological_ownership import (
    TopologicalCommitmentCatalog,
    inspect_topological_crossings,
)


_TERMINATION_OUTLET = np.uint8(1)
_MAXIMUM_EVENT_ROUNDS = 32
_EVENT_OUTLET = np.uint8(1)
_EVENT_TOPOLOGY = np.uint8(2)
_EVENT_WALL = np.uint8(3)


@dataclass(frozen=True, slots=True)
class RbcStochasticSweepResult:
    position_grid: np.ndarray
    alive: np.ndarray
    vessel_id: np.ndarray
    termination_reason: np.ndarray
    exit_time_s: np.ndarray
    reflection_count: np.ndarray
    requested_displacement_um: np.ndarray
    event_count: np.ndarray
    event_code: np.ndarray
    event_fraction: np.ndarray
    event_position_xz_um: np.ndarray
    topology_event_bubble_id: np.ndarray
    topology_event_from_vessel_id: np.ndarray
    topology_event_to_vessel_id: np.ndarray
    topology_event_section_index: np.ndarray
    topology_event_position_xz_um: np.ndarray


@dataclass(frozen=True, slots=True)
class RbcStochasticGeometryFailure:
    permanent_microbubble_id: int
    physical_time_s: float
    position_xz_um: tuple[float, float]
    requested_displacement_xz_um: tuple[float, float]
    remaining_displacement_xz_um: tuple[float, float]
    wall_gap_um: float
    inward_wall_normal_xz: tuple[float, float]
    event_codes: tuple[int, ...]
    event_fractions: tuple[float, ...]


class RbcStochasticGeometryError(RuntimeError):
    """A random geometric path could not be completed without approximation."""

    def __init__(self, message: str, failures: tuple[RbcStochasticGeometryFailure, ...]):
        self.failures = failures
        super().__init__(message)


def sweep_rbc_stochastic_displacement(
    positions_grid: np.ndarray,
    displacement_xz_um: np.ndarray,
    permanent_ids: np.ndarray,
    alive: np.ndarray,
    vessel_id: np.ndarray,
    radii_um: np.ndarray,
    termination_reason: np.ndarray,
    exit_time_s: np.ndarray,
    *,
    physical_time_s: float,
    boundary_geometry: ContinuousVesselGeometry,
    topological_ownership: TopologicalCommitmentCatalog,
    use_numba: bool,
) -> RbcStochasticSweepResult:
    """Complete one reflected random chord without consuming physical time."""

    positions = np.ascontiguousarray(positions_grid, dtype=np.float64).copy()
    displacement = np.ascontiguousarray(displacement_xz_um, dtype=np.float64)
    ids = np.ascontiguousarray(permanent_ids, dtype=np.int64)
    live = np.ascontiguousarray(alive, dtype=bool).copy()
    owners = np.ascontiguousarray(vessel_id, dtype=np.int32).copy()
    radii = np.ascontiguousarray(radii_um, dtype=np.float64)
    termination = np.ascontiguousarray(termination_reason, dtype=np.uint8).copy()
    exit_times = np.ascontiguousarray(exit_time_s, dtype=np.float64).copy()
    count = int(ids.size)
    if not (
        positions.shape == displacement.shape == (count, 2)
        and live.shape == owners.shape == radii.shape == termination.shape == exit_times.shape == (count,)
    ):
        raise ValueError("RBC stochastic sweep arrays have inconsistent shapes.")
    if np.any(~np.isfinite(positions)) or np.any(~np.isfinite(displacement)):
        raise ValueError("RBC stochastic positions and displacements must be finite.")
    if np.any(live & ((owners <= 0) | (radii <= 0.0))):
        raise ValueError("Every live stochastic lane needs a vessel owner and radius.")

    requested_norm = np.linalg.norm(displacement, axis=1)
    remaining = displacement.copy()
    remaining[~live] = 0.0
    initial_world = np.asarray(
        boundary_geometry.grid_to_world_xz(positions[live]),
        dtype=np.float64,
    )
    residual_tolerance_um = particle_geometry_roundoff_tolerance_um(
        initial_world,
        radii[live],
        displacement[live],
    )
    reflections = np.zeros(count, dtype=np.int32)
    event_count = np.zeros(count, dtype=np.int32)
    event_lanes_log: list[np.ndarray] = []
    event_code_log: list[np.uint8] = []
    event_fraction_log: list[np.ndarray] = []
    event_position_log: list[np.ndarray] = []
    topology_ids: list[np.ndarray] = []
    topology_from: list[np.ndarray] = []
    topology_to: list[np.ndarray] = []
    topology_section: list[np.ndarray] = []
    topology_position: list[np.ndarray] = []

    for event_round in range(_MAXIMUM_EVENT_ROUNDS):
        remaining_norm = np.linalg.norm(remaining, axis=1)
        roundoff_residual = live & (
            remaining_norm <= residual_tolerance_um
        )
        remaining[roundoff_residual] = 0.0
        pending = live & (remaining_norm > residual_tolerance_um)
        if not np.any(pending):
            break

        lanes = np.flatnonzero(pending)
        starts_grid = positions[lanes]
        starts_world = np.asarray(
            boundary_geometry.grid_to_world_xz(starts_grid), dtype=np.float64
        )
        ends_world = starts_world + remaining[lanes]
        ends_grid = np.asarray(
            boundary_geometry.world_xz_to_grid(ends_world), dtype=np.float64
        )

        inlet_fraction = np.full(lanes.size, np.inf, dtype=np.float64)
        inlet_values = (
            boundary_geometry.first_inlet_crossing_arrays_grid_accelerated(
                starts_grid, ends_grid
            )
        )
        inlet_fraction_raw = np.asarray(
            inlet_values[0], dtype=np.float64
        ).reshape(-1)
        inlet_fraction = np.where(
            np.isfinite(inlet_fraction_raw), inlet_fraction_raw, np.inf
        )

        outlet_values = boundary_geometry.first_outlet_crossing_arrays_grid_accelerated(
            starts_grid, ends_grid
        )
        outlet_fraction = np.asarray(outlet_values[0], dtype=np.float64).reshape(-1)
        outlet_position_grid = np.asarray(outlet_values[2], dtype=np.float64).reshape(-1, 2)
        outlet_fraction = np.where(np.isfinite(outlet_fraction), outlet_fraction, np.inf)

        topology = inspect_topological_crossings(
            starts_world,
            ends_world,
            owners[lanes],
            np.ones(lanes.size, dtype=bool),
            topological_ownership,
            use_numba=bool(use_numba),
        )
        topology_fraction = np.asarray(topology.fraction, dtype=np.float64)
        topology_fraction = np.where(
            np.isfinite(topology_fraction), topology_fraction, np.inf
        )

        wall_values = boundary_geometry.inspect_swept_solid_wall_paths_xz_um_accelerated(
            starts_world,
            ends_world,
            radii[lanes],
            tolerance_um=0.0,
        )
        wall_fraction = np.asarray(wall_values[1], dtype=np.float64).reshape(-1)
        wall_face = np.asarray(wall_values[2], dtype=np.int64).reshape(-1)
        multiple_wall = np.asarray(wall_values[3], dtype=bool).reshape(-1)
        wall_fraction = np.where(np.isfinite(wall_fraction), wall_fraction, np.inf)

        scale = np.maximum.reduce(
            [
                np.where(np.isfinite(outlet_fraction), np.abs(outlet_fraction), 0.0),
                np.where(
                    np.isfinite(topology_fraction),
                    np.abs(topology_fraction),
                    0.0,
                ),
                np.where(np.isfinite(wall_fraction), np.abs(wall_fraction), 0.0),
                np.where(np.isfinite(inlet_fraction), np.abs(inlet_fraction), 0.0),
                np.ones(lanes.size, dtype=np.float64),
            ]
        )
        tie_tolerance = 256.0 * np.finfo(np.float64).eps * scale
        first_valid_fraction = np.minimum.reduce(
            [outlet_fraction, topology_fraction, wall_fraction]
        )
        inlet_event = np.isfinite(inlet_fraction) & (
            inlet_fraction <= first_valid_fraction + tie_tolerance
        )
        if np.any(inlet_event):
            raise _geometry_error(
                "An RBC random displacement attempted to leave through an inlet.",
                lanes[inlet_event],
                ids,
                positions,
                displacement,
                remaining,
                radii,
                float(physical_time_s),
                boundary_geometry,
                event_lanes_log,
                event_code_log,
                event_fraction_log,
                event_position_log,
                event_count,
            )
        outlet_event = np.isfinite(outlet_fraction) & (
            outlet_fraction <= topology_fraction + tie_tolerance
        ) & (outlet_fraction <= wall_fraction + tie_tolerance)
        topology_event = (~outlet_event) & np.isfinite(topology_fraction) & (
            topology_fraction <= wall_fraction + tie_tolerance
        )
        wall_event = (~outlet_event) & (~topology_event) & np.isfinite(wall_fraction)
        retry_departure = np.zeros(lanes.size, dtype=bool)
        zero_wall = wall_event & (~multiple_wall) & (
            wall_fraction <= 256.0 * np.finfo(np.float64).eps
        )
        if np.any(zero_wall):
            local = np.flatnonzero(zero_wall)
            accepted = lanes[local]
            face = wall_face[local]
            inward_all = getattr(
                boundary_geometry, "solid_face_inward_normal_xz", None
            )
            inward = (
                np.asarray(inward_all[face], dtype=np.float64)
                if inward_all is not None
                else -np.asarray(
                    boundary_geometry.solid_face_outward_normal_xz[face],
                    dtype=np.float64,
                )
            )
            normal_component = np.sum(remaining[accepted] * inward, axis=1)
            direction_tolerance = (
                256.0
                * np.finfo(np.float64).eps
                * np.maximum(np.linalg.norm(remaining[accepted], axis=1), 1.0)
            )
            feasible_departure = normal_component >= -direction_tolerance
            if np.any(feasible_departure):
                depart_local = local[feasible_departure]
                depart_accepted = accepted[feasible_departure]
                depart_inward = inward[feasible_departure]
                depart_component = np.maximum(
                    normal_component[feasible_departure], 0.0
                )
                coordinate_scale = np.maximum.reduce(
                    [
                        np.max(np.abs(starts_world[depart_local]), axis=1),
                        radii[depart_accepted],
                        np.ones(depart_accepted.size, dtype=np.float64),
                    ]
                )
                nudge = 4096.0 * np.finfo(np.float64).eps * coordinate_scale
                compensated = np.minimum(nudge, 0.5 * depart_component)
                applied = np.where(depart_component > 0.0, compensated, nudge)
                departed_world = (
                    starts_world[depart_local]
                    + applied[:, None] * depart_inward
                )
                positions[depart_accepted] = boundary_geometry.world_xz_to_grid(
                    departed_world
                )
                remaining[depart_accepted] -= (
                    compensated[:, None] * depart_inward
                )
                retry_departure[depart_local] = True
                wall_event[depart_local] = False
        no_event = ~(
            outlet_event | topology_event | wall_event | retry_departure
        )

        if np.any(no_event):
            accepted = lanes[no_event]
            positions[accepted] = ends_grid[no_event]
            remaining[accepted] = 0.0

        if np.any(outlet_event):
            local = np.flatnonzero(outlet_event)
            accepted = lanes[local]
            positions[accepted] = outlet_position_grid[local]
            live[accepted] = False
            termination[accepted] = _TERMINATION_OUTLET
            exit_times[accepted] = float(physical_time_s)
            remaining[accepted] = 0.0
            _store_events(
                accepted,
                _EVENT_OUTLET,
                outlet_fraction[local],
                boundary_geometry.grid_to_world_xz(outlet_position_grid[local]),
                event_lanes_log,
                event_code_log,
                event_fraction_log,
                event_position_log,
                event_count,
            )

        if np.any(topology_event):
            local = np.flatnonzero(topology_event)
            accepted = lanes[local]
            fractions = topology_fraction[local]
            crossing_world = np.asarray(
                topology.position_xz_um[local], dtype=np.float64
            )
            previous = owners[accepted].copy()
            owners[accepted] = np.asarray(
                topology.new_vessel_id[local], dtype=np.int32
            )
            positions[accepted] = boundary_geometry.world_xz_to_grid(crossing_world)
            remaining[accepted] *= (1.0 - fractions)[:, None]
            _store_events(
                accepted,
                _EVENT_TOPOLOGY,
                fractions,
                crossing_world,
                event_lanes_log,
                event_code_log,
                event_fraction_log,
                event_position_log,
                event_count,
            )
            topology_ids.append(ids[accepted].copy())
            topology_from.append(previous)
            topology_to.append(owners[accepted].copy())
            topology_section.append(
                np.asarray(topology.section_index[local], dtype=np.int32).copy()
            )
            topology_position.append(crossing_world.copy())

        if np.any(wall_event):
            local = np.flatnonzero(wall_event)
            accepted = lanes[local]
            if np.any(multiple_wall[local]):
                failed = accepted[multiple_wall[local]]
                raise _geometry_error(
                    "An RBC random displacement reached simultaneous solid walls.",
                    failed,
                    ids,
                    positions,
                    displacement,
                    remaining,
                    radii,
                    float(physical_time_s),
                    boundary_geometry,
                    event_lanes_log,
                    event_code_log,
                    event_fraction_log,
                    event_position_log,
                    event_count,
                )
            fractions = wall_fraction[local]
            crossing_world = starts_world[local] + fractions[:, None] * (
                ends_world[local] - starts_world[local]
            )
            face = wall_face[local]
            inward_all = getattr(
                boundary_geometry, "solid_face_inward_normal_xz", None
            )
            inward = (
                np.asarray(inward_all[face], dtype=np.float64)
                if inward_all is not None
                else -np.asarray(
                    boundary_geometry.solid_face_outward_normal_xz[face],
                    dtype=np.float64,
                )
            )
            remainder = remaining[accepted] * (1.0 - fractions)[:, None]
            projected = np.sum(remainder * inward, axis=1)
            reflected = remainder - 2.0 * projected[:, None] * inward

            coordinate_scale = np.maximum.reduce(
                [
                    np.max(np.abs(crossing_world), axis=1),
                    radii[accepted],
                    np.ones(accepted.size, dtype=np.float64),
                ]
            )
            nudge = 4096.0 * np.finfo(np.float64).eps * coordinate_scale
            available_inward = np.maximum(np.sum(reflected * inward, axis=1), 0.0)
            compensated = np.minimum(nudge, 0.5 * available_inward)
            # A tangential remainder has no inward component to compensate.  One
            # numerical-scale inward offset prevents a zero-fraction re-hit; the
            # final exact-gap audit makes this explicit and fail-closed.
            has_remainder = (
                np.linalg.norm(reflected, axis=1) > residual_tolerance_um
            )
            applied_nudge = np.where(
                has_remainder,
                np.where(available_inward > 0.0, compensated, nudge),
                0.0,
            )
            next_world = crossing_world + applied_nudge[:, None] * inward
            reflected -= compensated[:, None] * inward
            positions[accepted] = boundary_geometry.world_xz_to_grid(next_world)
            remaining[accepted] = reflected
            reflections[accepted] += 1
            _store_events(
                accepted,
                _EVENT_WALL,
                fractions,
                crossing_world,
                event_lanes_log,
                event_code_log,
                event_fraction_log,
                event_position_log,
                event_count,
            )
    else:
        pending = live & (
            np.linalg.norm(remaining, axis=1) > residual_tolerance_um
        )
        if np.any(pending):
            raise _geometry_error(
                "An RBC random displacement exceeded the geometric event limit.",
                np.flatnonzero(pending),
                ids,
                positions,
                displacement,
                remaining,
                radii,
                float(physical_time_s),
                boundary_geometry,
                event_lanes_log,
                event_code_log,
                event_fraction_log,
                event_position_log,
                event_count,
            )

    final_live = np.flatnonzero(live)
    if final_live.size:
        final_world = np.asarray(
            boundary_geometry.grid_to_world_xz(positions[final_live]),
            dtype=np.float64,
        )
        strict_roundoff = particle_geometry_roundoff_tolerance_um(
            final_world, radii[final_live]
        )
        gaps = _project_roundoff_negative_final_gaps(
            positions,
            final_live,
            radii,
            boundary_geometry,
            strict_roundoff,
        )
        failed_mask = (~np.isfinite(gaps)) | (gaps < -strict_roundoff)
        # A tolerance can authorize a projection, never a stored negative gap.
        failed_mask |= gaps < 0.0
        if np.any(failed_mask):
            raise _geometry_error(
                "An RBC random displacement ended outside the feasible lumen.",
                final_live[failed_mask],
                ids,
                positions,
                displacement,
                remaining,
                radii,
                float(physical_time_s),
                boundary_geometry,
                event_lanes_log,
                event_code_log,
                event_fraction_log,
                event_position_log,
                event_count,
            )

    event_code, event_fraction, event_position = _pack_event_log(
        count,
        event_count,
        event_lanes_log,
        event_code_log,
        event_fraction_log,
        event_position_log,
    )

    return RbcStochasticSweepResult(
        position_grid=positions,
        alive=live,
        vessel_id=owners,
        termination_reason=termination,
        exit_time_s=exit_times,
        reflection_count=reflections,
        requested_displacement_um=requested_norm,
        event_count=event_count,
        event_code=event_code,
        event_fraction=event_fraction,
        event_position_xz_um=event_position,
        topology_event_bubble_id=_concatenate(topology_ids, np.int64),
        topology_event_from_vessel_id=_concatenate(topology_from, np.int32),
        topology_event_to_vessel_id=_concatenate(topology_to, np.int32),
        topology_event_section_index=_concatenate(topology_section, np.int32),
        topology_event_position_xz_um=(
            np.vstack(topology_position).astype(np.float64, copy=False)
            if topology_position
            else np.empty((0, 2), dtype=np.float64)
        ),
    )


def _project_roundoff_negative_final_gaps(
    positions_grid: np.ndarray,
    live_lanes: np.ndarray,
    radii_um: np.ndarray,
    boundary_geometry: ContinuousVesselGeometry,
    tolerance_um: float,
) -> np.ndarray:
    """Move roundoff-negative stochastic endpoints strictly into the lumen."""

    lanes = np.asarray(live_lanes, dtype=np.int64)
    radii = np.asarray(radii_um, dtype=np.float64)
    gaps = np.empty(lanes.size, dtype=np.float64)
    for _ in range(8):
        world = np.asarray(
            boundary_geometry.grid_to_world_xz(positions_grid[lanes]),
            dtype=np.float64,
        )
        exact = boundary_geometry.exact_solid_wall_state_xz_um_accelerated(world)
        gaps = (
            np.asarray(exact.distance_um, dtype=np.float64).reshape(-1)
            - radii[lanes]
        )
        eligible = np.isfinite(gaps) & (gaps < 0.0) & (gaps >= -tolerance_um)
        if not np.any(eligible):
            return gaps
        normals = np.asarray(
            exact.inward_normal_xz, dtype=np.float64
        ).reshape(-1, 2)[eligible]
        norm_squared = np.sum(normals * normals, axis=1)
        if np.any(~np.isfinite(norm_squared)) or np.any(norm_squared <= 0.0):
            return gaps
        correction = (
            (-gaps[eligible] + tolerance_um)[:, None]
            * normals
            / norm_squared[:, None]
        )
        positions_grid[lanes[eligible]] = boundary_geometry.world_xz_to_grid(
            world[eligible] + correction
        )
    world = np.asarray(
        boundary_geometry.grid_to_world_xz(positions_grid[lanes]),
        dtype=np.float64,
    )
    exact = boundary_geometry.exact_solid_wall_state_xz_um_accelerated(world)
    return (
        np.asarray(exact.distance_um, dtype=np.float64).reshape(-1)
        - radii[lanes]
    )


def _store_events(
    lanes,
    code,
    fractions,
    positions,
    event_lanes_log,
    event_code_log,
    event_fraction_log,
    event_position_log,
    event_count,
):
    accepted = np.asarray(lanes, dtype=np.int64).copy()
    event_lanes_log.append(accepted)
    event_code_log.append(np.uint8(code))
    event_fraction_log.append(np.asarray(fractions, dtype=np.float64).copy())
    event_position_log.append(
        np.asarray(positions, dtype=np.float64).reshape(-1, 2).copy()
    )
    event_count[lanes] += 1


def _pack_event_log(
    particle_count,
    event_count,
    event_lanes_log,
    event_code_log,
    event_fraction_log,
    event_position_log,
):
    """Pack only realized events; an event-free step allocates zero rows."""

    rows = int(np.max(event_count, initial=0))
    event_code = np.zeros((rows, particle_count), dtype=np.uint8)
    event_fraction = np.full((rows, particle_count), np.nan, dtype=np.float64)
    event_position = np.full(
        (rows, particle_count, 2), np.nan, dtype=np.float64
    )
    ordinal = np.zeros(particle_count, dtype=np.int32)
    for lanes, code, fractions, positions in zip(
        event_lanes_log,
        event_code_log,
        event_fraction_log,
        event_position_log,
        strict=True,
    ):
        row = ordinal[lanes]
        event_code[row, lanes] = code
        event_fraction[row, lanes] = fractions
        event_position[row, lanes] = positions
        ordinal[lanes] += 1
    return event_code, event_fraction, event_position


def _concatenate(values: list[np.ndarray], dtype: object) -> np.ndarray:
    return (
        np.concatenate(values).astype(dtype, copy=False)
        if values
        else np.empty(0, dtype=dtype)
    )


def _geometry_error(
    message,
    failed_lanes,
    permanent_ids,
    positions_grid,
    requested_displacement,
    remaining,
    radii,
    physical_time_s,
    geometry,
    event_lanes_log,
    event_code_log,
    event_fraction_log,
    event_position_log,
    event_count,
):
    event_code, event_fraction, _ = _pack_event_log(
        int(np.asarray(permanent_ids).size),
        event_count,
        event_lanes_log,
        event_code_log,
        event_fraction_log,
        event_position_log,
    )
    failures: list[RbcStochasticGeometryFailure] = []
    for lane_value in np.asarray(failed_lanes, dtype=np.int64):
        lane = int(lane_value)
        world = np.asarray(
            geometry.grid_to_world_xz(positions_grid[lane]), dtype=np.float64
        ).reshape(2)
        exact = geometry.exact_solid_wall_state_xz_um(world)
        gap = float(np.asarray(exact.distance_um)) - float(radii[lane])
        normal = np.asarray(exact.inward_normal_xz, dtype=np.float64).reshape(2)
        stored = np.asarray(event_code[:, lane] != 0, dtype=bool)
        failures.append(
            RbcStochasticGeometryFailure(
                permanent_microbubble_id=int(permanent_ids[lane]),
                physical_time_s=float(physical_time_s),
                position_xz_um=(float(world[0]), float(world[1])),
                requested_displacement_xz_um=(
                    float(requested_displacement[lane, 0]),
                    float(requested_displacement[lane, 1]),
                ),
                remaining_displacement_xz_um=(
                    float(remaining[lane, 0]),
                    float(remaining[lane, 1]),
                ),
                wall_gap_um=gap,
                inward_wall_normal_xz=(float(normal[0]), float(normal[1])),
                event_codes=tuple(
                    int(value) for value in event_code[stored, lane]
                ),
                event_fractions=tuple(
                    float(value) for value in event_fraction[stored, lane]
                ),
            )
        )
    detail = "; ".join(
        f"id={failure.permanent_microbubble_id}, position={failure.position_xz_um}, "
        f"requested={failure.requested_displacement_xz_um}, "
        f"remaining={failure.remaining_displacement_xz_um}, gap={failure.wall_gap_um:.17g}, "
        f"normal={failure.inward_wall_normal_xz}, events={failure.event_codes}"
        for failure in failures
    )
    return RbcStochasticGeometryError(f"{message} {detail}", tuple(failures))
