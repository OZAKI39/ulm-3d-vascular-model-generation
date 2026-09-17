"""Generate topology-based candidate vessel-bed masks after an accepted CFD solve."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy import ndimage

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover - the production environment includes Numba
    njit = None
    prange = range

from ulm_vascular_model_generator.utils.core.models import Vessel

from ..particles.particle_hydrodynamic_fields import ParticleHydrodynamicFields
from .molecular_target_wall_measure import build_molecular_target_wall_measure
from ..core.types import FlowField, GridDomain, RasterizedVessels
from ..geometry.vessel_bed_topology import VesselBedTopology, VesselBedUnit, build_vessel_bed_topology


@dataclass(frozen=True)
class MolecularTargetCandidate:
    """One selectable local vessel unit or downstream vessel-bed subtree."""

    candidate_id: str
    kind: str
    label: str
    root_unit_id: int
    member_unit_ids: tuple[int, ...]
    parent_candidate_id: str | None
    depth: int
    topology_depth: int
    inlet_flow_um3_s: float
    root_flow_fraction: float
    network_flow_fraction: float
    volume_um3: float
    residence_time_s: float
    mean_wall_shear_pa: float
    endothelial_wall_area_um2: float
    endothelial_wall_area_fraction: float
    wall_area_centroid_x_um: float
    wall_area_centroid_z_um: float
    radius_of_gyration_um: float
    expected_bubble_visits: float


@dataclass(frozen=True)
class MolecularTargetCandidateCatalog:
    """Compact candidate catalog backed by one unit-ownership grid."""

    x_coordinates_um: np.ndarray
    z_coordinates_um: np.ndarray
    unit_id_grid: np.ndarray
    candidate_support_mask: np.ndarray
    lumen_mask: np.ndarray
    solid_wall_mask: np.ndarray
    open_boundary_mask: np.ndarray
    topology: VesselBedTopology
    candidates: tuple[MolecularTargetCandidate, ...]
    unresolved_junction_cells: int
    network_endothelial_wall_area_um2: float = math.nan
    network_inlet_flow_um3_s: float = math.nan
    injection_rate_per_s: float = math.nan
    observation_time_s: float = math.nan
    wall_area_weight_um2: np.ndarray | None = None
    wall_segment_id_grid: np.ndarray | None = None
    accessible_wall_mask: np.ndarray | None = None
    expected_bubble_visits_by_unit: np.ndarray | None = None
    mapped_endothelial_wall_area_um2: float = math.nan
    unmapped_endothelial_wall_area_um2: float = math.nan
    schema_version: str = "v3"

    @property
    def shape(self) -> tuple[int, int]:
        return (int(self.x_coordinates_um.size), int(self.z_coordinates_um.size))

    def candidate_by_id(self, candidate_id: str) -> MolecularTargetCandidate:
        for candidate in self.candidates:
            if candidate.candidate_id == candidate_id:
                return candidate
        raise KeyError(f"Unknown molecular-target candidate: {candidate_id}")

    @property
    def automatic_metrics_available(self) -> bool:
        return bool(
            self.schema_version == "v3"
            and math.isfinite(self.network_endothelial_wall_area_um2)
            and self.network_endothelial_wall_area_um2 > 0.0
            and math.isfinite(self.network_inlet_flow_um3_s)
            and self.network_inlet_flow_um3_s > 0.0
            and math.isfinite(self.injection_rate_per_s)
            and self.injection_rate_per_s > 0.0
            and math.isfinite(self.observation_time_s)
            and self.observation_time_s > 0.0
            and self.wall_area_weight_um2 is not None
            and self.wall_area_weight_um2.shape == self.shape
            and self.accessible_wall_mask is not None
            and self.accessible_wall_mask.shape == self.shape
            and np.any(self.accessible_wall_mask)
        )

    def mask_for_candidate_ids(self, candidate_ids: list[str] | tuple[str, ...]) -> np.ndarray:
        """Return the Boolean union of selected candidate vessel beds."""

        selected_units: set[int] = set()
        for candidate_id in simplify_candidate_selection(self, candidate_ids):
            selected_units.update(self.candidate_by_id(candidate_id).member_unit_ids)
        if not selected_units:
            return np.zeros(self.shape, dtype=bool)
        unit_values = np.asarray(sorted(selected_units), dtype=np.int32)
        mask = np.isin(self.unit_id_grid, unit_values)
        mask &= self.candidate_support_mask
        mask[self.open_boundary_mask] = False
        return np.ascontiguousarray(mask, dtype=bool)


def build_molecular_target_candidates(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    hydrodynamic_fields: ParticleHydrodynamicFields,
    vessels: list[Vessel],
    *,
    injection_rate_per_s: float = math.nan,
    observation_time_s: float = math.nan,
) -> MolecularTargetCandidateCatalog:
    """Build candidate vessel beds and their steady-flow descriptors."""

    topology = build_vessel_bed_topology(vessels)
    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    open_boundary = np.asarray(hydrodynamic_fields.open_boundary_mask, dtype=bool)
    solid_wall = np.asarray(hydrodynamic_fields.solid_site_mask, dtype=bool)
    unit_grid = _initial_unit_ownership(raster, lumen, topology)
    unit_grid, unresolved_count = _resolve_junction_ownership(
        domain,
        raster,
        flow,
        lumen,
        unit_grid,
    )
    unit_grid = _extend_ownership_to_solid_wall(unit_grid, lumen, solid_wall)
    support = ((lumen & ~open_boundary) | solid_wall) & (unit_grid >= 0)

    _validate_automatic_selection_context(injection_rate_per_s, observation_time_s)
    network_wall_area = float(
        sum(unit.endothelial_wall_area_um2 for unit in topology.units)
    )
    network_inlet_flow = float(
        sum(topology.units[root_id].flow_rate_um3_s for root_id in topology.root_unit_ids)
    )
    unit_mean_wss = _unit_mean_wall_shear(raster, flow, unit_grid, topology)
    candidates = _build_candidates(
        topology,
        unit_mean_wss,
        network_endothelial_wall_area_um2=network_wall_area,
        network_inlet_flow_um3_s=network_inlet_flow,
        injection_rate_per_s=float(injection_rate_per_s),
        observation_time_s=float(observation_time_s),
    )
    if not candidates:
        raise ValueError("The accepted vessel tree produced no perfused target candidates.")
    wall_measure = build_molecular_target_wall_measure(
        raster,
        topology,
        unit_grid,
        solid_wall,
        open_boundary,
        vessels,
        injection_rate_per_s=float(injection_rate_per_s),
        observation_time_s=float(observation_time_s),
        network_inlet_flow_um3_s=network_inlet_flow,
    )

    return MolecularTargetCandidateCatalog(
        x_coordinates_um=np.ascontiguousarray(domain.x_coordinates_um, dtype=np.float64),
        z_coordinates_um=np.ascontiguousarray(domain.z_coordinates_um, dtype=np.float64),
        unit_id_grid=np.ascontiguousarray(unit_grid, dtype=np.int32),
        candidate_support_mask=np.ascontiguousarray(support, dtype=bool),
        lumen_mask=np.ascontiguousarray(lumen, dtype=bool),
        solid_wall_mask=np.ascontiguousarray(solid_wall, dtype=bool),
        open_boundary_mask=np.ascontiguousarray(open_boundary, dtype=bool),
        topology=topology,
        candidates=candidates,
        unresolved_junction_cells=int(unresolved_count),
        network_endothelial_wall_area_um2=network_wall_area,
        network_inlet_flow_um3_s=network_inlet_flow,
        injection_rate_per_s=float(injection_rate_per_s),
        observation_time_s=float(observation_time_s),
        wall_area_weight_um2=wall_measure.wall_area_weight_um2,
        wall_segment_id_grid=wall_measure.wall_segment_id_grid,
        accessible_wall_mask=wall_measure.accessible_wall_mask,
        expected_bubble_visits_by_unit=wall_measure.expected_bubble_visits_by_unit,
        mapped_endothelial_wall_area_um2=(
            wall_measure.mapped_endothelial_wall_area_um2
        ),
        unmapped_endothelial_wall_area_um2=(
            wall_measure.unmapped_endothelial_wall_area_um2
        ),
    )


def simplify_candidate_selection(
    catalog: MolecularTargetCandidateCatalog,
    candidate_ids: list[str] | tuple[str, ...],
) -> tuple[str, ...]:
    """Remove duplicates and candidates already contained by a selected parent."""

    unique = tuple(dict.fromkeys(str(value) for value in candidate_ids))
    candidates = [catalog.candidate_by_id(candidate_id) for candidate_id in unique]
    member_sets = {
        candidate.candidate_id: frozenset(candidate.member_unit_ids)
        for candidate in candidates
    }
    retained: list[str] = []
    for candidate in candidates:
        members = member_sets[candidate.candidate_id]
        covered = any(
            candidate.candidate_id != other.candidate_id
            and members < member_sets[other.candidate_id]
            for other in candidates
        )
        if not covered:
            retained.append(candidate.candidate_id)
    return tuple(retained)


def _initial_unit_ownership(
    raster: RasterizedVessels,
    lumen: np.ndarray,
    topology: VesselBedTopology,
) -> np.ndarray:
    vessel_ids = np.asarray(raster.vessel_id, dtype=np.int64)
    ownership = np.full(lumen.shape, -1, dtype=np.int32)
    for segment_id, unit_id in topology.segment_id_to_unit_id.items():
        ownership[lumen & (vessel_ids == int(segment_id))] = int(unit_id)
    if np.any(lumen & (ownership < 0)):
        raise ValueError("Some lumen cells have no vessel-unit ownership.")
    return ownership


def _resolve_junction_ownership(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    lumen: np.ndarray,
    ownership: np.ndarray,
    *,
    use_numba: bool = True,
) -> tuple[np.ndarray, int]:
    junction = raster.junction_core_mask
    if junction is None:
        return ownership, 0
    junction_cells = np.asarray(junction, dtype=bool) & lumen
    if not np.any(junction_cells):
        return ownership, 0

    resolved = ownership.copy()
    resolved[junction_cells] = -1
    velocity = np.asarray(flow.velocity_xz_um_s, dtype=np.float64)
    step_um = 0.5 * float(domain.spacing_um)
    maximum_radius_um = float(np.nanmax(np.asarray(raster.radius_um)[lumen]))
    maximum_steps = max(16, int(math.ceil(8.0 * maximum_radius_um / step_um)))
    maximum_speed = float(np.nanmax(np.linalg.norm(velocity[lumen], axis=-1)))
    minimum_speed = max(64.0 * np.finfo(float).eps * maximum_speed, np.finfo(float).tiny)

    junction_indices = np.ascontiguousarray(np.argwhere(junction_cells), dtype=np.int64)
    if use_numba and njit is not None:
        traced_unit_ids = _trace_junction_cells_numba(
            np.ascontiguousarray(domain.x_coordinates_um, dtype=np.float64),
            np.ascontiguousarray(domain.z_coordinates_um, dtype=np.float64),
            float(domain.spacing_um),
            np.ascontiguousarray(lumen, dtype=np.bool_),
            np.ascontiguousarray(junction_cells, dtype=np.bool_),
            np.ascontiguousarray(ownership, dtype=np.int32),
            np.ascontiguousarray(velocity, dtype=np.float64),
            junction_indices,
            step_um,
            maximum_steps,
            minimum_speed,
        )
    else:
        traced_unit_ids = _trace_junction_cells_python(
            domain,
            lumen,
            junction_cells,
            ownership,
            velocity,
            junction_indices,
            step_um,
            maximum_steps,
            minimum_speed,
        )

    traced_successfully = traced_unit_ids >= 0
    successful_indices = junction_indices[traced_successfully]
    resolved[successful_indices[:, 0], successful_indices[:, 1]] = traced_unit_ids[
        traced_successfully
    ]
    unresolved_indices = junction_indices[~traced_successfully]

    if unresolved_indices.shape[0] > 0:
        traced = junction_cells & (resolved >= 0)
        if not np.any(traced):
            raise ValueError(
                "The accepted flow could not assign any bifurcation cell to a downstream vessel."
            )
        nearest = ndimage.distance_transform_edt(
            ~traced,
            return_distances=False,
            return_indices=True,
        )
        unresolved_x = unresolved_indices[:, 0]
        unresolved_z = unresolved_indices[:, 1]
        resolved[unresolved_x, unresolved_z] = resolved[
            nearest[0, unresolved_x, unresolved_z],
            nearest[1, unresolved_x, unresolved_z],
        ]
    return resolved, int(unresolved_indices.shape[0])


def _trace_junction_cells_python(
    domain: GridDomain,
    lumen: np.ndarray,
    junction: np.ndarray,
    ownership: np.ndarray,
    velocity: np.ndarray,
    junction_indices: np.ndarray,
    step_um: float,
    maximum_steps: int,
    minimum_speed: float,
) -> np.ndarray:
    """Reference implementation used for exact accelerated-path comparisons."""

    unit_ids = np.full(junction_indices.shape[0], -1, dtype=np.int32)
    for lane, (ix, iz) in enumerate(junction_indices):
        unit_ids[lane] = _trace_to_downstream_unit(
            domain,
            lumen,
            junction,
            ownership,
            velocity,
            int(ix),
            int(iz),
            step_um,
            maximum_steps,
            minimum_speed,
        )
    return unit_ids


def _trace_junction_cells_kernel(
    x_coordinates_um: np.ndarray,
    z_coordinates_um: np.ndarray,
    spacing_um: float,
    lumen: np.ndarray,
    junction: np.ndarray,
    ownership: np.ndarray,
    velocity: np.ndarray,
    junction_indices: np.ndarray,
    step_um: float,
    maximum_steps: int,
    minimum_speed: float,
) -> np.ndarray:
    """Trace independent junction cells in one allocation-free numeric batch."""

    unit_ids = np.full(junction_indices.shape[0], -1, dtype=np.int32)
    x_origin_um = float(x_coordinates_um[0])
    z_origin_um = float(z_coordinates_um[0])
    nx, nz = lumen.shape
    for lane in prange(junction_indices.shape[0]):
        start_ix = int(junction_indices[lane, 0])
        start_iz = int(junction_indices[lane, 1])
        position_x_um = float(x_coordinates_um[start_ix])
        position_z_um = float(z_coordinates_um[start_iz])
        for _ in range(maximum_steps):
            grid_x = (position_x_um - x_origin_um) / spacing_um
            grid_z = (position_z_um - z_origin_um) / spacing_um
            i0 = int(math.floor(grid_x))
            j0 = int(math.floor(grid_z))
            i0 = min(max(i0, 0), velocity.shape[0] - 2)
            j0 = min(max(j0, 0), velocity.shape[1] - 2)
            tx = min(max(grid_x - i0, 0.0), 1.0)
            tz = min(max(grid_z - j0, 0.0), 1.0)
            w00 = (1.0 - tx) * (1.0 - tz)
            w10 = tx * (1.0 - tz)
            w01 = (1.0 - tx) * tz
            w11 = tx * tz
            velocity_x = (
                w00 * velocity[i0, j0, 0]
                + w10 * velocity[i0 + 1, j0, 0]
                + w01 * velocity[i0, j0 + 1, 0]
                + w11 * velocity[i0 + 1, j0 + 1, 0]
            )
            velocity_z = (
                w00 * velocity[i0, j0, 1]
                + w10 * velocity[i0 + 1, j0, 1]
                + w01 * velocity[i0, j0 + 1, 1]
                + w11 * velocity[i0 + 1, j0 + 1, 1]
            )
            speed = math.sqrt(velocity_x * velocity_x + velocity_z * velocity_z)
            if not math.isfinite(speed) or speed <= minimum_speed:
                break
            position_x_um += step_um * velocity_x / speed
            position_z_um += step_um * velocity_z / speed
            ix = int(np.rint((position_x_um - x_origin_um) / spacing_um))
            iz = int(np.rint((position_z_um - z_origin_um) / spacing_um))
            if ix < 0 or iz < 0 or ix >= nx or iz >= nz:
                break
            if not lumen[ix, iz]:
                break
            if not junction[ix, iz] and ownership[ix, iz] >= 0:
                unit_ids[lane] = ownership[ix, iz]
                break
    return unit_ids


if njit is not None:
    _trace_junction_cells_numba = njit(cache=True, parallel=True)(
        _trace_junction_cells_kernel
    )
else:  # pragma: no cover - retained for deliberately minimal environments
    _trace_junction_cells_numba = _trace_junction_cells_kernel


def _trace_to_downstream_unit(
    domain: GridDomain,
    lumen: np.ndarray,
    junction: np.ndarray,
    ownership: np.ndarray,
    velocity: np.ndarray,
    start_ix: int,
    start_iz: int,
    step_um: float,
    maximum_steps: int,
    minimum_speed: float,
) -> int:
    position = np.asarray(
        [domain.x_coordinates_um[start_ix], domain.z_coordinates_um[start_iz]],
        dtype=np.float64,
    )
    for _ in range(maximum_steps):
        local_velocity = _sample_bilinear_vector(domain, velocity, position)
        speed = float(np.linalg.norm(local_velocity))
        if not np.isfinite(speed) or speed <= minimum_speed:
            return -1
        position += step_um * local_velocity / speed
        ix, iz = _nearest_cell(domain, position)
        if ix < 0 or iz < 0 or ix >= lumen.shape[0] or iz >= lumen.shape[1]:
            return -1
        if not lumen[ix, iz]:
            return -1
        if not junction[ix, iz] and ownership[ix, iz] >= 0:
            return int(ownership[ix, iz])
    return -1


def _sample_bilinear_vector(
    domain: GridDomain,
    values: np.ndarray,
    position_xz_um: np.ndarray,
) -> np.ndarray:
    x = (float(position_xz_um[0]) - float(domain.x_coordinates_um[0])) / float(domain.spacing_um)
    z = (float(position_xz_um[1]) - float(domain.z_coordinates_um[0])) / float(domain.spacing_um)
    i0 = int(np.floor(x))
    j0 = int(np.floor(z))
    i0 = min(max(i0, 0), values.shape[0] - 2)
    j0 = min(max(j0, 0), values.shape[1] - 2)
    tx = min(max(x - i0, 0.0), 1.0)
    tz = min(max(z - j0, 0.0), 1.0)
    return (
        (1.0 - tx) * (1.0 - tz) * values[i0, j0]
        + tx * (1.0 - tz) * values[i0 + 1, j0]
        + (1.0 - tx) * tz * values[i0, j0 + 1]
        + tx * tz * values[i0 + 1, j0 + 1]
    )


def _nearest_cell(domain: GridDomain, position_xz_um: np.ndarray) -> tuple[int, int]:
    ix = int(
        np.rint(
            (float(position_xz_um[0]) - float(domain.x_coordinates_um[0]))
            / float(domain.spacing_um)
        )
    )
    iz = int(
        np.rint(
            (float(position_xz_um[1]) - float(domain.z_coordinates_um[0]))
            / float(domain.spacing_um)
        )
    )
    return ix, iz


def _extend_ownership_to_solid_wall(
    ownership: np.ndarray,
    lumen: np.ndarray,
    solid_wall: np.ndarray,
) -> np.ndarray:
    extended = ownership.copy()
    valid = lumen & (ownership >= 0)
    nearest = ndimage.distance_transform_edt(
        ~valid,
        return_distances=False,
        return_indices=True,
    )
    extended[solid_wall] = ownership[nearest[0][solid_wall], nearest[1][solid_wall]]
    return extended


def _unit_mean_wall_shear(
    raster: RasterizedVessels,
    flow: FlowField,
    ownership: np.ndarray,
    topology: VesselBedTopology,
) -> np.ndarray:
    wss = np.asarray(flow.wall_shear_stress_pa, dtype=np.float64)
    closed_wall = np.asarray(raster.wall_mask, dtype=bool)
    means = np.full(len(topology.units), np.nan, dtype=np.float64)
    for unit in topology.units:
        selected = closed_wall & (ownership == unit.unit_id) & np.isfinite(wss)
        if np.any(selected):
            means[unit.unit_id] = float(np.mean(wss[selected]))
    return means


def _build_candidates(
    topology: VesselBedTopology,
    unit_mean_wss: np.ndarray,
    *,
    network_endothelial_wall_area_um2: float,
    network_inlet_flow_um3_s: float,
    injection_rate_per_s: float,
    observation_time_s: float,
) -> tuple[MolecularTargetCandidate, ...]:
    candidates: list[MolecularTargetCandidate] = []
    subtree_ids = {
        unit.unit_id: f"subtree:{unit.segment_ids[0]}"
        for unit in topology.units
        if unit.perfused
        and unit.parent_unit_id >= 0
        and bool(topology.descendants(unit.unit_id))
    }

    for unit in topology.units:
        if not unit.perfused:
            continue
        depth = int(unit.topology_depth)
        ancestor_parent = _nearest_parent_subtree_id(topology.units, unit, subtree_ids)
        descendants = topology.descendants(unit.unit_id)
        if unit.unit_id in subtree_ids:
            members = (unit.unit_id, *descendants)
            candidates.append(
                _candidate_from_units(
                    topology,
                    unit_mean_wss,
                    candidate_id=subtree_ids[unit.unit_id],
                    kind="downstream_subtree",
                    label=f"Downstream bed from vessel {unit.segment_ids[0]}",
                    root_unit_id=unit.unit_id,
                    member_unit_ids=members,
                    parent_candidate_id=ancestor_parent,
                    depth=2 * depth,
                    network_endothelial_wall_area_um2=network_endothelial_wall_area_um2,
                    network_inlet_flow_um3_s=network_inlet_flow_um3_s,
                    injection_rate_per_s=injection_rate_per_s,
                    observation_time_s=observation_time_s,
                )
            )
        candidates.append(
            _candidate_from_units(
                topology,
                unit_mean_wss,
                candidate_id=f"segment:{unit.segment_ids[0]}",
                kind="single_vessel_unit",
                label=f"Vessel unit {unit.segment_ids[0]}",
                root_unit_id=unit.unit_id,
                member_unit_ids=(unit.unit_id,),
                parent_candidate_id=subtree_ids.get(unit.unit_id, ancestor_parent),
                depth=2 * depth + 1,
                network_endothelial_wall_area_um2=network_endothelial_wall_area_um2,
                network_inlet_flow_um3_s=network_inlet_flow_um3_s,
                injection_rate_per_s=injection_rate_per_s,
                observation_time_s=observation_time_s,
            )
        )

    candidates.sort(key=lambda item: (item.depth, item.root_unit_id, item.kind))
    return tuple(candidates)


def _candidate_from_units(
    topology: VesselBedTopology,
    unit_mean_wss: np.ndarray,
    *,
    candidate_id: str,
    kind: str,
    label: str,
    root_unit_id: int,
    member_unit_ids: tuple[int, ...],
    parent_candidate_id: str | None,
    depth: int,
    network_endothelial_wall_area_um2: float,
    network_inlet_flow_um3_s: float,
    injection_rate_per_s: float,
    observation_time_s: float,
) -> MolecularTargetCandidate:
    root_unit = topology.units[root_unit_id]
    tree_root = topology.units[root_unit.root_unit_id]
    volume = float(sum(topology.units[index].volume_um3 for index in member_unit_ids))
    flow = float(root_unit.flow_rate_um3_s)
    flow_fraction = flow / tree_root.flow_rate_um3_s if tree_root.flow_rate_um3_s > 0.0 else math.nan
    residence = volume / flow if flow > 0.0 else math.inf
    wall_areas = np.asarray(
        [topology.units[index].endothelial_wall_area_um2 for index in member_unit_ids],
        dtype=np.float64,
    )
    shear = np.asarray([unit_mean_wss[index] for index in member_unit_ids], dtype=np.float64)
    valid = np.isfinite(shear) & (wall_areas > 0.0)
    mean_wss = (
        float(np.sum(wall_areas[valid] * shear[valid]) / np.sum(wall_areas[valid]))
        if np.any(valid)
        else math.nan
    )
    wall_area = float(np.sum(wall_areas))
    centroid_x = float(
        sum(
            topology.units[index].endothelial_wall_area_um2
            * topology.units[index].wall_area_centroid_x_um
            for index in member_unit_ids
        )
        / wall_area
    )
    centroid_z = float(
        sum(
            topology.units[index].endothelial_wall_area_um2
            * topology.units[index].wall_area_centroid_z_um
            for index in member_unit_ids
        )
        / wall_area
    )
    second_moment = float(
        sum(topology.units[index].wall_area_second_moment_um4 for index in member_unit_ids)
    )
    radius_of_gyration = math.sqrt(
        max(second_moment / wall_area - centroid_x * centroid_x - centroid_z * centroid_z, 0.0)
    )
    wall_area_fraction = (
        wall_area / network_endothelial_wall_area_um2
        if network_endothelial_wall_area_um2 > 0.0
        else math.nan
    )
    network_flow_fraction = (
        flow / network_inlet_flow_um3_s if network_inlet_flow_um3_s > 0.0 else math.nan
    )
    expected_visits = (
        injection_rate_per_s * observation_time_s * network_flow_fraction
        if math.isfinite(injection_rate_per_s)
        and math.isfinite(observation_time_s)
        and math.isfinite(network_flow_fraction)
        else math.nan
    )
    return MolecularTargetCandidate(
        candidate_id=candidate_id,
        kind=kind,
        label=label,
        root_unit_id=root_unit_id,
        member_unit_ids=tuple(int(value) for value in member_unit_ids),
        parent_candidate_id=parent_candidate_id,
        depth=int(depth),
        topology_depth=int(root_unit.topology_depth),
        inlet_flow_um3_s=flow,
        root_flow_fraction=float(flow_fraction),
        network_flow_fraction=float(network_flow_fraction),
        volume_um3=volume,
        residence_time_s=float(residence),
        mean_wall_shear_pa=mean_wss,
        endothelial_wall_area_um2=wall_area,
        endothelial_wall_area_fraction=float(wall_area_fraction),
        wall_area_centroid_x_um=centroid_x,
        wall_area_centroid_z_um=centroid_z,
        radius_of_gyration_um=float(radius_of_gyration),
        expected_bubble_visits=float(expected_visits),
    )


def _validate_automatic_selection_context(
    injection_rate_per_s: float,
    observation_time_s: float,
) -> None:
    injection_is_missing = math.isnan(float(injection_rate_per_s))
    duration_is_missing = math.isnan(float(observation_time_s))
    if injection_is_missing != duration_is_missing:
        raise ValueError(
            "Automatic target accessibility requires both injection_rate_per_s and "
            "observation_time_s, or neither value."
        )
    if injection_is_missing:
        return
    if not math.isfinite(float(injection_rate_per_s)) or injection_rate_per_s <= 0.0:
        raise ValueError("injection_rate_per_s must be finite and greater than zero.")
    if not math.isfinite(float(observation_time_s)) or observation_time_s <= 0.0:
        raise ValueError("observation_time_s must be finite and greater than zero.")


def _nearest_parent_subtree_id(
    units: tuple[VesselBedUnit, ...],
    unit: VesselBedUnit,
    subtree_ids: dict[int, str],
) -> str | None:
    parent_id = unit.parent_unit_id
    while parent_id >= 0:
        if parent_id in subtree_ids:
            return subtree_ids[parent_id]
        parent_id = units[parent_id].parent_unit_id
    return None
