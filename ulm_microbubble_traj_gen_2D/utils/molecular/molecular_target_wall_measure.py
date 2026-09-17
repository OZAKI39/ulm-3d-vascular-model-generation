"""Physical endothelial-area weights and accessibility on raster wall sites."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy import ndimage

from ulm_vascular_model_generator.utils.core.models import Vessel

from ..core.types import RasterizedVessels
from ..geometry.vessel_bed_topology import VesselBedTopology


MINIMUM_EXPECTED_BUBBLE_VISITS = 1.0


@dataclass(frozen=True)
class MolecularTargetWallMeasure:
    """One physical wall-area quadrature weight and access flag per grid site."""

    wall_area_weight_um2: np.ndarray
    wall_segment_id_grid: np.ndarray
    accessible_wall_mask: np.ndarray
    expected_bubble_visits_by_unit: np.ndarray
    mapped_endothelial_wall_area_um2: float
    unmapped_endothelial_wall_area_um2: float


def build_molecular_target_wall_measure(
    raster: RasterizedVessels,
    topology: VesselBedTopology,
    unit_id_grid: np.ndarray,
    solid_wall_mask: np.ndarray,
    open_boundary_mask: np.ndarray,
    vessels: list[Vessel],
    *,
    injection_rate_per_s: float,
    observation_time_s: float,
    network_inlet_flow_um3_s: float,
) -> MolecularTargetWallMeasure:
    """Map cylindrical segment area to closed raster walls without double counting."""

    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    solid_wall = np.asarray(solid_wall_mask, dtype=bool)
    open_boundary = np.asarray(open_boundary_mask, dtype=bool)
    ownership = np.asarray(unit_id_grid, dtype=np.int32)
    if lumen.shape != solid_wall.shape or lumen.shape != open_boundary.shape:
        raise ValueError("Molecular-target wall masks must share one grid shape.")
    if ownership.shape != lumen.shape:
        raise ValueError("Molecular-target unit ownership must match the wall grid shape.")
    closed_wall = solid_wall & ~open_boundary & (ownership >= 0)
    if not np.any(closed_wall):
        raise ValueError("Automatic molecular-target preparation found no closed wall sites.")

    nearest_lumen = ndimage.distance_transform_edt(
        ~lumen,
        return_distances=False,
        return_indices=True,
    )
    lumen_segment_ids = np.asarray(raster.vessel_id, dtype=np.int64)
    wall_segment_ids = np.full(lumen.shape, -1, dtype=np.int32)
    wall_segment_ids[closed_wall] = lumen_segment_ids[
        nearest_lumen[0][closed_wall],
        nearest_lumen[1][closed_wall],
    ].astype(np.int32)

    vessel_by_id = {int(vessel.vid): vessel for vessel in vessels}
    if len(vessel_by_id) != len(vessels):
        raise ValueError("Vessel IDs must be unique before wall area is distributed.")
    weights = np.zeros(lumen.shape, dtype=np.float64)
    missing_area_by_unit = np.zeros(len(topology.units), dtype=np.float64)
    for segment_id, vessel in sorted(vessel_by_id.items()):
        if segment_id not in topology.segment_id_to_unit_id:
            raise ValueError(f"Vessel {segment_id} has no molecular-target unit ownership.")
        segment_area = 2.0 * np.pi * float(vessel.radius) * float(vessel.length())
        selected = closed_wall & (wall_segment_ids == segment_id)
        count = int(np.count_nonzero(selected))
        if count:
            weights[selected] += segment_area / count
        else:
            unit_id = topology.segment_id_to_unit_id[segment_id]
            missing_area_by_unit[unit_id] += segment_area

    unmapped_area = 0.0
    for unit in topology.units:
        missing_area = float(missing_area_by_unit[unit.unit_id])
        if missing_area <= 0.0:
            continue
        selected = closed_wall & (ownership == unit.unit_id)
        count = int(np.count_nonzero(selected))
        if count:
            weights[selected] += missing_area / count
        else:
            unmapped_area += missing_area

    expected_visits = np.full(len(topology.units), np.nan, dtype=np.float64)
    accessible = np.zeros(lumen.shape, dtype=bool)
    context_available = bool(
        math.isfinite(injection_rate_per_s)
        and injection_rate_per_s > 0.0
        and math.isfinite(observation_time_s)
        and observation_time_s > 0.0
        and math.isfinite(network_inlet_flow_um3_s)
        and network_inlet_flow_um3_s > 0.0
    )
    if context_available:
        for unit in topology.units:
            flow_fraction = unit.flow_rate_um3_s / network_inlet_flow_um3_s
            expected = injection_rate_per_s * observation_time_s * flow_fraction
            expected_visits[unit.unit_id] = expected
            if unit.perfused and expected >= MINIMUM_EXPECTED_BUBBLE_VISITS:
                accessible |= closed_wall & (ownership == unit.unit_id)

    mapped_area = float(np.sum(weights[closed_wall]))
    return MolecularTargetWallMeasure(
        wall_area_weight_um2=np.ascontiguousarray(weights, dtype=np.float64),
        wall_segment_id_grid=np.ascontiguousarray(wall_segment_ids, dtype=np.int32),
        accessible_wall_mask=np.ascontiguousarray(accessible, dtype=bool),
        expected_bubble_visits_by_unit=np.ascontiguousarray(
            expected_visits,
            dtype=np.float64,
        ),
        mapped_endothelial_wall_area_um2=mapped_area,
        unmapped_endothelial_wall_area_um2=float(unmapped_area),
    )
