"""Spatial propagation of a periodic cardiac-flow multiplier on the vessel grid."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy import ndimage

try:
    from numba import njit
except ImportError:  # pragma: no cover - the production environment includes Numba
    njit = None

from ulm_vascular_model_generator.utils.core.models import Vessel

from .cardiac_waveform import PeriodicCardiacWaveform, build_periodic_ecg_flow_waveform
from ..core.config import CardiacPulsatilityConfig
from ..particles.particle_field_sampling import sample_bilinear
from ..core.types import GridDomain, RasterizedVessels


@dataclass(frozen=True)
class CardiacSample:
    """Local multiplier and its spatial gradient at one simulation time."""

    multiplier: np.ndarray
    gradient_per_um: np.ndarray


@dataclass(frozen=True)
class CardiacPulsatility:
    """Compact cardiac waveform plus root-to-grid propagation delays."""

    waveform: PeriodicCardiacWaveform
    phase_offset_s: float
    path_distance_um: np.ndarray
    delay_s: np.ndarray
    delay_gradient_s_per_um: np.ndarray
    waveform_name: str
    preserve_cycle_mean_flow: bool
    modulation_strength: float
    pulse_propagation_velocity_um_s: float

    def sample(
        self,
        positions_grid: np.ndarray,
        time_s: float,
        *,
        use_numba: bool = False,
    ) -> CardiacSample:
        """Evaluate retarded cardiac phase continuously at fractional grid positions."""

        positions = np.asarray(positions_grid, dtype=np.float64)
        if use_numba and njit is not None:
            multiplier, gradient = _sample_cardiac_numba(
                np.ascontiguousarray(positions),
                float(time_s),
                float(self.phase_offset_s),
                np.ascontiguousarray(self.delay_s, dtype=np.float64),
                np.ascontiguousarray(self.delay_gradient_s_per_um, dtype=np.float64),
                float(self.waveform.period_s),
                np.ascontiguousarray(self.waveform.multiplier, dtype=np.float64),
            )
            return CardiacSample(multiplier=multiplier, gradient_per_um=gradient)
        delay = sample_bilinear(self.delay_s, positions)
        delay_gradient = sample_bilinear(self.delay_gradient_s_per_um, positions)
        phase_time = float(time_s) + self.phase_offset_s - delay
        multiplier = np.asarray(self.waveform.evaluate(phase_time), dtype=np.float64)
        time_derivative = np.asarray(
            self.waveform.derivative_s_inv(phase_time), dtype=np.float64
        )
        gradient = -time_derivative[:, None] * delay_gradient
        return CardiacSample(multiplier=multiplier, gradient_per_um=gradient)

    def integrate_inlet_multiplier_s(
        self, start_time_s: float | np.ndarray, end_time_s: float | np.ndarray
    ) -> float | np.ndarray:
        """Integrate inlet flow multiplier for a non-uniform perfusion schedule."""

        return self.waveform.integrate_s(
            start_time_s,
            end_time_s,
            phase_offset_s=self.phase_offset_s,
        )


def build_cardiac_pulsatility(
    domain: GridDomain,
    raster: RasterizedVessels,
    vessels: list[Vessel],
    cfg: CardiacPulsatilityConfig | None,
) -> CardiacPulsatility | None:
    """
    Build the enabled cardiac model or return ``None`` for steady transport.
    """

    # First check whether the user actually wants a heartbeat effect.
    if cfg is None or not bool(cfg.enabled):
        return None
    
    # A heartbeat signal needs a real vessel tree along which it can travel.
    if not vessels:
        raise ValueError("Enabled cardiac pulsatility requires a non-empty vessel tree.")

    # Build one complete heartbeat curve that can be repeated throughout the simulation.
    waveform = build_periodic_ecg_flow_waveform(
        cfg.bpm,
        cfg.waveform_samples_per_cycle,
        preserve_cycle_mean_flow=cfg.preserve_cycle_mean_flow,
        modulation_strength=cfg.modulation_strength,
    )

    # Measure how far every grid position lies from a root inlet when travelling through the vessel tree.
    path_distance       = _root_path_distance_field(domain, raster, vessels)

    # Convert distance into travel time with the familiar rule: time = distance / speed.
    # Farther positions therefore feel the same heartbeat later than positions close to the root inlet.
    delay               = path_distance / float(cfg.pulse_propagation_velocity_um_s)

    # Read the real physical width of one grid cell so nearby delay values can be compared in micrometres.
    spacing             = float(domain.spacing_um)

    # Find how much the arrival delay changes when moving in the X or Z direction.
    # A "gradient" is simply a small map showing which way a value increases and how quickly it changes.
    # This lets later velocity and shear calculations respect the different heartbeat timing of nearby cells.
    delay_dx, delay_dz  = np.gradient(delay, spacing, spacing, edge_order=1)

    # Join the X and Z delay changes into one two-number direction for every grid cell.
    delay_gradient      = np.stack((delay_dx, delay_dz), axis=-1)
    
    return CardiacPulsatility(
        waveform=waveform,
        phase_offset_s=float(cfg.initial_phase_fraction) * waveform.period_s,
        path_distance_um=np.asarray(path_distance, dtype=np.float32),
        delay_s=np.asarray(delay, dtype=np.float64),
        delay_gradient_s_per_um=np.asarray(delay_gradient, dtype=np.float64),
        waveform_name=str(cfg.waveform),
        preserve_cycle_mean_flow=bool(cfg.preserve_cycle_mean_flow),
        modulation_strength=float(cfg.modulation_strength),
        pulse_propagation_velocity_um_s=float(cfg.pulse_propagation_velocity_um_s),
    )


def _root_path_distance_field(
    domain: GridDomain,
    raster: RasterizedVessels,
    vessels: list[Vessel],
) -> np.ndarray:
    vessel_by_id = {int(vessel.vid): vessel for vessel in vessels}
    if len(vessel_by_id) != len(vessels):
        raise ValueError("Cardiac propagation requires unique vessel IDs.")
    prefix_cache: dict[int, float] = {}
    visiting: set[int] = set()

    def proximal_distance(vessel_id: int) -> float:
        if vessel_id in prefix_cache:
            return prefix_cache[vessel_id]
        if vessel_id in visiting:
            raise ValueError("Cardiac propagation found a cycle in the vessel parent graph.")
        vessel = vessel_by_id.get(vessel_id)
        if vessel is None:
            raise ValueError(f"Cardiac propagation cannot find vessel ID {vessel_id}.")
        visiting.add(vessel_id)
        parent_id = int(vessel.parent_id)
        if parent_id < 0:
            distance = 0.0
        else:
            parent = vessel_by_id.get(parent_id)
            if parent is None:
                raise ValueError(
                    f"Cardiac propagation cannot find parent vessel ID {parent_id}."
                )
            distance = proximal_distance(parent_id) + _projected_vessel_length(parent)
        visiting.remove(vessel_id)
        prefix_cache[vessel_id] = distance
        return distance

    x, z = np.meshgrid(
        np.asarray(domain.x_coordinates_um, dtype=np.float64),
        np.asarray(domain.z_coordinates_um, dtype=np.float64),
        indexing="ij",
    )
    distance = np.full(domain.shape, np.nan, dtype=np.float64)
    vessel_ids = np.asarray(raster.vessel_id, dtype=np.int64)
    attributed = vessel_ids >= 0
    unique_ids = np.unique(vessel_ids[attributed])
    proximal_xz = np.empty((unique_ids.size, 2), dtype=np.float64)
    direction_xz = np.empty((unique_ids.size, 2), dtype=np.float64)
    segment_length = np.empty(unique_ids.size, dtype=np.float64)
    prefix_distance = np.empty(unique_ids.size, dtype=np.float64)
    for row, vessel_id in enumerate(unique_ids):
        vessel = vessel_by_id.get(int(vessel_id))
        if vessel is None:
            raise ValueError(
                f"Rasterized vessel ID {int(vessel_id)} is absent from the vessel tree."
            )
        proximal = np.asarray(vessel.x_p, dtype=np.float64)[[0, 2]]
        distal = np.asarray(vessel.x_d, dtype=np.float64)[[0, 2]]
        axis = distal - proximal
        length = _projected_vessel_length(vessel)
        if not math.isfinite(length) or length <= 0.0:
            raise ValueError(f"Vessel {int(vessel_id)} has non-positive projected length.")
        proximal_xz[row] = proximal
        direction_xz[row] = axis / length
        segment_length[row] = length
        prefix_distance[row] = proximal_distance(int(vessel_id))

    flat_ids = vessel_ids[attributed]
    rows = np.searchsorted(unique_ids, flat_ids)
    local = (
        (x[attributed] - proximal_xz[rows, 0]) * direction_xz[rows, 0]
        + (z[attributed] - proximal_xz[rows, 1]) * direction_xz[rows, 1]
    )
    distance[attributed] = prefix_distance[rows] + np.clip(
        local, 0.0, segment_length[rows]
    )

    source = np.isfinite(distance) & np.asarray(raster.lumen_mask, dtype=bool)
    if not np.any(source):
        raise ValueError("Cardiac propagation could not attribute any lumen grid cell.")
    nearest = ndimage.distance_transform_edt(~source, return_distances=False, return_indices=True)
    filled = distance[tuple(nearest)]
    if not np.all(np.isfinite(filled)):
        raise ValueError("Cardiac root-path distance field contains non-finite values.")
    return np.asarray(filled, dtype=np.float64)


def _projected_vessel_length(vessel: Vessel) -> float:
    proximal = np.asarray(vessel.x_p, dtype=np.float64)[[0, 2]]
    distal = np.asarray(vessel.x_d, dtype=np.float64)[[0, 2]]
    return float(np.linalg.norm(distal - proximal))


def _sample_cardiac_kernel(
    positions,
    time_s,
    phase_offset_s,
    delay_field,
    delay_gradient_field,
    period_s,
    waveform_multiplier,
):
    count = positions.shape[0]
    nx = delay_field.shape[0]
    nz = delay_field.shape[1]
    sample_count = waveform_multiplier.size
    sample_dt_s = period_s / sample_count
    multiplier = np.empty(count, dtype=np.float64)
    gradient = np.empty((count, 2), dtype=np.float64)
    for lane in range(count):
        gx = positions[lane, 0]
        gz = positions[lane, 1]
        if gx < 0.0:
            gx = 0.0
        elif gx > nx - 1.0:
            gx = nx - 1.0
        if gz < 0.0:
            gz = 0.0
        elif gz > nz - 1.0:
            gz = nz - 1.0
        i0 = int(np.floor(gx))
        j0 = int(np.floor(gz))
        i1 = min(i0 + 1, nx - 1)
        j1 = min(j0 + 1, nz - 1)
        tx = gx - i0
        tz = gz - j0
        w00 = (1.0 - tx) * (1.0 - tz)
        w10 = tx * (1.0 - tz)
        w01 = (1.0 - tx) * tz
        w11 = tx * tz
        delay = (
            w00 * delay_field[i0, j0]
            + w10 * delay_field[i1, j0]
            + w01 * delay_field[i0, j1]
            + w11 * delay_field[i1, j1]
        )
        phase_time = time_s + phase_offset_s - delay
        cycles = np.floor(phase_time / period_s)
        remainder = phase_time - cycles * period_s
        if remainder < 0.0:
            remainder = 0.0
        elif remainder >= period_s:
            remainder = np.nextafter(period_s, 0.0)
        coordinate = remainder / sample_dt_s
        left = int(np.floor(coordinate))
        if left < 0:
            left = 0
        elif left >= sample_count:
            left = sample_count - 1
        right = (left + 1) % sample_count
        fraction = coordinate - left
        left_value = waveform_multiplier[left]
        right_value = waveform_multiplier[right]
        multiplier[lane] = left_value + fraction * (right_value - left_value)
        derivative = (right_value - left_value) / sample_dt_s
        for component in range(2):
            delay_gradient = (
                w00 * delay_gradient_field[i0, j0, component]
                + w10 * delay_gradient_field[i1, j0, component]
                + w01 * delay_gradient_field[i0, j1, component]
                + w11 * delay_gradient_field[i1, j1, component]
            )
            gradient[lane, component] = -derivative * delay_gradient
    return multiplier, gradient


if njit is not None:
    _sample_cardiac_numba = njit(cache=True)(_sample_cardiac_kernel)
else:  # pragma: no cover - retained for deliberately minimal environments
    _sample_cardiac_numba = _sample_cardiac_kernel
