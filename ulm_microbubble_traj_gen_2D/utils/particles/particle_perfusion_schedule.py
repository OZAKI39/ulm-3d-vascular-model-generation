"""Time-step-independent deterministic events for continuous MB perfusion."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING

import numpy as np

from .particle_inlet_flux import InletFluxModel

if TYPE_CHECKING:
    from ..cardiac.cardiac_pulsatility import CardiacPulsatility


@dataclass(frozen=True)
class PerfusionSchedule:
    """All planned events up to a known numerical safety horizon."""

    planned_time_s: np.ndarray
    radius_um: np.ndarray
    position_grid: np.ndarray
    initial_vessel_id: np.ndarray

    @property
    def count(self) -> int:
        return int(self.planned_time_s.size)


def build_perfusion_schedule(
    model: InletFluxModel,
    horizon_s: float,
    cardiac: CardiacPulsatility | None = None,
) -> PerfusionSchedule:
    """Discretize ``j(s,R)`` with midpoint event times and a 2-D Halton sequence."""

    horizon = float(horizon_s)
    if not math.isfinite(horizon) or horizon < 0.0:
        raise ValueError("The perfusion schedule horizon must be finite and non-negative.")
    rate = float(model.injection_rate_per_s)
    if cardiac is None:
        cumulative_events = rate * horizon
    else:
        cumulative_events = rate * float(
            cardiac.integrate_inlet_multiplier_s(0.0, horizon)
        )
    count = max(0, int(math.floor(cumulative_events + 0.5)))
    ids = np.arange(count, dtype=np.int64)
    if cardiac is None:
        planned = (ids.astype(np.float64) + 0.5) / rate
    else:
        planned = _invert_cardiac_event_times(
            ids.astype(np.float64) + 0.5,
            rate,
            horizon,
            cardiac,
        )
    keep = planned <= horizon + 16.0 * np.finfo(float).eps * max(horizon, 1.0)
    planned = planned[keep]
    ids = ids[keep]
    radii = np.empty(ids.size, dtype=np.float64)
    positions = np.empty((ids.size, 2), dtype=np.float64)
    initial_vessel_id = np.empty(ids.size, dtype=np.int32)
    for row, event_id in enumerate(ids):
        eta = radical_inverse(int(event_id) + 1, 2)
        xi = radical_inverse(int(event_id) + 1, 3)
        radius = model.sample_radius_um(eta)
        radii[row] = radius
        topological_sampler = getattr(model, "sample_position_and_vessel_id", None)
        if callable(topological_sampler):
            positions[row], initial_vessel_id[row] = topological_sampler(radius, xi)
        else:
            positions[row] = model.sample_position_grid(radius, xi)
            initial_vessel_id[row] = 1
    return PerfusionSchedule(
        planned_time_s=planned,
        radius_um=radii,
        position_grid=positions,
        initial_vessel_id=initial_vessel_id,
    )


def _invert_cardiac_event_times(
    targets: np.ndarray,
    mean_rate_per_s: float,
    horizon_s: float,
    cardiac: CardiacPulsatility,
) -> np.ndarray:
    """Invert cumulative pulsatile number flux without tying events to a time step."""

    if targets.size == 0:
        return np.empty(0, dtype=np.float64)
    lower = np.zeros(targets.size, dtype=np.float64)
    upper = np.full(targets.size, float(horizon_s), dtype=np.float64)
    for _ in range(64):
        middle = 0.5 * (lower + upper)
        cumulative = mean_rate_per_s * np.asarray(
            cardiac.integrate_inlet_multiplier_s(0.0, middle), dtype=np.float64
        )
        before = cumulative < targets
        lower = np.where(before, middle, lower)
        upper = np.where(before, upper, middle)
    return 0.5 * (lower + upper)


def radical_inverse(index: int, base: int) -> float:
    """Return one unscrumbled Halton coordinate strictly inside ``(0, 1)``."""

    value = int(index)
    radix = int(base)
    if value <= 0:
        raise ValueError("Halton indices must start at one.")
    if radix < 2:
        raise ValueError("The radical-inverse base must be at least two.")
    inverse = 1.0 / radix
    factor = inverse
    result = 0.0
    while value:
        value, digit = divmod(value, radix)
        result += digit * factor
        factor *= inverse
    return float(result)
