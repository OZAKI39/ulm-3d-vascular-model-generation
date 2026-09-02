"""Periodic cardiac-flow waveform shared by both trajectory pipelines."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .signal import generate_ecg_normalized


@dataclass(frozen=True)
class PeriodicCardiacWaveform:
    """One positive, periodically interpolated cardiac-flow multiplier."""

    bpm: float
    period_s: float
    sample_time_s: np.ndarray
    multiplier: np.ndarray
    cumulative_integral_s: np.ndarray

    @property
    def sample_dt_s(self) -> float:
        return float(self.period_s / self.multiplier.size)

    @property
    def cycle_mean(self) -> float:
        return float(self.cumulative_integral_s[-1] / self.period_s)

    def evaluate(self, time_s: float | np.ndarray) -> float | np.ndarray:
        """Evaluate the periodic piecewise-linear waveform at arbitrary times."""

        values, _, scalar = self._interval_coordinates(time_s)
        left, right, fraction = values
        result = self.multiplier[left] + fraction * (
            self.multiplier[right] - self.multiplier[left]
        )
        return float(result) if scalar else result

    def derivative_s_inv(self, time_s: float | np.ndarray) -> float | np.ndarray:
        """Return the time derivative of the piecewise-linear multiplier."""

        values, _, scalar = self._interval_coordinates(time_s)
        left, right, _ = values
        result = (self.multiplier[right] - self.multiplier[left]) / self.sample_dt_s
        return float(result) if scalar else result

    def primitive_s(self, time_s: float | np.ndarray) -> float | np.ndarray:
        """Integrate the periodic multiplier from zero to any signed time."""

        values, cycles, scalar = self._interval_coordinates(time_s)
        left, right, fraction = values
        v0 = self.multiplier[left]
        v1 = self.multiplier[right]
        local = self.cumulative_integral_s[left] + self.sample_dt_s * (
            v0 * fraction + 0.5 * (v1 - v0) * fraction * fraction
        )
        result = cycles * self.cumulative_integral_s[-1] + local
        return float(result) if scalar else result

    def integrate_s(
        self,
        start_time_s: float | np.ndarray,
        end_time_s: float | np.ndarray,
        *,
        phase_offset_s: float = 0.0,
    ) -> float | np.ndarray:
        """Integrate the multiplier between two times with a fixed phase offset."""

        start = np.asarray(start_time_s, dtype=np.float64) + float(phase_offset_s)
        end = np.asarray(end_time_s, dtype=np.float64) + float(phase_offset_s)
        result = np.asarray(self.primitive_s(end)) - np.asarray(self.primitive_s(start))
        return float(result) if result.ndim == 0 else result

    def _interval_coordinates(
        self, time_s: float | np.ndarray
    ) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray, bool]:
        time = np.asarray(time_s, dtype=np.float64)
        scalar = time.ndim == 0
        cycles = np.floor(time / self.period_s).astype(np.int64)
        remainder = time - cycles * self.period_s
        remainder = np.clip(remainder, 0.0, np.nextafter(self.period_s, 0.0))
        coordinate = remainder / self.sample_dt_s
        left = np.floor(coordinate).astype(np.int64)
        left = np.clip(left, 0, self.multiplier.size - 1)
        right = (left + 1) % self.multiplier.size
        fraction = coordinate - left
        return (left, right, fraction), cycles, scalar


def build_periodic_ecg_flow_waveform(
    bpm: float,
    samples_per_cycle: int,
    *,
    preserve_cycle_mean_flow: bool = True,
    modulation_strength: float = 1.0,
) -> PeriodicCardiacWaveform:
    """Build the legacy ECG-shaped surrogate as a compact periodic flow waveform."""

    heart_rate = float(bpm)
    if not math.isfinite(heart_rate) or heart_rate <= 0.0:
        raise ValueError("cardiac_pulsatility.bpm must be finite and positive.")
    if isinstance(samples_per_cycle, bool) or not isinstance(samples_per_cycle, int):
        raise ValueError("cardiac_pulsatility.waveform_samples_per_cycle must be an integer.")
    count = int(samples_per_cycle)
    if count < 32:
        raise ValueError(
            "cardiac_pulsatility.waveform_samples_per_cycle must be at least 32."
        )
    strength = float(modulation_strength)
    if not math.isfinite(strength) or not 0.0 <= strength <= 1.0:
        raise ValueError("cardiac_pulsatility.modulation_strength must be in [0, 1].")

    period = 60.0 / heart_rate
    sample_dt = period / count
    sample_time = np.arange(count, dtype=np.float64) * sample_dt
    legacy = np.asarray(
        generate_ecg_normalized(heart_rate, sample_dt, 3.0 * period), dtype=np.float64
    )
    if legacy.size < 2 * count:
        raise RuntimeError("The synthetic ECG generator returned an incomplete cardiac cycle.")
    values = np.ascontiguousarray(legacy[count : 2 * count], dtype=np.float64)
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("The cardiac-flow multiplier must remain finite and positive.")
    if preserve_cycle_mean_flow:
        values = values / float(np.mean(values))
    # Scale only the deviation from the steady multiplier of one. This leaves the
    # cycle-mean CFD reference unchanged when mean preservation is enabled, while
    # strength zero gives a strictly steady flow and strength one keeps the full
    # synthetic waveform.
    values = 1.0 + strength * (values - 1.0)

    following = np.roll(values, -1)
    interval_area = 0.5 * (values + following) * sample_dt
    cumulative = np.concatenate(
        (np.asarray([0.0], dtype=np.float64), np.cumsum(interval_area))
    )
    return PeriodicCardiacWaveform(
        bpm=heart_rate,
        period_s=period,
        sample_time_s=sample_time,
        multiplier=values,
        cumulative_integral_s=cumulative,
    )
