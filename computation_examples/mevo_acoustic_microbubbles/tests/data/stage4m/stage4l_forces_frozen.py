"""Harmonic primary and quasistatic far-field secondary Bjerknes forces in N."""

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


def primary_bjerknes(volume_response_m3: np.ndarray, pressure_gradient_Pa_m: np.ndarray) -> np.ndarray:
    return -0.5 * np.real(np.asarray(volume_response_m3)[..., None] * np.conj(pressure_gradient_Pa_m))


class PairInteractionBackend(ABC):
    @abstractmethod
    def forces(self, positions_m, radii_m, volume_response_m3, frequency_Hz, density_kg_m3):
        """Net force in N for every bubble; interaction is reciprocal."""


@dataclass
class LegacyRegularizedBjerknes(PairInteractionBackend):
    """DEPRECATED / DEMO ONLY: Stage-1 distance clipping, no near-field physics."""

    minimum_separation_m: float
    collision_regularization_m: float
    cutoff_m: float

    def forces(self, positions_m, radii_m, volume_response_m3, frequency_Hz, density_kg_m3):
        delta = positions_m[None, :, :] - positions_m[:, None, :]
        d = np.linalg.norm(delta, axis=-1)
        valid = (d > 0) & (d <= self.cutoff_m)
        minimum = np.maximum(self.minimum_separation_m, radii_m[:, None] + radii_m[None, :])
        safe = np.maximum(d, minimum + self.collision_regularization_m)
        dotV = (2j * np.pi * frequency_Hz) * volume_response_m3
        correlation = 0.5 * np.real(dotV[:, None] * np.conj(dotV[None, :]))
        magnitude = density_kg_m3 / (4 * np.pi * safe**2) * correlation * valid
        direction = delta / np.maximum(d[..., None], np.finfo(float).tiny)
        return np.sum(magnitude[..., None] * direction, axis=1)


# Import compatibility for the preserved Stage-1 scripts and regression tests.
NumpyPairBackend = LegacyRegularizedBjerknes


class ClassicalFarFieldBjerknes(PairInteractionBackend):
    """Unregularized cycle-average point-monopole force, SI, reciprocal.

    F_i = rho/(4*pi*d²) <dot(V_i) dot(V_j)> e_ij; complex convention exp(+iwt).
    Optional N x H volume harmonics retain the h² derivative weighting.
    Near-field policy and oscillation/retardation limits are checked by the
    shared monitor. Strict use rejects close pairs; warning-only population use
    delegates nonpenetration to the geometric constraints, without changing force.
    """

    def __init__(self, backend="numpy", min_distance_over_sum=5.0, max_bubbles=2000, enforce_far_field=True):
        if backend == "cupy":
            import cupy as xp
        elif backend == "numpy":
            xp = np
        else:
            raise ValueError("Pair backend must be numpy or cupy")
        if not np.isfinite(min_distance_over_sum) or min_distance_over_sum < 1:
            raise ValueError("Distance criterion must be at least non-overlap")
        self.xp, self.min_ratio, self.max_bubbles = xp, min_distance_over_sum, max_bubbles
        self.enforce_far_field = enforce_far_field

    def forces(self, positions_m, radii_m, volume_response_m3, frequency_Hz, density_kg_m3):
        xp = self.xp
        x, r, v = xp.asarray(positions_m), xp.asarray(radii_m), xp.asarray(volume_response_m3)
        if len(x) > self.max_bubbles:
            raise ValueError("All-pairs size budget exceeded; benchmark before increasing max_bubbles")
        if x.shape != (len(r), 3) or v.shape[0] != len(r) or v.ndim not in (1, 2):
            raise ValueError("Expected positions[N,3], radii[N], volumes[N] or [N,H]")
        if (
            not bool(xp.all(xp.isfinite(x)))
            or not bool(xp.all(xp.isfinite(v)))
            or bool(xp.any(~xp.isfinite(r) | (r <= 0)))
        ):
            raise ValueError("Finite positions/volumes and positive radii required")
        if not np.isfinite([frequency_Hz, density_kg_m3]).all() or frequency_Hz <= 0 or density_kg_m3 <= 0:
            raise ValueError("Positive frequency/density required")
        delta = x[None, :, :] - x[:, None, :]
        d2 = xp.sum(delta * delta, axis=-1)
        diagonal = xp.eye(len(r), dtype=bool)
        if bool(xp.any((~diagonal) & (d2 == 0))):
            raise ValueError("Coincident centers; force undefined")
        if self.enforce_far_field and bool(xp.any((~diagonal) & (d2 < (self.min_ratio * (r[:, None] + r[None, :])) ** 2))):
            raise ValueError("Outside far-field separation; force was not evaluated")
        d2 = xp.where(diagonal, xp.inf, d2)
        if v.ndim == 1:
            v = v[:, None]
        dotv = (2j * np.pi * frequency_Hz) * v * xp.arange(1, v.shape[-1] + 1)
        correlation = 0.5 * xp.real(dotv @ xp.conj(dotv.T))
        return density_kg_m3 / (4 * np.pi) * xp.sum((correlation / d2**1.5)[..., None] * delta, axis=1)


class CuPyPairBackend(ClassicalFarFieldBjerknes):
    def __init__(self, **kwargs):
        super().__init__(backend="cupy", **kwargs)
