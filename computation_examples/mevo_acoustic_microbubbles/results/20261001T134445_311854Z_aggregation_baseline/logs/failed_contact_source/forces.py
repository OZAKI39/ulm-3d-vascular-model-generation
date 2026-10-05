"""Harmonic primary and quasistatic far-field secondary Bjerknes forces in N."""
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


def primary_bjerknes(volume_response_m3: np.ndarray, pressure_gradient_Pa_m: np.ndarray) -> np.ndarray:
    return -0.5*np.real(np.asarray(volume_response_m3)[..., None]*np.conj(pressure_gradient_Pa_m))


class PairInteractionBackend(ABC):
    @abstractmethod
    def forces(self, positions_m, radii_m, volume_response_m3, frequency_Hz, density_kg_m3):
        """Net force in N for every bubble; interaction is reciprocal."""


@dataclass
class NumpyPairBackend(PairInteractionBackend):
    minimum_separation_m: float
    collision_regularization_m: float
    cutoff_m: float

    def forces(self, positions_m, radii_m, volume_response_m3, frequency_Hz, density_kg_m3):
        delta = positions_m[None, :, :] - positions_m[:, None, :]
        d = np.linalg.norm(delta, axis=-1)
        valid = (d > 0) & (d <= self.cutoff_m)
        minimum = np.maximum(self.minimum_separation_m, radii_m[:, None] + radii_m[None, :])
        safe = np.maximum(d, minimum+self.collision_regularization_m)
        dotV = (2j*np.pi*frequency_Hz) * volume_response_m3
        correlation = 0.5*np.real(dotV[:, None] * np.conj(dotV[None, :]))
        magnitude = density_kg_m3/(4*np.pi*safe**2)*correlation*valid
        direction = delta / np.maximum(d[..., None], np.finfo(float).tiny)
        return np.sum(magnitude[..., None]*direction, axis=1)
