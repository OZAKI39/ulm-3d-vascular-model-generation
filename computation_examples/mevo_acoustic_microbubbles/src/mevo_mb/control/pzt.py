"""Acoustic source controls, independent of a particular discretization backend."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class PZTElement:
    id: str
    position_m: list[float]
    normal: list[float]
    geometry: str
    size_m: list[float]
    amplitude_Pa: float
    phase_rad: float
    frequency_Hz: float
    enabled: bool = True

    def transfer_function(self, frequency_Hz: float) -> complex:
        """Unity demo response; replace with measured hydrophone calibration."""
        return 1.0 + 0j


class PZTArray:
    def __init__(self, elements: list[PZTElement]):
        self.elements = elements

    @classmethod
    def from_config(cls, config: dict) -> PZTArray:
        frequency = config["acoustics"]["frequency_Hz"]
        return cls([PZTElement(**element, frequency_Hz=frequency) for element in config["pzt"]])

    def _set(self, name: str, values) -> None:
        if len(values) != len(self.elements):
            raise ValueError("One control per element required")
        for element, value in zip(self.elements, values, strict=True):
            setattr(element, name, value)

    def set_amplitudes(self, amplitudes_Pa) -> None:
        if np.any(np.asarray(amplitudes_Pa) < 0):
            raise ValueError("Amplitude must be nonnegative")
        self._set("amplitude_Pa", amplitudes_Pa)

    def set_phases(self, phases_rad) -> None:
        self._set("phase_rad", phases_rad)

    def set_frequency(self, frequency_Hz: float) -> None:
        self._set("frequency_Hz", [frequency_Hz] * len(self.elements))

    def set_enabled(self, enabled) -> None:
        self._set("enabled", enabled)

    def get_control_vector(self) -> np.ndarray:
        return np.array([e.enabled * e.amplitude_Pa * np.exp(1j * e.phase_rad) * e.transfer_function(e.frequency_Hz) for e in self.elements])
