"""Radial response abstraction; no MHz ODE solved in the slow translation loop."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


class BubbleResponseModel(ABC):
    @abstractmethod
    def response(self, R0_m, frequency_Hz, local_complex_pressure_Pa, medium_properties) -> np.ndarray:
        """Return fundamental complex volume response V_hat in m^3."""


@dataclass
class LinearResonantBubbleModel(BubbleResponseModel):
    """Low-pressure architecture surrogate, not a calibrated SonoVue shell model.

    xi''+2*zeta*w0*xi'+w0^2*xi=-p/(rho*R0); V_hat=4*pi*R0^2*xi_hat.
    Positive acoustic pressure compresses a below-resonance bubble.
    """
    resonance_frequency_Hz: float
    damping_ratio: float

    def __post_init__(self):
        if self.resonance_frequency_Hz <= 0 or self.damping_ratio <= 0:
            raise ValueError("Positive resonance and damping required")

    def response(self, R0_m, frequency_Hz, local_complex_pressure_Pa, medium_properties):
        w, w0 = 2*np.pi*frequency_Hz, 2*np.pi*self.resonance_frequency_Hz
        denominator = medium_properties["density_kg_m3"] * (w0*w0 - w*w + 2j*self.damping_ratio*w0*w)
        return -4*np.pi*np.asarray(R0_m)*local_complex_pressure_Pa / denominator


@dataclass(frozen=True)
class MarmottantParameters:
    """All values must come from a separately verified calibration, never defaults."""
    shell_elasticity_N_m: float
    shell_viscosity_kg_s: float
    buckling_radius_m: float
    rupture_radius_m: float
    water_surface_tension_N_m: float
    gas_polytropic_exponent: float
    ambient_pressure_Pa: float


@dataclass
class MarmottantBubbleModel(BubbleResponseModel):
    """Reserved nonlinear response/table interface; explicit, tested refusal to guess."""
    parameters: MarmottantParameters

    def response(self, R0_m, frequency_Hz, local_complex_pressure_Pa, medium_properties):
        raise NotImplementedError("TODO: validated Keller-Miksis + Marmottant or APECSS response table, with calibrated shell parameters and independent ODE tests")
