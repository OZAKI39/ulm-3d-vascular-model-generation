"""Particle-1: one overdamped sphere. SI only, no inertial state."""
from dataclasses import dataclass
import numpy as np
from .field import FlowSample


def vector3(value, name):
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (3,) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be finite float64 shape (3,)")
    # Immutable bytes backing also prevents callers from resetting WRITEABLE.
    return np.frombuffer(array.tobytes(), dtype=np.float64)


def positive_scalar(value, name):
    if isinstance(value, (bool, np.bool_)) or np.ndim(value) != 0:
        raise ValueError(f"{name} must be a positive finite scalar")
    result = np.float64(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a positive finite scalar")
    return result


@dataclass(frozen=True, slots=True)
class MicrobubbleState:
    particle_id: int
    position_m: np.ndarray
    radius_m: np.float64
    velocity_m_s: np.ndarray
    angular_velocity_s_inv: np.ndarray

    def __post_init__(self):
        if isinstance(self.particle_id, (bool, np.bool_)) or not isinstance(self.particle_id, (int, np.integer)) or self.particle_id < 0:
            raise ValueError("particle_id must be a non-negative integer")
        object.__setattr__(self, "particle_id", int(self.particle_id))
        for name in ["position_m", "velocity_m_s", "angular_velocity_s_inv"]:
            object.__setattr__(self, name, vector3(getattr(self, name), name))
        object.__setattr__(self, "radius_m", positive_scalar(self.radius_m, "radius_m"))


def stokes_drag_force_n(velocity_m_s, background_velocity_m_s, viscosity_pa_s, radius_m):
    """F_hydro=-6*pi*mu*a*(V-u_inf). Diagnostic drag, no force integration."""
    mu = positive_scalar(viscosity_pa_s, "viscosity_pa_s")
    a = positive_scalar(radius_m, "radius_m")
    slip = vector3(velocity_m_s, "velocity_m_s") - vector3(background_velocity_m_s, "background_velocity_m_s")
    return vector3(-6 * np.pi * mu * a * slip, "drag_force_n")


class FieldUnavailableError(ValueError):
    """Sampling outside the closed fluid domain cannot produce an active state."""


def equilibrium_state(state, sample: FlowSample):
    """No other forces: V=u_inf, Omega=one half of curl(u_inf)."""
    if not sample.inside_lumen or sample.tetra_id < 0:
        raise FieldUnavailableError("microbubble center has no frozen fluid field")
    return MicrobubbleState(state.particle_id, state.position_m, state.radius_m,
                            sample.velocity_m_s, 0.5 * sample.vorticity_s_inv)
