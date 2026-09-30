"""One velocity/pressure/gradient field supplies both traction and transport."""
from dataclasses import dataclass
from typing import Protocol
import numpy as np
from numba import njit
from ..streaming import pipe_field, fluid_traction
from ..fragment_fluid import FragmentFluid


@dataclass
class FlowSample:
    velocity: np.ndarray
    pressure: np.ndarray
    gradient: np.ndarray
    du_dt: np.ndarray | None = None


class StreamingFlowProvider(Protocol):
    def sample(self, points, time_s) -> FlowSample: ...
    def traction(self, points, normals, time_s): ...


class ManufacturedStreamingProvider:
    """Curl of a localized Gaussian vector potential plus Poiseuille pipe flow.

    A_y = -U L exp(-r²/(2L²)); u = curl(A).
    Therefore div(u) = 0 analytically. Prescribed p is smooth, but the pair (u,p)
    is a manufactured verification field and NOT a solved microbubble CFD field.
    """
    def __init__(self, config):
        self.c = config
        s = config['streaming']
        self.center = np.asarray(s['bubble_center_m'], float)
        self.L = s['streaming_length_scale_m']
        self.U = s['streaming_velocity_scale_m_s']
        self.mu = config['pipe']['viscosity_Pa_s']
        if self.L <= 0 or self.mu <= 0:
            raise ValueError('Positive length scale and viscosity required')
        self.M = np.array([[0., 0., -1.], [0., 0., 0.], [1., 0., 0.]])

    def localized(self, points, time_s):
        points = np.asarray(points, float)
        d = points-self.center
        E = np.exp(-.5*np.sum(d*d, axis=1)/self.L**2)
        b = d@self.M.T/self.L
        s = self.c['streaming']
        omega = 2*np.pi*self.c['simulation']['representative_frequency_Hz']
        phase = omega*time_s+s['phase_rad']
        factor = s['mean_fraction']+s['oscillatory_amplitude']*np.sin(phase)
        factor_dt = s['oscillatory_amplitude']*omega*np.cos(phase)
        velocity = self.U*E[:, None]*b
        gradient = self.U*E[:, None, None]*(self.M[None, :, :]/self.L-b[:, :, None]*d[:, None, :]/self.L**2)
        return FlowSample(factor*velocity, factor*s.get('pressure_amplitude_Pa', 0.)*E,
                          factor*gradient, factor_dt*velocity)

    def sample(self, points, time_s):
        local = self.localized(points, time_s)
        velocity, pressure, gradient = pipe_field(np.asarray(points, float), self.c['pipe'])
        return FlowSample(velocity+local.velocity, pressure+local.pressure,
                          gradient+local.gradient, local.du_dt)

    def traction(self, points, normals, time_s):
        value = self.sample(points, time_s)
        return fluid_traction(value.gradient, value.pressure, normals, self.mu)


class ConsistentFragmentFluid(FragmentFluid):
    """Reuse legacy surface geometry/load mask, replace only its field source."""
    def __init__(self, config, provider=None):
        self.c = config
        self.flow_provider = provider or ManufacturedStreamingProvider(config)
        self.traction_provider = self.flow_provider

    def velocity(self, positions, time_s):
        return self.flow_provider.sample(positions, time_s).velocity

    def surface_force(self, cloud, positions, integrity, exposed, attached, time_s):
        normals, valid = surface_normals(positions, cloud.pairs, cloud.weight*integrity,
                                          self.c['surface']['normal_moment_tolerance'])
        load = exposed & attached & valid & ~cloud.fixed
        traction = np.zeros_like(positions)
        traction[load] = self.traction_provider.traction(positions[load], normals[load], time_s)
        area = cloud.spacing**2*self.c['surface']['area_factor']
        return traction*area, normals, traction


@njit(cache=True)
def surface_normals(x, pairs, weights, tolerance):
    moment = np.zeros_like(x)
    for k in range(len(pairs)):
        i,j=pairs[k];distance=0.
        for a in range(3):distance+=(x[j,a]-x[i,a])**2
        length=max(np.sqrt(distance),np.finfo(np.float64).tiny)
        for a in range(3):
            value=weights[k]*(x[j,a]-x[i,a])/length
            moment[i,a]+=value;moment[j,a]-=value
    normal=np.zeros_like(x);valid=np.zeros(len(x),np.bool_)
    for i in range(len(x)):
        norm=np.sqrt(np.sum(moment[i]*moment[i]))
        if norm>tolerance:
            valid[i]=True
            normal[i]=-moment[i]/norm
    return normal,valid
