"""Independent tensile-strain dD/dt model: DEVELOPMENT / UNCALIBRATED.

This module has no fracture coupon scale, fatigue law or cycle-count clock.
"""
import numpy as np


def tensile_strain(stretch):
    value=np.asarray(stretch,dtype=float)
    if not np.isfinite(value).all() or np.any(value<0):
        raise ValueError('Finite nonnegative stretch required')
    return np.maximum(value-1.,0.)


class StreamingErosionDamage:
    def __init__(self,count,config,initial_damage=None):
        self.c=dict(config)
        if self.c.get('D_break')!=1.:raise ValueError('STREAMING_EROSION requires D_break = 1')
        if self.c.get('tension_only') is not True:raise ValueError('Only tension_only=true is implemented')
        if self.c.get('spatial_feedback') is not False:raise ValueError('Spatial feedback is not implemented')
        if self.c.get('demonstration_parameter') is not True:raise ValueError('Erosion parameters must be labeled uncalibrated demonstration inputs')
        for key in ['C_t_s_inv','m_t','epsilon_ref','threshold_ratio']:
            v=self.c[key]
            if not np.isfinite(v) or v<0 or (key in ['m_t','epsilon_ref'] and v==0):
                raise ValueError(f'Invalid erosion parameter {key}')
        self.D=np.zeros(count) if initial_damage is None else np.asarray(initial_damage,dtype=float).copy()
        if self.D.shape!=(count,) or not np.isfinite(self.D).all() or np.any((self.D<0)|(self.D>1)):
            raise ValueError('Initial D must be finite and in [0,1]')
        self.active=self.D<1.;self.g=1.-self.D
        self.last_rate_s_inv=np.zeros(count);self.last_epsilon_plus=np.zeros(count)

    def rate(self,epsilon_plus):
        eps=np.asarray(epsilon_plus,dtype=float)
        if eps.shape!=self.D.shape or not np.isfinite(eps).all() or np.any(eps<0):
            raise ValueError('Finite nonnegative bond tensile strain required')
        with np.errstate(over='raise',invalid='raise'):
            rate=self.c['C_t_s_inv']*np.maximum(eps/self.c['epsilon_ref']-self.c['threshold_ratio'],0.)**self.c['m_t']
        return np.where(self.active,rate,0.)

    def advance(self,epsilon_plus,exposure_step_s):
        if not np.isfinite(exposure_step_s) or exposure_step_s<=0:
            raise ValueError('Positive finite slow exposure step required')
        self.last_epsilon_plus=np.asarray(epsilon_plus,dtype=float).copy()
        self.last_rate_s_inv=self.rate(self.last_epsilon_plus)
        old=self.D.copy();old_active=self.active.copy()
        with np.errstate(over='raise',invalid='raise'):
            requested=exposure_step_s*self.last_rate_s_inv
        self.D=np.minimum(1.,old+requested)
        self.active=self.D<1.;self.g=1.-self.D
        if np.any(self.D<old) or np.any(self.active & ~old_active):raise AssertionError('Damage healed')
        return np.flatnonzero(old_active & ~self.active),old,self.D-old


class SettlingMonitor:
    """End-window strain range on the active graph; kinetic energy also recorded.

    No temporal range is ever passed to the damage law. The damage driver is
    the *last configuration's tensile strain*. Range is only a settling flag.
    """
    def __init__(self,active,absolute_floor=1e-12):
        self.active=np.asarray(active,bool).copy();self.floor=absolute_floor
        if not np.isfinite(absolute_floor) or absolute_floor<=0:raise ValueError('Positive drift floor required')
        self.rows=[];self.low=None;self.high=None

    def add(self,time_s,window_time_s,epsilon_plus,kinetic_energy_J):
        eps=np.asarray(epsilon_plus)[self.active]
        if not np.isfinite(eps).all() or not np.isfinite(kinetic_energy_J):raise FloatingPointError('Nonfinite settling sample')
        if self.low is None:self.low=eps.copy();self.high=eps.copy()
        else:self.low=np.minimum(self.low,eps);self.high=np.maximum(self.high,eps)
        self.rows.append(dict(mechanical_time_s=float(time_s),window_time_s=float(window_time_s),
                              kinetic_energy_J=float(kinetic_energy_J),maximum_epsilon_plus=float(eps.max()) if len(eps) else 0.,
                              mean_epsilon_plus=float(eps.mean()) if len(eps) else 0.))

    def summary(self,tolerance):
        if len(self.rows)<2:raise ValueError('At least two end-window samples required')
        peak=max(r['maximum_epsilon_plus'] for r in self.rows)
        vector=float(np.max(self.high-self.low)/max(peak,self.floor)) if len(self.high) else 0.
        def drift(key):
            a=np.array([r[key] for r in self.rows]);return float(np.ptp(a)/max(float(np.max(np.abs(a))),self.floor))
        max_drift=drift('maximum_epsilon_plus');mean_drift=drift('mean_epsilon_plus')
        worst=max(vector,max_drift,mean_drift)
        return dict(settling_status='DRIVER_SETTLING_CHECK_PASSED' if worst<=tolerance else 'MECHANICAL_STEADY_STATE_NOT_ESTABLISHED',
                    driver_relative_drift=worst,bond_vector_relative_drift=vector,maximum_driver_relative_drift=max_drift,
                    mean_driver_relative_drift=mean_drift,maximum_relative_driver_drift_allowed=tolerance,
                    settling_sample_count=len(self.rows),settling_window_start_mechanical_time_s=self.rows[0]['mechanical_time_s'],
                    settling_window_end_mechanical_time_s=self.rows[-1]['mechanical_time_s'],
                    settling_peak_kinetic_energy_J=max(r['kinetic_energy_J'] for r in self.rows),
                    settling_final_kinetic_energy_J=self.rows[-1]['kinetic_energy_J'])
