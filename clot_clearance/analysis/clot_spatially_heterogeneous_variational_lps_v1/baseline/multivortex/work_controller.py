"""Frozen scalar power-budget servo; no outcome-dependent spatial controls."""
from common import *

SPEC=dict(W_TARGET_J=W_TARGET,W_DRIVE_REF_J=27.443461e-9,W_DRAG_REF_J=-1.745264e-9,
    relative_tolerance=.05,nominal_amplitude_cap_Pa=6000.,amplitude_slew_Pa_s=8e6,
    initial_nominal_amplitude_Pa=3200.,tracking_time_s=.00015,
    observation_drag_forecast='K*(1-exp(P_drag/K*(END-OFF))), frozen instantaneous linear-drag decay surrogate',
    schedule='Normalized integral of the common half-cosine envelope squared',
    controller='P_goal = budget*dS/dt + (budget*S-W_net)/tracking_time; common amplitude solves (P_goal-P_drag)/unit_drive_power, bounded/slew limited',
    scope='Both new modes same algorithm, one formal run each. Historical A is budget only, old surface model.',
    safeguards='No node-wise control, no field relocation or shape tuning, no posthoc work scaling; a miss is WORK_MATCH_NOT_ACHIEVED')

class WorkController:
    def __init__(self,field):
        self.field=field
        self.t=np.linspace(0,OFF,28001);q=field.envelope(self.t)**2
        integ=np.r_[0,np.cumsum((q[1:]+q[:-1])*.5*np.diff(self.t))]
        self.S=integ/integ[-1];self.dS=q/integ[-1];self.last={}
    def prepare(self,s,W,h):
        t=s['time'];old=float(self.field.amplitude(t))
        pd=float(np.sum(s['external_drive_force']*s['velocity']))
        pr=float(np.sum(s['external_drag_force']*s['velocity']))
        K=max(float(s['K']),1e-30)
        forecast=K*(1-np.exp(np.clip(pr/K*(END-OFF),-100,0)))
        budget=W_TARGET+forecast
        power=budget*float(np.interp(t,self.t,self.dS))+(budget*float(np.interp(t,self.t,self.S))-W)/SPEC['tracking_time_s']
        coeff=pd/max(old,1e-30)
        wanted=(power-pr)/coeff if coeff>1e-14 and t>1e-4 else 3200.
        if t>=OFF:wanted=0.
        bounded=float(np.clip(wanted,0,SPEC['nominal_amplitude_cap_Pa']))
        next_amp=float(np.clip(bounded,old-SPEC['amplitude_slew_Pa_s']*h,old+SPEC['amplitude_slew_Pa_s']*h))
        self.field.segment=(t,t+h,old,next_amp)
        self.last=dict(nominal_amplitude_Pa=next_amp,nominal_amplitude_requested_Pa=float(wanted),
            controller_saturated=bool(wanted>SPEC['nominal_amplitude_cap_Pa'] or wanted<0),
            controller_slew_limited=bool(abs(next_amp-bounded)>1e-9),forecast_observation_drag_J=forecast,
            controller_budget_J=budget,controller_goal_power_W=power)
        return self.last
