"""Frequency-specific opt-in interface; original Segment and kernels unchanged."""
from dataclasses import dataclass,replace
import numpy as np,h5py,json
from mevo_mb.control.scheduled_pzt import Segment,Schedule
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.models import ResponseSection
from scripts.straight_tube_exp1.f10s.common import context as old_context
from scripts.straight_tube_exp1.f10ac_qd.common import *
FREQUENCIES=(2400000.,2500000.,2600000.)
@dataclass(frozen=True)
class FrequencySegment(Segment):
 def __post_init__(self):
  f=self.frequency_Hz
  if f not in FREQUENCIES:raise ValueError('Unregistered discrete frequency; no extrapolation')
  object.__setattr__(self,'frequency_Hz',2500000.)
  super().__post_init__()
  object.__setattr__(self,'frequency_Hz',f)
 def event(self):return dict(super().event(),frequency_Hz=self.frequency_Hz)
def make_schedule(amplitudes,phases,duration=.002,frequency=2500000.):
 return Schedule((FrequencySegment(0,0,duration,'FOUR_CHANNEL_INDEPENDENT',tuple(amplitudes),tuple(phases),tuple(float(a)>0 for a in amplitudes),frequency_Hz=frequency),))
class FixedResponse:
 """Validated linear pressure intervals at one fixed R,f; no frequency interpolation."""
 def __init__(self,path,physics):
  with h5py.File(path) as h:
   self.knots=h['pressure_Pa'][:];self.knot_values=h['harmonics_m3'][:];self.metadata=json.loads(h.attrs['metadata_json'])
  if self.metadata['frequency_Hz']!=physics.frequency_Hz or self.metadata['radius_m']!=physics.radius_m:raise ValueError('Response frequency/radius mismatch')
  if not self.metadata['validated']:raise ValueError('Unvalidated response')
  self.values=np.stack([self.knot_values[:-1],self.knot_values[1:]],1)
  self.valid=np.ones(len(self.knots)-1,np.int32);self.knot_valid=np.ones(len(self.knots),np.int32)
  self.lower=self.knots[:-1];self.upper=self.knots[1:]
 def query(self,p):
  p=np.asarray(p)
  if not np.isfinite(p).all():raise ValueError("Nonfinite pressure is outside the model domain")
  a=np.abs(p);nz=a>0
  if np.any(nz & ((a<self.knots[0])|(a>self.knots[-1]))):raise ValueError('Pressure outside validated table; no clipping')
  out=np.zeros(p.shape+(8,),complex)
  for i in range(8):out[...,i]=np.where(nz,np.interp(a,self.knots,self.knot_values[:,i].real)+1j*np.interp(a,self.knots,self.knot_values[:,i].imag),0)*np.exp(1j*(i+1)*np.angle(p))
  return out

def context(frequency=2500000.,shell=1.,require_admission=True):
 c=old_context()
 if frequency==2500000. and shell==1.:return c
 if frequency not in FREQUENCIES:raise ValueError('Unregistered frequency')
 key=str(round(frequency));folder=AD/key
 check=read(folder/'PHYSICS_GATE.json')
 if require_admission and not check['PASS']:raise ValueError('Frequency physics not admitted')
 c['physics']=replace(c['physics'],frequency_Hz=frequency)
 c['basis']=CachedWaterBasis(folder/'coarse.h5')
 c['response']=FixedResponse(folder/(f'response_shell_{shell:.1f}.h5'),c['physics'])
 return c
