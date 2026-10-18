"""Frozen uniform +X tangential slip surrogate, no spatial window or controller."""
from common import *
class ShearField:
 def __init__(self,mode,theta=None,surface_z=0.):
  assert mode=='uniform' and theta is None
  self.mode=mode;self.beta=BETA;self.U=1000./BETA;self.theta=dict(nominal_shear_Pa=1000.,direction=[1.,0.,0.],drive_duration_s=OFF)
 def envelope(self,t):
  t=np.asarray(t);u=np.clip(t/RISE,0.,1.);v=np.clip((t-HOLD)/(OFF-HOLD),0.,1.)
  return np.where(t<0,0.,np.where(t<RISE,3*u*u-2*u*u*u,np.where(t<HOLD,1.,1-3*v*v+2*v*v*v)))
 def spatial(self,x,n):
  raw=np.zeros_like(x);raw[:,0]=1.;return raw-n*n[:,0,None]
 def __call__(self,x,n,t):return self.U*float(self.envelope(t))*self.spatial(x,n)
 def amplitude(self,t):return 1000.
 def applicability(self):
  return dict(label=LABEL,mode='UNIFORM_SURFACE_SHEAR_X',nominal_stationary_shear_cap_Pa=1000.,outer_speed_bound_m_s=self.U,mu_Pa_s=MU,h_eff_m=H_EFF,Re_h_bound=RHO*self.U*H_EFF/MU,acoustic_realizability='UNESTABLISHED',meaning='Prescribed current-surface tangential slip proxy. No solved acoustic streaming or calibrated bubble power. Passive drag remains after drive-off.')
