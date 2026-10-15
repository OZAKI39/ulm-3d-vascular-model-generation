"""Frozen-reference LPS force and expanded energy, SI physical nodal force.
Reuses the exact independent algebra validated in native_state_v1/audits.py.
"""
import numpy as np
class LPS:
 def __init__(self,reference):
  self.X=reference['X'];self.b=reference['bonds'];self.r=reference['r0'];self.h=reference['taper'];self.m=reference['m'];self.V=float(reference['V'][0]);self.K=1e5/3;self.G=940.;self.n=len(self.X)
 def eval(self,x,alive,theta_override=None):
  a,j=self.b[alive].T;r=self.r[alive];h=self.h[alive];m=self.m;V=self.V;K=self.K;G=self.G
  y=x[j]-x[a];length=np.linalg.norm(y,axis=1);e=length-r;e[abs(e)<2.2204e-16]=0
  theta=3*V/m*np.bincount(np.r_[a,j],weights=np.tile(h*e,2),minlength=self.n)
  U=.5*V*(K-5*G/3)*np.sum(theta**2)+7.5*G*V**2*np.sum(h/r*e**2*(1/m[a]+1/m[j]))
  tt=theta if theta_override is None else theta_override
  mag=V**2*h*((3*K-5*G)*(tt[a]/m[a]+tt[j]/m[j])+15*G*(1/m[a]+1/m[j])*e/r)
  force=np.zeros_like(x);f=mag[:,None]*y/length[:,None];np.add.at(force,a,f);np.add.at(force,j,-f)
  return theta,force,float(U)
