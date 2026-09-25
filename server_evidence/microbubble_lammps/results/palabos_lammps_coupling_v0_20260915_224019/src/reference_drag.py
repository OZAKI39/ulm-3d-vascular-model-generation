"""Analytical/Python reference independent of the C++ coupling adapter."""
import numpy as np
from scipy.integrate import solve_ivp
def drag(mu,diameter,fluid,velocity):return 3*np.pi*mu*np.asarray(diameter)[...,None]*(np.asarray(fluid)-np.asarray(velocity))
def tau(d,rho=1000.,mu=.001):return rho*np.asarray(d)**2/(18*mu)
def uniform_trajectory(times,x0,v0,fluid,diameter):
 t=np.asarray(times)[:,None];tp=tau(diameter);u=np.asarray(fluid);v=np.asarray(v0);x=np.asarray(x0)
 decay=np.exp(-t/tp);return x+u*t+(v-u)*tp*(-np.expm1(-t/tp)),u+(v-u)*decay

def linear_trajectory(times,x0,v0,diameter,intercept,matrix):
 tp=float(tau(diameter));A=np.asarray(matrix);b=np.asarray(intercept)
 def rhs(t,y):return np.r_[y[3:],(b+A@y[:3]-y[3:])/tp]
 sol=solve_ivp(rhs,(0,float(times[-1])),np.r_[x0,v0],method='DOP853',t_eval=times,rtol=1e-12,atol=np.r_[[1e-20]*3,[1e-15]*3])
 if not sol.success:raise RuntimeError(sol.message)
 return sol.y[:3].T,sol.y[3:].T,{'method':'DOP853','rtol':1e-12,'atol_position_m':1e-20,'atol_velocity_m_s':1e-15,'nfev':sol.nfev}
