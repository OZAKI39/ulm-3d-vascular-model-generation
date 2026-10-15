"""Offline variational series-compliance damage of the FULL frozen-reference LPS.
No production imports/mutation. SI: extension m, edge force N, energy J.
The internal elastic extension v minimizes a positive quadratic functional;
it is not an independently prescribed geometric bond extension.
"""
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.linalg import cho_factor, cho_solve

class SofteningLPS:
 def __init__(self,ref,K=1e5/3,G=940.):
  self.X=ref['X'].copy(); self.b=ref['bonds'].copy(); self.r=ref['r0'].copy()
  self.m=ref['m'].copy();self.h=ref['taper'].copy();self.V=float(ref['V'][0]);self.K=K;self.G=G
  assert np.all(ref['V']==self.V) and K>5*G/3
  self.n=len(self.X);self.nb=len(self.b);self.i,self.j=self.b.T
  Lvals=np.r_[3*self.V*self.h/self.m[self.i],3*self.V*self.h/self.m[self.j]]
  self.L=coo_matrix((Lvals,(np.r_[self.i,self.j],np.tile(np.arange(self.nb),2))),shape=(self.n,self.nb)).tocsr()
  self.D=15*G*self.V**2*self.h/self.r*(1/self.m[self.i]+1/self.m[self.j])
  self.A=self.V*(K-5*G/3)
  self.kref=(self.D*self.r+self.A*np.asarray(self.L.T@(self.L@self.r)))/self.r
  assert np.all(self.D>0) and np.all(self.kref>0)
 def geometry(self,x):
  dy=x[self.j]-x[self.i];length=np.linalg.norm(dy,axis=1)
  assert np.all(length>0)
  return length-self.r,dy/length[:,None]
 def assemble(self,t,direction):
  v=t[:,None]*direction;f=np.zeros_like(self.X)
  for k in range(3):f[:,k]=np.bincount(self.i,weights=v[:,k],minlength=self.n)-np.bincount(self.j,weights=v[:,k],minlength=self.n)
  return f
 def baseline(self,x,alive=None):
  if alive is None:alive=np.ones(self.nb)
  alive=np.asarray(alive,float);e,di=self.geometry(x);theta=self.L@(e*alive)
  t=alive*(self.D*e+self.A*np.asarray(self.L.T@theta))
  U=.5*np.dot(self.D*alive,e*e)+.5*self.A*np.dot(theta,theta)
  return dict(U=float(U),force=self.assemble(t,di),theta=theta,t=t,e=e)
 def evaluate(self,x,d=None,mode='CANDIDATE_SOFTENING'):
  if mode=='BASELINE_INTACT_LPS':return self.baseline(x)
  if mode!='CANDIDATE_SOFTENING':raise ValueError(mode)
  d=np.zeros(self.nb) if d is None else np.asarray(d,float)
  if np.any((d<0)|(d>1)):raise ValueError('damage must be in [0,1]')
  e,di=self.geometry(x);q=1-d
  # Stable exact limits: d=0 -> v=e; d=1 -> edge force identically zero.
  den=self.D*d+self.kref*q
  w=self.kref*q/den;h=d/den
  M=(self.L@diags(h)@self.L.T).toarray();M.flat[::self.n+1]+=1/self.A
  p=cho_solve(cho_factor(M,lower=True,check_finite=False),self.L@(w*e),check_finite=False)
  coupling=np.asarray(self.L.T@p);v=w*e-h*coupling
  t=w*(self.D*e+coupling)
  theta=self.L@v
  penalty=.5*np.sum(self.kref*q*d*((self.D*e+coupling)/den)**2)
  U_local=.5*np.dot(self.D,v*v);U_volume=.5*self.A*np.dot(theta,theta)
  U=U_local+U_volume+penalty
  # Y=-dU/dd >= 0, including the d=0 one-sided limit. No abs/clipping.
  Y=.5*self.kref*((self.D*e+coupling)/den)**2
  return dict(U=float(U),force=self.assemble(t,di),theta=theta,t=t,e=e,v=v,Y=Y,
              U_local=float(U_local),U_volume=float(U_volume),U_compliance=float(penalty),
              dual_energy_error_J=float(U-.5*np.dot(e,t)),d=d)
 def update(self,x,kappa,s_onset=.15,s_complete=.24):
  if not s_complete>s_onset>0:raise ValueError('invalid thresholds')
  e,_=self.geometry(x);new=np.maximum(kappa,np.maximum(e/self.r,0.))
  d=np.zeros(self.nb);part=(new>s_onset)&(new<s_complete);full=new>=s_complete
  d[part]=s_complete*(new[part]-s_onset)/(new[part]*(s_complete-s_onset));d[full]=1.
  return new,d
