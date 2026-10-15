"""Offline variational series-compliance damage of the FULL frozen-reference LPS.
No production imports/mutation. SI: extension m, edge force N, energy J.
The internal elastic extension v minimizes a positive quadratic functional;
it is not an independently prescribed geometric bond extension.
"""
import numpy as np
from scipy.sparse import coo_matrix

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
  # Independent endpoint-star elastic states are essential: sharing the free
  # internal extension across a failed edge would couple disconnected bodies.
  self.center=np.r_[self.i,self.j];self.edge=np.tile(np.arange(self.nb),2)
  self.l=Lvals
  self.Ddir=15*G*self.V**2*np.tile(self.h/self.r,2)/self.m[self.center]
  self.krefdir=9*K*self.V**2*np.tile(self.h/self.r,2)/self.m[self.center]
  self.kref=self.krefdir[:self.nb]+self.krefdir[self.nb:]
  assert np.all(self.Ddir>0) and np.all(self.krefdir>0)
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
  e,di=self.geometry(x);dd=d[self.edge];q=1-dd;ee=e[self.edge]
  den=self.Ddir*dd+self.krefdir*q
  w=self.krefdir*q/den;h=dd/den
  # Exact rank-one local solve for each original endpoint star (no global solve).
  numerator=np.bincount(self.center,weights=self.l*w*ee,minlength=self.n)
  denominator=1/self.A+np.bincount(self.center,weights=self.l*self.l*h,minlength=self.n)
  p=numerator/denominator
  coupling=self.l*p[self.center];v=w*ee-h*coupling
  tdir=w*(self.Ddir*ee+coupling);t=tdir[:self.nb]+tdir[self.nb:]
  theta=np.bincount(self.center,weights=self.l*v,minlength=self.n)
  penalty=.5*np.sum(self.krefdir*q*dd*((self.Ddir*ee+coupling)/den)**2)
  U_local=.5*np.dot(self.Ddir,v*v);U_volume=.5*self.A*np.dot(theta,theta)
  U=U_local+U_volume+penalty
  # Shared irreversible edge damage; Y is the SUM of endpoint-star derivatives.
  Ydir=.5*self.krefdir*((self.Ddir*ee+coupling)/den)**2
  Y=Ydir[:self.nb]+Ydir[self.nb:]
  return dict(U=float(U),force=self.assemble(t,di),theta=theta,t=t,e=e,v=v,Y=Y,
              U_local=float(U_local),U_volume=float(U_volume),U_compliance=float(penalty),
              dual_energy_error_J=float(U-.5*np.dot(e,t)),d=d)
 def update(self,x,kappa,s_onset=.15,s_complete=.24):
  if not s_complete>s_onset>0:raise ValueError('invalid thresholds')
  e,_=self.geometry(x);new=np.maximum(kappa,np.maximum(e/self.r,0.))
  d=np.zeros(self.nb);part=(new>s_onset)&(new<s_complete);full=new>=s_complete
  d[part]=s_complete*(new[part]-s_onset)/(new[part]*(s_complete-s_onset));d[full]=1.
  return new,d
