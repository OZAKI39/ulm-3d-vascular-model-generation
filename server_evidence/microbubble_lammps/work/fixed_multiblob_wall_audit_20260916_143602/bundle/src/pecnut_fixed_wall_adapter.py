"""Thin prescribed-motion adapter around unmodified Pecnut operators.

No new hydrodynamic kernel. Target F/T rows of M_inf^{-1}+R_2B,excess
are obtained with standard Schur-complement linear algebra. All fixed wall
U/O/E entries are zero. Full native-matrix equality is audited separately.
"""
from pathlib import Path
import os,sys,time,hashlib
import numpy as np
from scipy.linalg import cho_factor,cho_solve
code=Path(os.environ['PECNUT_CODE']);sys.path.insert(0,str(code))
os.environ.setdefault('NUMBA_CACHE_DIR',str(code.parent/('numba_cache_'+hashlib.sha256((code/'resistance_scalars/scalars_general_resistance_d.npy').read_bytes()).hexdigest()[:16])))
saved_argv=sys.argv[:];sys.argv=['fixed_wall_adapter']
import settings
from numba import config
config.DISABLE_JIT=False
from functions.generate_Minfinity import generate_Minfinity
from functions.generate_R2Bexact import generate_R2Bexact
from functions.generate_grand_resistance_matrix import generate_grand_resistance_matrix
sys.argv=saved_argv

def posdata(x,a):
 x=np.asarray(x,dtype=float);a=np.asarray(a,dtype=float)
 return (a,x,np.zeros((len(a),2,3)),np.empty(0),np.empty((0,3)),np.empty((0,3)))
def index(N,i):return np.r_[np.arange(3*i,3*i+3),np.arange(3*N+3*i,3*N+3*i+3),np.arange(6*N+5*i,6*N+5*i+5)]
def drag(a,mu=1.):return np.array([6*np.pi*mu*a]*3+[8*np.pi*mu*a**3]*3+[20/3*np.pi*mu*a**3]*5)
def dense(x):return x.toarray() if hasattr(x,'toarray') else np.asarray(x)
def matrix(x,a,mu=1.):return generate_Minfinity(posdata(x,a),mu=mu)[0]
def pair_excess(x,a,mu=1.):return dense(generate_R2Bexact(posdata(x,a),cutoff_factor=2,mu=mu)[0])
def native_effective(x,a,target=0,mu=1.):
 t=time.perf_counter();R,_,_,times=generate_grand_resistance_matrix(posdata(x,a),[],regenerate_Minfinity=True,cutoff_factor=2,mu=mu);ids=index(len(a),target)[:6];r=dense(R)[np.ix_(ids,ids)];return r,{'total_seconds':time.perf_counter()-t,'native_times':times}

class FixedWall:
 def __init__(self,x,a):
  self.x=np.asarray(x,float);self.a=np.asarray(a,float);self.N=len(self.a);self.factor=None;self.assembly_seconds=0.;self.factor_seconds=0.;self.peak_operator_bytes=0
  if not self.N:return
  t=time.perf_counter();M=matrix(self.x,self.a);idx=np.concatenate([index(self.N,i) for i in range(self.N)]);self.sqrt=np.sqrt(np.concatenate([drag(aa) for aa in self.a]));self.C=M[np.ix_(idx,idx)]*np.outer(self.sqrt,self.sqrt);self.assembly_seconds=time.perf_counter()-t;del M
  symmetry=np.max(abs(self.C-self.C.T))/np.max(abs(self.C));assert symmetry<=1e-12
  t=time.perf_counter();self.factor=cho_factor(self.C,lower=True,check_finite=True);self.factor_seconds=time.perf_counter()-t;self.peak_operator_bytes=self.C.nbytes+self.factor[0].nbytes
 def query(self,targets,radii,store_forces=False):
  t=time.perf_counter();x=np.atleast_2d(np.asarray(targets,float));a=np.atleast_1d(radii).astype(float);nt=len(a);Qt=np.sqrt(np.concatenate([drag(aa) for aa in a]));ids=np.concatenate([index(nt,i) for i in range(nt)]);A=matrix(x,a)[np.ix_(ids,ids)]*np.outer(Qt,Qt)
  B=np.zeros((11*nt,11*self.N));lub=np.zeros((11*nt,11*nt));ii=index(2,0);jj=index(2,1)
  min_gap=float('inf')
  for j,(w,aw) in enumerate(zip(self.x,self.a)):
   for i,(p,aa) in enumerate(zip(x,a)):
    gap=np.linalg.norm(p-w)-aa-aw;min_gap=min(min_gap,gap)
    if gap<-1e-12:raise ValueError('PHYSICAL_SPHERE_OVERLAP')
    pp=np.vstack([p,w]);r=np.array([aa,aw]);m=matrix(pp,r);sl=slice(11*i,11*i+11);wl=slice(11*j,11*j+11);B[sl,wl]=m[np.ix_(ii,jj)]*np.outer(Qt[sl],self.sqrt[wl])
    if 2*np.linalg.norm(p-w)/(aa+aw)<4:
     delta=pair_excess(pp,r)[np.ix_(ii,ii)];lub[sl,sl]+=delta/np.outer(Qt[sl],Qt[sl])
  if nt>1:
   delta=pair_excess(x,a)[np.ix_(ids,ids)];lub+=delta/np.outer(Qt,Qt)
  assembly=time.perf_counter()-t;t=time.perf_counter()
  Z=cho_solve(self.factor,B.T,check_finite=False) if self.N else np.empty((0,11*nt));S=A-B@Z;Ssym=.5*(S+S.T);ev=np.linalg.eigvalsh(Ssym)
  if ev[0]<=0:raise ValueError('NONPOSITIVE_TARGET_SCHUR')
  Y=np.linalg.solve(S,np.eye(11*nt));total=Y+lub;keep=np.concatenate([np.arange(11*i,11*i+6) for i in range(nt)]);R=total[np.ix_(keep,keep)];solve=time.perf_counter()-t
  residual=float(np.linalg.norm(S@Y-np.eye(11*nt),ord=np.inf));data={'R_scaled':R,'R_SI_units_mu1_a1':R*np.outer(Qt[keep],Qt[keep]),'M_scaled':np.linalg.solve(R,np.eye(6*nt)),'assembly_seconds':assembly,'solve_seconds':solve,'total_seconds':time.perf_counter()-t+assembly,'wall_wall_assembly_seconds':self.assembly_seconds,'wall_factor_seconds':self.factor_seconds,'operator_bytes':self.peak_operator_bytes,'target_schur_condition':float(ev[-1]/ev[0]),'linear_residual':residual,'iterations':0,'solver':'dense Cholesky + exact Schur complement; no pseudoinverse','minimum_target_wall_gap_over_a':None if self.N==0 else min_gap}
  if store_forces:
   # Scaled generalized test velocities are the columns of the identity.
   test=np.zeros((11*nt,6*nt));test[keep,np.arange(6*nt)]=1;data['wall_farfield_generalized_forces_scaled']=-Z@Y@test;data['target_generalized_velocities_scaled']=np.eye(6*nt)
  return data

def lattice(layout,beta,layers,extent,spacing_ratio=1.01,count=None):
 s=2*beta*spacing_ratio;b1=np.array([s,0]);b2=np.array([.5*s,np.sqrt(3)*s/2]) if layout=='HEX' else np.array([0,s]);dz=s*np.sqrt(2/3) if layout=='HEX' else s/np.sqrt(2)
 shift=(b1+b2)/3 if layout=='HEX' else (b1+b2)/2
 m=int(np.ceil(2*extent/s))+3;points=[];layerids=[]
 for layer in range(layers):
  xy=np.array([i*b1+j*b2+(layer%2)*shift for i in range(-m,m+1) for j in range(-m,m+1)]);xy=xy[np.linalg.norm(xy,axis=1)<=extent+1e-12];order=np.lexsort((xy[:,1],xy[:,0],np.linalg.norm(xy,axis=1)));xy=xy[order]
  for p in xy:points.append([*p,-beta-layer*dz]);layerids.append(layer)
 points=np.array(points);layerids=np.array(layerids)
 if count is not None:points=points[:count];layerids=layerids[:count]
 return points,np.full(len(points),beta),layerids,{'b1':b1.tolist(),'b2':b2.tolist(),'spacing':s,'dz':dz,'layer_shift':shift.tolist()}
def phases(layout,beta,spacing_ratio=1.01):
 s=2*beta*spacing_ratio;b1=np.array([s,0]);b2=np.array([.5*s,np.sqrt(3)*s/2]) if layout=='HEX' else np.array([0,s]);special=[('BEAD',np.zeros(2)),('BRIDGE',b1/2),('PORE',(b1+b2)/(3 if layout=='HEX' else 2))];regular=[(f'U{i}_{j}',(i+.5)/4*b1+(j+.5)/4*b2) for i in range(4) for j in range(4)];return special+regular
