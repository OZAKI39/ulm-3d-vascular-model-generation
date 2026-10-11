from analysis.dual_track_v1.core_a import *
import ctypes,time
from scipy.sparse.linalg import LinearOperator,cg

def risk(state):
 X=state['X'];b=state['bonds'];x=state['x_m'];alive=state['alive'].astype(bool);r=state['r'];s=(np.linalg.norm(delta_y(x,b),axis=1)-r)/r;d=delta_y(X,b);mid=X[b[:,0]]+.5*d;mid[:,1]=(mid[:,1]+PERIOD/2)%PERIOD-PERIOD/2
 tip=float(state.get('tip_m',-.001));front=np.hypot(mid[:,0]-tip,mid[:,2])<=DELTA;xwall=.0025-abs(mid[:,0])<=DELTA;grip=.000475-abs(mid[:,2])<=DELTA
 labels=np.full(len(b),'free_interior',dtype='U40');labels[xwall]='other_X_free_surface';labels[grip]='grip';labels[front]='notch_front_periodic_layer';labels[front&xwall]='notch_finite_X_intersection'
 regions={}
 for name in ['notch_front_periodic_layer','notch_finite_X_intersection','grip','other_X_free_surface','free_interior']:
  ids=np.flatnonzero((labels==name)&alive);regions[name]=dict(bonds=len(ids),max_utilization=float(s[ids].max()/.15) if len(ids) else None)
 ids=np.flatnonzero(alive);order=ids[np.argsort(s[ids])[-12:][::-1]];j=order[0];return dict(max_stretch=float(s[j]),max_utilization=float(s[j]/.15),leading_region=str(labels[j]),regions=regions,top_bonds=[dict(reference_bond_id=int(k+1),atom_ids=(b[k]+1).tolist(),stretch=float(s[k]),utilization=float(s[k]/.15),reference_midpoint_m=mid[k].tolist(),region=str(labels[k])) for k in order],Y_scope='Periodic layer; no free Y surface',natural_failure='Not run; exact frozen-topology utilization criterion >1')

def correction(x,e,alive,free,max_hvp=600):
 lib=ctypes.CDLL(str(ROOT/'analysis/dual_track_v1/hvp_y.so'));hv=lib.hessian_vector;hv.restype=None;hv.argtypes=[ctypes.c_int]*2+[ctypes.c_void_p]*6+[ctypes.c_double]*3+[ctypes.c_void_p]*2;calls=0;t=time.monotonic()
 def H(v):
  nonlocal calls
  calls+=1
  if calls>max_hvp+3:raise RuntimeError('HVP_BUDGET')
  v=np.ascontiguousarray(v);z=np.zeros_like(x);hv(len(x),len(e.b),*[a.ctypes.data for a in [x,e.b,alive,e.r,e.h,e.m]],e.V,K,G,v.ctypes.data,z.ctypes.data);return z
 rng=np.random.default_rng(928);v=rng.normal(size=x.shape);v[~free]=0;v/=np.linalg.norm(v);h=1e-9;fd=-(e.eval(x+h*v,alive,True)[1]-e.eval(x-h*v,alive,True)[1])/(2*h);err=float(np.linalg.norm(fd-H(v))/np.linalg.norm(fd));assert err<1e-6
 modes=np.zeros((2,len(x),3));modes[0,:,0]=1;modes[1,:,1]=1;Q=np.linalg.qr(np.stack([v[free] for v in modes],axis=1))[0];project=lambda z:z-Q@(Q.T@z)
 def matvec(z):
  v=np.zeros_like(x);v[free]=project(z);return project(H(v)[free])+Q@(Q.T@z)
 U,f=e.eval(x,alive,True);rhs=project(f[free]);dq,code=cg(LinearOperator((int(free.sum()),)*2,matvec=matvec,dtype=float),rhs,rtol=1e-7,atol=1e-18,maxiter=max_hvp);assert code==0,('HVP_CG_FAILED',code)
 for alpha in [1.,.5,.25,.125,.0625]:
  trial=x.copy();trial[free]+=alpha*project(dq);ut,ft=e.eval(trial,alive,True)
  if np.linalg.norm(ft[free])<np.linalg.norm(f[free]) and ut<=U+100*np.finfo(float).eps*U:return trial,dict(alpha=alpha,HVP_calls=calls,HVP_FD_relative=err,wall_s=time.monotonic()-t,correction_max_m=float(abs(trial-x).max()))
 raise RuntimeError('NO_RESIDUAL_REDUCING_CORRECTION')
