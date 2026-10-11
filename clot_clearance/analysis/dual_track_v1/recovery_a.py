from analysis.dual_track_v1.common import *
from analysis.lps_pvce_v1.core import model
from analysis.lps_pvce_v1.run import contact
from analysis.lps_pvce_v2.native import Native
import ctypes,time,resource,subprocess,traceback,csv
from scipy.sparse.linalg import LinearOperator,cg

def main():
 assert json.loads((A/'PVCE_ENERGY_RELEASE_RESULTS.json').read_text())['status']=='PASS'
 out=A/'old_high_recovery';out.mkdir();save(out/'ATTEMPT.json',dict(formal_equilibrium=1,kind='Single targeted local Hessian-vector residual correction',started_unix=time.time()));start=time.monotonic();cpu=time.process_time();src=ROOT/'analysis/dual_track_v1/hvp.cpp';r=subprocess.run(['g++','-O3','-shared','-fPIC',str(src),'-o',str(src.with_suffix('.so'))],capture_output=True,text=True);assert r.returncode==0
 old=V2/'lps/A0_search_01/equilibrium.npz';state=np.load(old);g,e,alive,mask,lo,hi,_,_=model(-100e-6);x=state['x_m'].copy();X=g['X'];free=~mask;criteria=json.loads((A/'DESIGN_FREEZE.json').read_text())['criteria'];hist=list(zip(state['last_positions'].copy(),state['last_U_J'].copy()));hist=hist[-10:];records=[];hvp_calls=0
 lib=ctypes.CDLL(str(src.with_suffix('.so')));hv=lib.hessian_vector;hv.restype=None;hv.argtypes=[ctypes.c_int]*2+[ctypes.c_void_p]*6+[ctypes.c_double]*3+[ctypes.c_void_p]*2
 def H(v,x):
  nonlocal hvp_calls
  hvp_calls+=1
  if hvp_calls>1500 or time.monotonic()-start>300:raise RuntimeError('RECOVERY_BUDGET')
  v=np.ascontiguousarray(v);z=np.zeros_like(x);hv(len(x),len(e.b),*[a.ctypes.data for a in [x,e.b,alive,e.r,e.h,e.m]],e.V,K,G,v.ctypes.data,z.ctypes.data);return z
 # One small directional finite-difference certification before local correction.
 rng=np.random.default_rng(391);v=rng.normal(size=x.shape);v[mask]=0;v/=np.linalg.norm(v);step=1e-9;fd=-(e.eval(x+step*v,alive,True)[1]-e.eval(x-step*v,alive,True)[1])/(2*step);herr=float(np.linalg.norm(fd-H(v,x))/np.linalg.norm(fd));assert herr<1e-6
 def measure(x,U,f):
  top=-f[hi,2].sum();bottom=-f[lo,2].sum();q=(top-bottom)/2;return dict(U_J=U,relative_force=float(np.linalg.norm(f[free])/np.linalg.norm(f[mask])),max_force_N=float(abs(f[free]).max()),reaction_imbalance=float(abs(top+bottom)/abs(q)),top_reaction_N=float(top),bottom_reaction_N=float(bottom),energy_span_relative=float((max(v[1] for v in hist)-min(v[1] for v in hist))/U),position_change_m=float(abs(x-hist[0][0]).max()))
 def gates(row):return row['relative_force']<=criteria['relative_free_force'] and row['max_force_N']<=criteria['max_free_force_N'] and row['reaction_imbalance']<=criteria['reaction_imbalance'] and row['energy_span_relative']<=criteria['energy_plateau_relative'] and row['position_change_m']<=criteria['position_plateau_m']
 U,f=e.eval(x,alive,True);initial=measure(x,U,f)
 for iteration in range(10):
  # Remove only exact rigid translations x/y and rotation about z; no grip changes.
  modes=np.zeros((3,len(x),3));modes[0,:,0]=1;modes[1,:,1]=1;modes[2,:,0]=-x[:,1];modes[2,:,1]=x[:,0];Q=np.linalg.qr(np.stack([v[free] for v in modes],axis=1))[0]
  project=lambda v:v-Q@(Q.T@v)
  def matvec(z):
   zz=np.zeros_like(x);zz[free]=project(z);return project(H(zz,x)[free])+Q@(Q.T@z)
  operator=LinearOperator((int(free.sum()),)*2,matvec=matvec,dtype=float);rhs=project(f[free]);dq,code=cg(operator,rhs,rtol=1e-7,atol=1e-18,maxiter=min(600,1490-hvp_calls));assert code==0,('Krylov did not converge',code)
  before=float(np.linalg.norm(f[free]));accepted=False
  for alpha in [1.,.5,.25,.125,.0625]:
   trial=x.copy();trial[free]+=alpha*project(dq);ut,ft=e.eval(trial,alive,True)
   if np.linalg.norm(ft[free])<before and ut<=U+100*np.finfo(float).eps*U:accepted=True;break
  assert accepted,'No residual-reducing local Newton correction'
  x,U,f=trial,ut,ft;hist.append((x.copy(),U));hist=hist[-10:];row=measure(x,U,f);row.update(correction=iteration+1,alpha=alpha,HVP_calls=hvp_calls,step_max_m=float(abs(alpha*dq).max()));records.append(row)
  if gates(row):break
 native=Native(out/'native',g,alive);fn,diag,theta=native.evaluate(x);native.close();Uc,fc,cp=contact(x,g,alive);nrel=float(np.linalg.norm(fn[free])/np.linalg.norm(fn[mask]));err=float(np.linalg.norm(fn-f-fc)/np.linalg.norm(fn));q=(row['top_reaction_N']-row['bottom_reaction_N'])/2;checks=dict(original_equilibrium=gates(row),native_free=nrel<=criteria['relative_free_force'] and abs(fn[free]).max()<=criteria['max_free_force_N'],native_energy=err<=criteria['native_force_relative'],contact=Uc/U<=criteria['contact_energy_fraction'] and np.linalg.norm(fc)/abs(q)<=criteria['contact_force_fraction'],topology=diag['natural_broken']==0,material_load_unchanged=np.array_equal(x[mask],state['x_m'][mask]))
 np.savez_compressed(out/'equilibrium.npz',**{k:state[k] for k in state.files if k not in ['x_m','displacement_m','independent_force_N','native_force_N','native_theta','last_positions','last_U_J']},x_m=x,displacement_m=x-X,independent_force_N=f,native_force_N=fn,native_theta=theta,last_positions=np.stack([v[0] for v in hist]),last_U_J=np.array([v[1] for v in hist]))
 audit=dict(status='PASS' if all(checks.values()) else 'FAIL',formal_attempts=1,source_state=str(old),source_sha=sha(old),opening_m=73.56055769588697e-6,initial=initial,corrections=records,checks=checks,HVP_FD_relative=herr,HVP_calls=hvp_calls,native_free_force_relative=nrel,max_native_free_force_N=float(abs(fn[free]).max()),native_force_energy_relative=err,native=diag,U_J=U,wall_s=time.monotonic()-start,cpu_s=time.process_time()-cpu,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,interpretation='Restores equilibrium of prior subcritical candidate only; does not prove notch-first onset or material Gc')
 save(A/'STATIC_SOLVER_RECOVERY.json',audit);print((A/'STATIC_SOLVER_RECOVERY.json').read_text())
if __name__=='__main__':
 try:main()
 except BaseException as e:save(A/'STATIC_SOLVER_RECOVERY.json',dict(status='FAIL',formal_attempts=int((A/'old_high_recovery/ATTEMPT.json').exists()),error=repr(e),traceback=traceback.format_exc()));raise
