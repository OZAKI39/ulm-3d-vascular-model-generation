"""Independent single-bubble equations integrated as a vectorized ensemble.
No F2; SciPy RK45 dense output, each member retains its own ID. A shared terminal
field buffer censors every remaining member at the first exit, reported explicitly.
"""
import time,json
import numpy as np,h5py
from scipy.integrate import solve_ivp
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis,cubic_query
from scripts.straight_tube_exp1.f10af.field import response_context
from .common import *
from .capture import capture

def integrate(label,f,phase,x,rtol=1e-8,atol=1e-13,max_step=.001,cache='fine',duration=.25):
 begin=time.perf_counter();c=response_context(f);basis=CachedWaterBasis(BASE/'fields'/str(int(f))/(cache+'.h5'));a=np.array([0,0,1000,1000]);p=np.array([0,0,0,np.radians(phase)]);u=a*np.exp(1j*p);pressure=np.einsum('i,ixyz->xyz',u,basis.pressures);zeta=6*np.pi*.001*R;Pmin=np.inf;Pmax=0.
 def rhs(t,y):
  nonlocal Pmin,Pmax
  q=y.reshape(-1,3);P,G=cubic_query(pressure,basis.origin,basis.spacing,q);V=c['response'].query(P);Pmin=min(Pmin,float(abs(P).min()));Pmax=max(Pmax,float(abs(P).max()));F=-.5*np.real(V[:,0,None]*G.conj());return (F/zeta).ravel()
 def event(t,y):
  q=y.reshape(-1,3);return min(float((q-(basis.origin+2*basis.spacing)).min()),float(((basis.origin+(np.array(basis.shape)-3)*basis.spacing)-q).min()))
 event.terminal=True;event.direction=-1
 sol=solve_ivp(rhs,[0,duration],x.ravel(),rtol=rtol,atol=atol,max_step=max_step,dense_output=True,events=event);ts=np.unique(np.r_[np.arange(0,sol.t[-1],.0005),sol.t[-1]]);xx=sol.sol(ts).T.reshape(len(ts),len(x),3);vel=np.array([rhs(t,z.ravel()).reshape(-1,3) for t,z in zip(ts,xx)])
 roots=np.array([r['position_m'][0] for rec in read(BASE/'STATIC_RESULTS.json')['rows'] if rec['frequency_Hz']==f and rec['phase_deg']==phase for r in rec['roots']]);rows={str(int(th*1e6)):capture(ts,xx,roots,th) for th in read(CONFIG)['capture']['thresholds_m']};path=BASE/'single'/(label+'.h5');path.parent.mkdir(exist_ok=True)
 with h5py.File(path,'x') as h:
  for k,v in dict(time_s=ts,positions_m=xx,velocity_m_s=vel,primary_force_N=vel*zeta,particle_id=np.arange(len(x)),accepted_time_s=sol.t,accepted_positions_m=sol.y.T.reshape(-1,len(x),3)).items():h[k]=v
  h.attrs.update(scope='SINGLE_BUBBLE_F1_STOKES_DIAGNOSTIC; mutually independent particles, no F2',sampling='SciPy RK45 dense output; raw accepted nodes saved separately',identity_json=json.dumps(identity()))
 r=dict(id=label,frequency=f,phase_deg=phase,N_independent=len(x),completed_s=float(sol.t[-1]),status='PASS' if sol.success and sol.status==0 else 'CENSORED',reason=sol.message,shared_boundary_censor=bool(sol.status==1),pressure_range_Pa=[Pmin,Pmax],capture=rows,nfev=sol.nfev,accepted_steps=len(sol.t)-1,wall_s=time.perf_counter()-begin,trajectory_path=str(path.relative_to(BASE)),sha256=sha(path),rtol=rtol,atol_m=atol,max_step_s=max_step,cache=cache,identity=identity());write(BASE/'single'/(label+'.json'),r);return r

def main():
 begin=time.perf_counter();verify();static=read(BASE/'STATIC_RESULTS.json');ref=next(v for v in static['rows'] if v['frequency_Hz']==2500000 and v['phase_deg']==0);roots=np.array([r['position_m'][0] for r in ref['roots']]);x=[]
 for root in roots:
  for d in [15e-6,45e-6,100e-6]:
   for sign in [-1,1]:
    for y in [-350e-6,0,350e-6]:
     for z in [-60e-6,0,60e-6]:x.append([root+sign*d,y,z])
 x=np.array(x);write(BASE/'SINGLE_PLAN.json',dict(positions_m=x,frequency=2500000,phase_deg=0,duration_s=.25,thresholds=read(CONFIG)['capture'],identity=identity()))
 out=[integrate('atlas',2500000,0,x),integrate('atlas_refined',2500000,0,x,rtol=1e-9,atol=1e-14,max_step=.0005)]
 # Identical absolute initial positions, selected before frequency/phase dynamics.
 xx=np.array([[45e-6,y,z] for y in [-200e-6,0,200e-6] for z in [-60e-6,0,60e-6]])
 write(BASE/'SINGLE_FREQUENCY_PHASE_PLAN.json',dict(positions_m=xx,frequencies=[2400000,2500000,2600000],phases_deg=[0,90,180,270],identity=identity()))
 for f in [2400000,2500000,2600000]:
  for phase in [0,90,180,270]:out.append(integrate(f'freq_{f}_phase_{phase}',f,phase,xx))
 with h5py.File(BASE/'single/atlas.h5') as a,h5py.File(BASE/'single/atlas_refined.h5') as b:error=float(np.linalg.norm(a['positions_m'][:]-b['positions_m'][:],axis=2).max())
 write(BASE/'SINGLE_RESULTS.json',dict(rows=out,max_refinement_difference_m=error,PASS=error<50e-9,wall_s=time.perf_counter()-begin,identity=identity()))
 print('SINGLE_COMPLETE',time.perf_counter()-begin,error,flush=True)
if __name__=='__main__':main()
