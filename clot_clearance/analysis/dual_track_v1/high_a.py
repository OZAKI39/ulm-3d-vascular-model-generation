"""One fixed-topology equilibrium; above-threshold states are diagnostic only."""
from analysis.dual_track_v1.core_a import *
from analysis.dual_track_v1.high_helpers import risk,correction
BASE=A
import time
from analysis.dual_track_v1.native_a import Native

from scipy.optimize import minimize
import sys,csv,resource

def run(name,amp,case='A0',warm=None,warm_amp=2e-6):
 out=BASE/'lps'/name;out.mkdir(parents=True,exist_ok=True);assert not (out/'ATTEMPT.json').exists();save(out/'ATTEMPT.json',dict(opening_m=amp,case=case,utc=time.strftime('%FT%TZ',time.gmtime()),equilibrium_attempt=1));start=time.perf_counter();cpu=time.process_time();cfg=json.loads((BASE/'HIGH_NUMERICAL_FREEZE.json').read_text());c=cfg['criteria'];g,e,alive,mask,lo,hi,_,_=model([-.001,-.00095,-.00090][int(case[-1])]);X=g['X'];free=~mask;scale=1e-6;escale=1e-12
 old=np.load(warm if warm else A/case/'equilibrium.npz');x0=X+(old['x_m']-X)*(amp/warm_amp);fixed=X.copy();fixed[lo,2]-=amp/2;fixed[hi,2]+=amp/2;x0[mask]=fixed[mask];calls=0;history=[];states=[];last=None;trial_max=-np.inf
 def unpack(q):
  x=fixed.copy();x[free]=X[free]+q*scale;return np.ascontiguousarray(x)
 def evaluate(q):
  nonlocal calls,last,trial_max
  calls+=1
  if time.perf_counter()-start>900:raise RuntimeError("WALL_BUDGET")
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>8*1024**2:raise RuntimeError("MEMORY_BUDGET")
  if calls>cfg['optimizer']['maxfun']:raise RuntimeError('FUNCTION_BUDGET')
  x=unpack(q);U,f=e.eval(x,alive,True);s=(np.linalg.norm(delta_y(x,e.b),axis=1)-e.r)/e.r;trial_max=max(trial_max,float(s[alive.astype(bool)].max()));last=(q.copy(),x,U,f)
  return U/escale,-f[free]*scale/escale
 def record(q):
  if last is None or not np.array_equal(last[0],q):evaluate(q)
  _,x,U,f=last;top=float(-f[hi,2].sum());bottom=float(-f[lo,2].sum());Q=(top-bottom)/2;states.append((x.copy(),U));states[:]=states[-10:]
  row=dict(iteration=len(history),evaluations=calls,U_J=U,relative_force=float(np.linalg.norm(f[free])/np.linalg.norm(f[mask])),max_force_N=float(abs(f[free]).max()),top_reaction_N=top,bottom_reaction_N=bottom,reaction_imbalance=abs(top+bottom)/abs(Q),energy_span_relative=(max(a[1] for a in states)-min(a[1] for a in states))/U if len(states)==10 else None,position_change_m=float(abs(x-states[0][0]).max()) if len(states)==10 else None)
  history.append(row)
  return row['relative_force']<=c['relative_free_force_L2_over_controlled_force_L2'] and row['max_force_N']<=c['max_free_force_N'] and row['reaction_imbalance']<=c['boundary_reaction_imbalance'] and row['energy_span_relative'] is not None and row['energy_span_relative']<=c['energy_span_last10_relative'] and row['position_change_m']<=c['max_position_change_last10_m']
 q0=(x0[free]-X[free])/scale;evaluate(q0);record(q0)
 def callback(intermediate_result):
  if record(intermediate_result.x):raise StopIteration
 result=minimize(evaluate,q0,jac=True,method='L-BFGS-B',callback=callback,options=cfg['optimizer']);x=unpack(result.x);U,f=e.eval(x,alive,True)
 if not np.array_equal(states[-1][0],x):record(result.x)
 corrections=[]
 row=history[-1]
 # Same verified local residual solver; limited one correction per wide candidate.
 if not (row['relative_force']<=1e-6 and row['max_force_N']<=1e-13):
  x,cr=correction(x,e,alive,free);corrections.append(cr);q=(x[free]-X[free])/scale;evaluate(q);record(q);U,f=e.eval(x,alive,True);row=history[-1]
 with (out/'history.csv').open('w') as fp:w=csv.DictWriter(fp,fieldnames=list(row));w.writeheader();w.writerows(history)
 n=Native(out/'native',g,alive);fn,diag,theta=n.evaluate(x);n.close();Uc,fc,cp=contact(x,g,alive);err=float(np.linalg.norm(fn-(f+fc))/np.linalg.norm(fn));nrel=float(np.linalg.norm(fn[free])/np.linalg.norm(fn[mask]));Q=(row['top_reaction_N']-row['bottom_reaction_N'])/2;cr=float(np.linalg.norm(fc)/abs(Q));eps=2*max(max(a[1] for a in states)-min(a[1] for a in states),np.linalg.norm(f[free])*np.linalg.norm((x-X)[free]),100*np.finfo(float).eps*U)
 gates=dict(force=row['relative_force']<=c['relative_free_force_L2_over_controlled_force_L2'] and row['max_force_N']<=c['max_free_force_N'],energy_plateau=row['energy_span_relative'] is not None and row['energy_span_relative']<=c['energy_span_last10_relative'],coordinate_plateau=row['position_change_m'] is not None and row['position_change_m']<=c['max_position_change_last10_m'],reaction_balance=row['reaction_imbalance']<=c['boundary_reaction_imbalance'],native_force=nrel<=c['relative_free_force_L2_over_controlled_force_L2'] and abs(fn[free]).max()<=c['max_free_force_N'],force_energy=err<=c['native_energy_force_relative'],contact=Uc/U<=c['max_contact_energy_over_U'] and cr<=c['max_contact_force_over_reaction'],topology=diag['natural_broken']==0)
 # Central FD at eight high-force/free components, step1pm; no extra equilibrium solve.
 testids=np.unique(np.r_[np.argsort(abs(f.ravel()))[-4:],np.flatnonzero(free.ravel())[np.argsort(abs(f[free]))[-4:]]]);fd=[];h=1e-12
 for j in testids:
  xp=x.copy();xm=x.copy();xp.ravel()[j]+=h;xm.ravel()[j]-=h;fdv=-(e.eval(xp,alive,False)-e.eval(xm,alive,False))/(2*h);fd.append(dict(dof=int(j),finite_difference_N=fdv,gradient_N=float(f.ravel()[j]),absolute_error_N=abs(fdv-f.ravel()[j])))
 fdscaled=max(a['absolute_error_N'] for a in fd)/abs(fn).max();gates['finite_difference']=fdscaled<=c['finite_difference_force_scaled_absolute']
 state=dict(tip_m=[-.001,-.00095,-.00090][int(case[-1])],X=X,x_m=x,displacement_m=x-X,alive=alive,bonds=e.b,m=e.m,r=e.r,h=e.h,mask=mask,low=lo,high=hi,independent_force_N=f,native_force_N=fn,native_theta=theta,contact_force_N=fc,contact_pairs=cp,last_positions=np.stack([a[0] for a in states]),last_U_J=np.array([a[1] for a in states]));np.savez_compressed(out/'equilibrium.npz',**state);r=risk(state)
 audit=dict(corrections=corrections,name=name,case=case,opening_m=amp,accepted=bool(all(gates.values())),gates={k:bool(v) for k,v in gates.items()},U_J=U,empirical_uncertainty_J=float(eps),last=row,native=diag,native_force_energy_relative=err,native_free_force_relative=nrel,max_native_free_force_N=float(abs(fn[free]).max()),contact_U_J=Uc,contact_force_over_Q=cr,contact_pairs=len(cp),risk=r,branch=('UNACCEPTED_STATIC_CANDIDATE' if not all(gates.values()) else 'SUBCRITICAL_FROZEN_TOPOLOGY_ELASTIC_EQUILIBRIUM' if r['max_utilization']<1 else 'FROZEN_TOPOLOGY_HYPOTHETICAL_ELASTIC_BRANCH'),native_natural_breaks='DISABLED_ONLY_FOR_DIAGNOSTIC; criterion read from unchanged source and evaluated on all alive bonds',max_stretch_all_trials=trial_max,finite_difference=fd,finite_difference_scaled_error=float(fdscaled),optimizer_message=str(result.message),iterations=result.nit,evaluations=calls,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,threads=1,MPI_ranks=1,GPU=False)
 save(out/'RESULT.json',audit);print((out/'RESULT.json').read_text(),flush=True);return audit
