"""Exactly one bounded static equilibrium for the requested frozen case."""
from analysis.lps_pvce_v1.core import *
import sys,time,csv,resource,traceback
from scipy.optimize import minimize
from scipy.spatial import cKDTree

def contact(x,g,alive):
 b=g['bonds'];n=len(x);p=cKDTree(x).query_pairs(DX,output_type='ndarray');p=p[~np.isin(p[:,0]*n+p[:,1],b[alive.astype(bool),0]*n+b[alive.astype(bool),1])]
 d=x[p[:,0]]-x[p[:,1]];r=np.linalg.norm(d,axis=1);gap=DX-r;c=15*(18*K/(math.pi*DELTA**4))*DX**6/DELTA
 f=np.zeros_like(x)
 if len(p):
  ff=c*gap[:,None]*d/r[:,None];np.add.at(f,p[:,0],ff);np.add.at(f,p[:,1],-ff)
 return float(.5*c*np.sum(gap**2)),f,p

def run(case):
 start=time.perf_counter();cpu=time.process_time();out=BASE/case;freeze=json.loads((BASE/'PVCE_DESIGN_FREEZE.json').read_text());cfg=json.loads((out/'config.json').read_text());assert not (out/'EXECUTION.json').exists()
 g,e,alive,mask,low,high,_,_=model(cfg['tip_m']);X=g['X'];free=~mask;scale=1e-6;Escale=1e-12;amp=freeze['amplitude_m'];H=float(X[high,2].mean()-X[low,2].mean());fixed=X.copy();fixed[low,2]-=amp/2;fixed[high,2]+=amp/2
 x0=X.copy();x0[:,2]*=1+amp/H
 n=Native(out/'native',g,alive);fn,diag,tn=n.evaluate(x0);u0,f0=e.eval(np.ascontiguousarray(x0),alive,True);err=np.linalg.norm(fn-f0)/np.linalg.norm(fn);assert err<1e-9 and diag['contact_J']==0
 save(out/'START_NATIVE_CHECK.json',dict(force_relative=float(err),U_J=u0,**diag))
 history=[];states=[];calls=0;maxstretch=-np.inf;last=None
 def unpack(q):
  x=fixed.copy();x[free]=X[free]+q*scale;return np.ascontiguousarray(x)
 def evaluate(q):
  nonlocal calls,maxstretch,last
  calls+=1
  if calls>freeze['budget']['per_state_max_function_evaluations']:raise RuntimeError('FUNCTION_BUDGET')
  x=unpack(q);U,f=e.eval(x,alive,True);b=e.b;st=(np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)-e.r)/e.r;maxstretch=max(maxstretch,float(st[alive.astype(bool)].max()))
  if maxstretch>=.15:raise RuntimeError('UNPRESCRIBED_FAILURE_CRITERION_REACHED')
  last=(q.copy(),x,U,f)
  return U/Escale,-f[free]*scale/Escale
 def record(q):
  if last is None or not np.array_equal(last[0],q):evaluate(q)
  _,x,U,f=last;den=float(np.linalg.norm(f[mask]));rel=float(np.linalg.norm(f[free])/den) if den>0 else None;top=float(-f[high,2].sum());bottom=float(-f[low,2].sum());Q=(top-bottom)/2
  row=dict(iteration=len(history),evaluations=calls,U_J=U,relative_force=rel,max_force_N=float(abs(f[free]).max()),top_reaction_N=top,bottom_reaction_N=bottom,reaction_imbalance=abs(top+bottom)/abs(Q) if Q!=0 else None,energy_span_relative=None,position_change_m=None,kinetic_J=None)
  states.append((x.copy(),U));states[:]=states[-10:]
  if len(states)==10:
   row['energy_span_relative']=(max(t[1] for t in states)-min(t[1] for t in states))/U;row['position_change_m']=float(np.max(np.abs(x-states[0][0])))
  history.append(row)
  with (out/'history.csv').open('w') as fp:w=csv.DictWriter(fp,fieldnames=list(row));w.writeheader();w.writerows(history)
  c=freeze['criteria'];passed=rel is not None and rel<=c['relative_free_force_L2_over_controlled_force_L2'] and row['max_force_N']<=c['max_free_force_N'] and row['energy_span_relative'] is not None and row['energy_span_relative']<=c['energy_span_last10_relative'] and row['position_change_m']<=c['max_position_change_last10_m'] and row['reaction_imbalance']<=c['boundary_reaction_imbalance']
  if len(history)%100==0:print(case,row,flush=True)
  return passed
 q0=(x0[free]-X[free])/scale;evaluate(q0);record(q0)
 def callback(intermediate_result):
  if record(intermediate_result.x):raise StopIteration
 result=minimize(evaluate,q0,method='L-BFGS-B',jac=True,callback=callback,options=freeze['optimizer']);x=unpack(result.x);U,f=e.eval(x,alive,True)
 if last is None or not np.array_equal(last[0],result.x):evaluate(result.x)
 if not np.array_equal(states[-1][0],x):record(result.x)
 fn,diag,theta=n.evaluate(x);native_error=float(np.linalg.norm(fn-f)/np.linalg.norm(fn));Uc,fc,cp=contact(x,g,alive)
 native_total_error=float(np.linalg.norm(fn-(f+fc))/np.linalg.norm(fn));row=history[-1];epsU=2*max(max(t[1] for t in states)-min(t[1] for t in states),float(np.linalg.norm(f[free])*np.linalg.norm((x-X)[free])),100*np.finfo(float).eps*U)
 Q=(row['top_reaction_N']-row['bottom_reaction_N'])/2;contact_ratio=float(np.linalg.norm(fc)/abs(Q));c=freeze['criteria'];gates=dict(relative_force=row['relative_force']<=c['relative_free_force_L2_over_controlled_force_L2'],max_force=row['max_force_N']<=c['max_free_force_N'],energy_plateau=row['energy_span_relative'] is not None and row['energy_span_relative']<=c['energy_span_last10_relative'],position_stable=row['position_change_m'] is not None and row['position_change_m']<=c['max_position_change_last10_m'],reaction_balance=row['reaction_imbalance']<=c['boundary_reaction_imbalance'],natural_breaks=diag['natural_broken']==0 and maxstretch<.15,contact=Uc/U<=c['max_contact_energy_over_U'] and contact_ratio<=c['max_contact_force_over_reaction'],force_energy=native_total_error<=c['native_energy_force_relative'])
 # Actual native free forces also satisfy same gates; no denominator floor.
 native_den=float(np.linalg.norm(fn[mask]));native_rel=float(np.linalg.norm(fn[free])/native_den);gates['native_free_force']=native_rel<=c['relative_free_force_L2_over_controlled_force_L2'] and float(abs(fn[free]).max())<=c['max_free_force_N']
 np.savez_compressed(out/'equilibrium.npz',X=X,x_m=x,displacement_m=x-X,alive=alive,bonds=e.b,m=e.m,r=e.r,h=e.h,mask=mask,low=low,high=high,independent_force_N=f,native_force_N=fn,native_theta=theta,contact_force_N=fc,contact_pairs=cp,last_positions=np.stack([v[0] for v in states]),last_U_J=np.array([v[1] for v in states]))
 audit=dict(case=case,accepted=all(gates.values()),gates=gates,last=row,U_J=U,empirical_uncertainty_J=epsU,uncertainty_scope=freeze['uncertainty'],native=diag,native_total_force_relative=native_total_error,native_bulk_force_relative=native_error,native_free_force_relative=native_rel,max_native_free_force_N=float(abs(fn[free]).max()),contact_independent_J=Uc,contact_native_energy_difference_J=Uc-diag['contact_J'],contact_force_over_reaction=contact_ratio,contact_pairs=len(cp),max_stretch_all_optimizer_evaluations=maxstretch,kinetic='NA_STATIC_MINIMIZATION',optimizer_message=str(result.message),optimizer_status=int(result.status),iterations=int(result.nit),function_evaluations=calls,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,mpi_ranks=1,OMP_threads=1,GPU=False)
 save(out/'EQUILIBRIUM.json',audit);save(out/'EXECUTION.json',dict(status='ACCEPTED' if audit['accepted'] else 'STOPPED_UNCONVERGED',wall_seconds=audit['wall_seconds'],formal_equilibrium_attempts=1));n.close();print(json.dumps(audit,indent=2),flush=True)
if __name__=='__main__':
 try:run(sys.argv[1])
 except BaseException:
  (BASE/sys.argv[1]/'ERROR.txt').write_text(traceback.format_exc());raise
