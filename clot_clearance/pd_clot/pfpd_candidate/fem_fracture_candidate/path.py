"""One sequential displacement-controlled path with accepted-state irreversibility."""
import time,traceback,itertools,csv
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from .specimen import Specimen
from ..curved_clot_audit.common import read,dump,digest

def run_case(r):
 O=Path(r['task_configuration']).parent;out=Path(r['output_dir']);a=Specimen(O);s=a.s;C=a.C['fracture'];G=a.G;resume=r.get('action')=='RESUME_IMPLEMENTATION_BUG'
 if not resume:out.joinpath('accepted').mkdir();out.joinpath('attempts').mkdir();a.static_checkpoint(out)
 else:
  identity=np.load(out/'native_identity.npz');assert np.array_equal(identity['u_coords'],a.coords) and np.array_equal(identity['d_coords'],s.dcoords) and np.array_equal(identity['u_dofmap'],a.V.dofmap.list) and np.array_equal(identity['d_dofmap'],a.dc)
 elapsed_offset=time.time()-read(out/'PATH_STARTED.json')['epoch'] if resume else 0.
 start=time.monotonic()-elapsed_offset;deadline=start+C['scientific_wall_seconds']
 if not resume:dump(out/'PATH_STARTED.json',dict(epoch=time.time(),scientific_wall_seconds=C['scientific_wall_seconds'],accepted_seed_initialization_in_budget=True,parameters=a.C['material'],degree_u=a.C['fracture']['degree_u']))
 else:dump(out/'RESUMED_BUDGET.json',dict(epoch=time.time(),elapsed_offset_s=elapsed_offset,remaining_wall_s=deadline-time.monotonic(),same_budget_origin=True))
 def timer(*args):
  if time.monotonic()>deadline:raise TimeoutError('Declared scientific wall budget exhausted')
 s.mechanics.solver.setMonitor(timer);s.mechanics.solver.getKSP().setMonitor(timer);s.phase.solver.setMonitor(timer)
 pairs=list(itertools.combinations(range(4),2));edges=np.unique(np.sort(a.dc[:,pairs].reshape(-1,2),axis=1),axis=0);N=len(s.d.x.array);history=[];attempts=[];strain=0.;step=C['initial_step'];fail_count=0;clear_at=None;mode='loading';stop=None;initial_d=None;excluded_front=None;initial_energy=None;accepted_u=None;accepted_d=None;last_quality=None
 def evaluate_quality(pres):
  q,react,_=s.mechanical_quality(G['force_floor_N'],G['force_scale_half_N'],G['length_scale_m'],pres);kkt,_=s.kkt();bounds=dict(lower=float(np.min(s.d.x.array-s.old.x.array)),minimum=float(s.d.x.array.min()),maximum=float(s.d.x.array.max()));E=s.energies();gates=dict(mechanics=q['free_residual_relative']<=1e-6,phase=kkt<=1e-6,bounds=bounds['lower']>=-1e-8 and bounds['minimum']>=-1e-8 and bounds['maximum']<=1+1e-8,J=q['J_min']>0,finite=q['finite'] and all(np.isfinite(v) for v in E.values()),force=q['global_force_relative']<=1e-5,moment=q['current_moment_relative']<=1e-5,BC=q['BC_error_over_L']<=1e-10,mechanical_reason=q['mechanical_reason']>0,phase_reason=s.phase.solver.getConvergedReason()>0)
  return q,react,kkt,bounds,E,gates
 def crack_metrics():
  d=s.d.x.array;high=d>=.95;selected=high[edges].all(1);e=edges[selected];graph=coo_matrix((np.ones(len(e)*2,dtype=np.int8),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(N,N)).tocsr();_,labels=connected_components(graph);seedlabels=np.unique(labels[a.seed]);connected=high&np.isin(labels,seedlabels);front=float(s.dcoords[connected,0].max()) if connected.any() else 0.;outside=initial_d<.01 if initial_d is not None else np.zeros(N,bool);new=connected&outside&(s.dcoords[:,0]>excluded_front if excluded_front is not None else False);length=max(0.,front-(excluded_front if excluded_front is not None else front));many=int(np.sum((connected[a.dc].sum(1)>=2)&(s.dcoords[a.dc,0].mean(1)>(excluded_front if excluded_front is not None else 1.))))
  return dict(seed_connected_high_damage_front_m=front,initial_exclusion_front_m=excluded_front,new_extension_m=length,new_high_damage_nodes=int(new.sum()),multi_element_count=many,d_max_outside_initial_band=float(d[outside].max()) if outside.any() else 0.,d_increment_max_outside_initial_band=float(np.max(d[outside]-initial_d[outside])) if outside.any() else 0.,diagnostic='P1 threshold graph, not detached fragment count')
 def accept(target,phase,am,q,react,kkt,bounds,E,gates):
  nonlocal initial_d,excluded_front,initial_energy,accepted_u,accepted_d,strain
  idx=len(history);d=s.d.x.array.copy();u=s.u.x.array.copy();old=s.old.x.array.copy();R=react.reshape(-1,3);F=float(R[a.top,1].sum()-R[a.bottom,1].sum()) # 2 thickness factor times half-model generalized reaction
  fullE={k:2*v for k,v in E.items()};Delta=target*.01
  if idx==0:
   initial_d=d.copy();excluded_front=float(s.dcoords[d>=.01,0].max());initial_energy=fullE['U_frac_J'];work=0.
  else:work=history[-1]['boundary_path_work_J']+.5*(history[-1]['grip_force_full_N']+F)*(Delta-history[-1]['Delta_m'])
  crack=crack_metrics();balance=work-(fullE['U_total_J']-initial_energy);row=dict(index=idx,stage=phase,Delta_over_H=target,Delta_m=Delta,grip_force_full_N=F,energy_full_specimen=fullE,initial_crack_energy_full_J=initial_energy,boundary_path_work_J=work,energy_balance_residual_J=balance,energy_balance_residual_fixed_scale=balance/G['energy_scale_full_J'],energy_scale_full_J=G['energy_scale_full_J'],quality=q,phase_KKT=kkt,bounds=bounds,gates=gates,alternate_iterations=am,crack=crack,elapsed_s=time.monotonic()-start,increment_elapsed_s=(time.monotonic()-start-history[-1]['elapsed_s']) if history else time.monotonic()-start)
  state=out/'accepted'/f'{idx:04d}.npz';np.savez_compressed(state,u=u,d=d,d_old=old,reaction_N=react,prescribed_u=a.boundary.x.array,Delta_m=Delta,Delta_over_H=target,history_lower_bound_next=d);row['checkpoint']=dict(path=str(state),sha256=digest(state));history.append(row);dump(out/'HISTORY.json',history);dump(out/'accepted'/f'{idx:04d}.json',row);s.old.x.array[:]=d;accepted_u=u;accepted_d=d;strain=target;print('ACCEPT',idx,phase,target,'F',F,'extension',crack['new_extension_m'],'res',q['free_residual_relative'],'KKT',kkt,'wall',row['elapsed_s'],flush=True)
  return row
 try:
  if resume:
   history=read(out/'HISTORY.json');attempts=read(out/'ATTEMPTS.json');base=np.load(history[0]['checkpoint']['path']);initial_d=base['d'].copy();excluded_front=history[0]['crack']['initial_exclusion_front_m'];initial_energy=history[0]['initial_crack_energy_full_J'];saved=np.load(history[-1]['checkpoint']['path']);accepted_u=saved['u'].copy();accepted_d=saved['d'].copy();strain=history[-1]['Delta_over_H'];s.u.x.array[:]=accepted_u;s.d.x.array[:]=accepted_d;s.old.x.array[:]=accepted_d;a.prescribe(strain);step=min(C['maximum_step'],(history[-1]['Delta_over_H']-history[-2]['Delta_over_H'])*C['accepted_step_growth_factor'])
   outside=initial_d<.01
   for h in history:
    oldstate=np.load(h['checkpoint']['path']);h['crack']['d_increment_max_outside_initial_band']=float(np.max(oldstate['d'][outside]-initial_d[outside]))
   dump(out/'HISTORY.json',history)
  else:
   pres=a.prescribe(0.);s.update();s.phase_solve();s.mechanics.solve();q,react,kkt,bounds,E,gates=evaluate_quality(pres)
   if not all(gates.values()):raise RuntimeError('Zero-load initialization failed: '+str(gates))
   accept(0.,'zero_load_precrack',0,q,react,kkt,bounds,E,gates)
  for attempt in range(len(attempts),C['maximum_attempts']):
   timer()
   if len(history)-1>=C['maximum_accepted_steps']:stop='ACCEPTED_STEP_BUDGET';break
   target=max(0.,strain*(1-C['unloading_fraction'])) if mode=='unloading' else min(C['load_max_Delta_over_H'],strain+step)
   if mode=='loading' and target<=strain+1e-14:stop='DECLARED_LOADING_REACHED';break
   dest=out/'attempts'/f'{attempt:04d}';dest.mkdir();s.u.x.array[:]=accepted_u;s.d.x.array[:]=accepted_d;s.old.x.array[:]=accepted_d
   dump(dest/'TRIAL_STARTED.json',dict(attempt=attempt,from_accepted=len(history)-1,target_strain=target,stage=mode,lower_bound_checkpoint=history[-1]['checkpoint']))
   # An initial guess only; accepted-state history is held fixed until gates pass.
   if strain>0:s.u.x.array[:]*=target/strain
   else:s.u.x.array.reshape(-1,3)[:,1]=target*a.coords[:,1]
   pres=a.prescribe(target);iteration_log=[];passed=False;error=None;t0=time.monotonic()
   try:
    for alternate in range(C['alternate_max_iterations']):
     timer();s.mechanics.solve();s.update();reason=int(s.mechanics.solver.getConvergedReason())
     if reason<=0:raise RuntimeError('Mechanical SNES reason '+str(reason)+'; KSP reason '+str(s.mechanics.solver.getKSP().getConvergedReason()))
     phase_reason=s.phase_solve();q,react,kkt,bounds,E,gates=evaluate_quality(pres);iteration_log.append(dict(alternate=alternate+1,mechanical_reason=reason,mechanical_iterations=int(s.mechanics.solver.getIterationNumber()),phase_reason=phase_reason,quality_after_damage=q,phase_KKT=kkt,energy_half_J=E,gates=gates));dump(dest/'ALTERNATE_ITERATIONS.json',iteration_log)
     if phase_reason<=0:raise RuntimeError('Phase VI reason '+str(phase_reason))
     if all(gates.values()):passed=True;break
    if not passed:error='ALTERNATE_ITERATION_LIMIT_OR_QUALITY_FAILURE'
   except Exception as exc:error=str(exc);(dest/'failure.log').write_text(traceback.format_exc())
   rec=dict(attempt=attempt,from_accepted=len(history)-1,from_strain=strain,target_strain=target,stage=mode,step=target-strain,PASS=passed,error=error,wall_s=time.monotonic()-t0,alternate_iterations=len(iteration_log),consecutive_bisections_before=fail_count);attempts.append(rec);dump(dest/'ATTEMPT.json',rec);dump(out/'ATTEMPTS.json',attempts)
   if passed:
    row=accept(target,mode,len(iteration_log),q,react,kkt,bounds,E,gates);fail_count=0
    if mode=='unloading':stop='POST_PATH_UNLOAD_COMPLETE';break
    loading=[h for h in history if h['stage']=='loading'];peak=max(h['grip_force_full_N'] for h in loading);row['accepted_peak_force_N']=peak;row['load_drop_fraction']=1-row['grip_force_full_N']/peak if peak>0 else 0.;peak_index=max(loading,key=lambda h:h['grip_force_full_N'])['index'];eligible=[h for h in loading if h['index']>peak_index and h['grip_force_full_N']<=.7*peak and h['crack']['new_extension_m']>G['crack_extension_gate_m'] and h['crack']['multi_element_count']>=3]
    if len(eligible)>=2 and clear_at is None:clear_at=row['index']
    dump(out/'HISTORY.json',history)
    if clear_at is not None and row['index']>=clear_at+C['post_failure_additional_accepted_steps']:mode='unloading'
    elif strain>=C['load_max_Delta_over_H']-1e-14:mode='unloading'
    step=min(C['maximum_step'],step*C['accepted_step_growth_factor'])
   else:
    np.savez_compressed(dest/'REJECTED_TRIAL.npz',u=s.u.x.array,d=s.d.x.array,d_old=accepted_d);s.u.x.array[:]=accepted_u;s.d.x.array[:]=accepted_d;s.old.x.array[:]=accepted_d;a.prescribe(strain)
    if time.monotonic()>deadline:stop='WALL_BUDGET';break
    if fail_count>=C['consecutive_bisections'] or mode=='unloading':stop='NUMERICAL_FAILURE';break
    fail_count+=1;step*=.5;print('REJECT / BISECT',attempt,error,'next step',step,flush=True)
  else:stop='ATTEMPT_BUDGET'
 except Exception as exc:
  stop='WALL_BUDGET' if isinstance(exc,TimeoutError) else 'NUMERICAL_FAILURE';(out/'PATH_EXCEPTION.log').write_text(traceback.format_exc())
 loading=[h for h in history if h['stage']=='loading'];peak=max((h['grip_force_full_N'] for h in loading),default=0.);drop=max((1-h['grip_force_full_N']/peak for h in loading if peak>0 and h['index']>=max(loading,key=lambda x:x['grip_force_full_N'])['index']),default=0.);extension=max((h['crack']['new_extension_m'] for h in history),default=0.);damage=max((h['crack']['d_increment_max_outside_initial_band'] for h in history),default=0.)
 physical='CLEAR_LOAD_BEARING_FAILURE' if clear_at is not None else ('RESOLVED_CRACK_GROWTH' if extension>0 and any(h['crack']['multi_element_count']>=3 for h in history) else ('DAMAGE_INITIATED' if damage>1e-3 else 'NO_DAMAGE'))
 outcome='NUMERICAL_FAILURE' if stop=='NUMERICAL_FAILURE' else ('NOT_OBSERVED_WITHIN_BUDGET' if 'BUDGET' in stop and physical!='CLEAR_LOAD_BEARING_FAILURE' else physical)
 result=dict(success=bool(history) and stop!='NUMERICAL_FAILURE',solve_status='NUMERICAL_QUALITY_FAILED' if stop=='NUMERICAL_FAILURE' else 'ACCEPTED_PATH_'+stop,fracture_outcome=outcome,accepted_state_physical_outcome=physical,stop_reason=stop,accepted_states=len(history),attempts=len(attempts),maximum_Delta_over_H=max((h['Delta_over_H'] for h in history),default=0.),new_extension_m=extension,peak_full_force_N=peak,post_peak_drop_fraction=drop,initial_crack_energy_full_J=initial_energy,maximum_abs_energy_balance_scaled=max((abs(h['energy_balance_residual_fixed_scale']) for h in history),default=None),wall_s=time.monotonic()-start,material_seconds=s.material_seconds,material_calls=s.material_calls,clear_at=clear_at,full_native_checkpoints=True,method_id=C['method_id'],discretization_scope=a.C['fracture']['discretization_scope'],failed_gates=['path_continuation'] if stop=='NUMERICAL_FAILURE' else [])
 dump(out/'PATH_SUMMARY.json',result)
 if history:
  fields=['index','stage','Delta_over_H','Delta_m','grip_force_full_N','boundary_path_work_J','energy_balance_residual_J','energy_balance_residual_fixed_scale','phase_KKT','elapsed_s','increment_elapsed_s']
  with (out/'CURVES.csv').open('w',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=fields+list(history[0]['energy_full_specimen'])+['new_extension_m','free_residual_relative','J_min','J_max']);writer.writeheader()
   for h in history:writer.writerow({**{k:h[k] for k in fields},**h['energy_full_specimen'],'new_extension_m':h['crack']['new_extension_m'],**{k:h['quality'][k] for k in ['free_residual_relative','J_min','J_max']}})
 result['response']=dict(path_summary=str(out/'PATH_SUMMARY.json'),history=str(out/'HISTORY.json'),outcome=outcome,accepted_state_physical_outcome=physical) if result['success'] else None
 if not result['success']:result['diagnostic_response']=dict(last_accepted_state=len(history)-1,accepted_physical_outcome=physical,summary=str(out/'PATH_SUMMARY.json'))
 return result
