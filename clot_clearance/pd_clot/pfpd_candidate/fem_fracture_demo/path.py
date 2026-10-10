"""A single irreversible path; only fully checked states become history."""
import time,traceback,itertools,csv
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from petsc4py import PETSc
from .specimen import Specimen
from .monolithic import CoupledVI
from .api import create_case
from .loads import displacement_schedule
from .native import vector
from ..curved_clot_audit.common import read,dump,digest

def run(O,request,hemisphere=False,resume=False):
 O=Path(O);C=read(O/'CONFIGURATION.json');cfg=C['fracture'];schedule=displacement_schedule(request,hemisphere);out=O/('HEMISPHERE_WALL_CASE' if hemisphere else 'MODE_I_CASE')
 if resume:
  assert read(out/'REQUEST.json')==request and not (out/'PATH_SUMMARY.json').exists(), 'Only an interrupted, unfinished identical case may resume'
 else:out=create_case(request,out)
 assert read(O/'MIXED_FORMULATION_CHECK.json')['PASS'] and read(O/'SMALL_COUPLED_CHECKS.json')['PASS'] and read(O/'MONOLITHIC_SMALL_CHECK.json')['PASS']
 if hemisphere:
  assert read(O/'MODE_I_CASE/PATH_SUMMARY.json')['numerical_gate_pass']
  from .wall import Hemisphere
  a=Hemisphere(O);maxload=1.;step=C['hemisphere']['initial_step'];maxstep=C['hemisphere']['maximum_step'];L=a.R;force_scale=1040*a.R*a.R;efactor=1.
 else:a=Specimen(O,options={'ksp_type':'preonly','pc_type':'lu','pc_factor_mat_solver_type':'mumps'});maxload=cfg['max_Delta_over_H'];step=cfg['initial_step'];maxstep=cfg['maximum_step'];L=.04;force_scale=.0624;efactor=4.
 maxload=schedule['maximum_parameter'];maximum_accepted_steps=min(cfg['maximum_accepted_steps'],schedule['max_steps']);s=a.s
 if not resume:(out/'accepted').mkdir();(out/'attempts').mkdir();a.checkpoint_identity(out)
 else:
  identity=np.load(out/'native_identity.npz');assert np.array_equal(identity['u_coords'],a.coords) and np.array_equal(identity['d_coords'],s.dcoords)
 vi=None;start=time.perf_counter();timers=dict(mechanics_s=0.,phase_s=0.,coupled_VI_s=0.,quality_s=0.,IO_s=0.)
 PETSc.Log.begin();stages={k:PETSc.Log.Stage(k) for k in ['mechanics','phase','coupled_VI','quality','IO']}
 def timed(key,fun):
  st=time.perf_counter()
  with stages[key]:value=fun()
  timers[key+'_s']+=time.perf_counter()-st;return value
 def mechanical_solve():
  # Align the absolute Newton tolerance with the declared relative force gate.
  # A final independent residual check after damage still decides acceptance.
  s.update(tangent=False);res=vector(s.Ru);reference=max(np.linalg.norm(res[s.fixed_flat]),force_scale*1e-6)
  s.mechanics.solver.setTolerances(atol=.5e-6*reference,rtol=1e-10,max_it=80)
  return s.mechanics.solve()
 edges=np.unique(np.sort(a.dc[:,list(itertools.combinations(range(4),2))].reshape(-1,2),axis=1),axis=0);N=len(s.d.x.array)
 history=[];attempts=[];initial_d=None;initial_front=.013 if not hemisphere else 0.;accepted_u=s.u.x.array.copy();accepted_d=s.d.x.array.copy();load=0.;fails=0;clear_at=None;stop=None
 def quality(pres):
  q,react,res=s.mechanical_quality(force_scale*1e-6,force_scale,L,pres);kkt,_=s.kkt();E=s.energies();d=s.d.x.array
  bounds=dict(minimum=float(d.min()),maximum=float(d.max()),irreversibility_min=float((d-s.old.x.array).min()))
  gates=dict(mechanics=q['free_residual_relative']<=1e-6,KKT=kkt<=1e-6,bounds=bounds['minimum']>=-1e-8 and bounds['maximum']<=1+1e-8 and bounds['irreversibility_min']>=-1e-8,J=q['J_min']>0,finite=q['finite'] and all(np.isfinite(v) for v in E.values()),BC=q['BC_error_over_L']<=1e-8,force=q['global_force_relative']<=1e-5,moment=q['current_moment_relative']<=1e-5)
  return q,react,kkt,bounds,E,gates
 def metrics():
  d=s.d.x.array;high=d>=.95;e=edges[high[edges].all(1)];graph=coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(N,N)).tocsr();_,labels=connected_components(graph)
  if len(a.seed):connected=high&np.isin(labels,labels[a.seed])
  elif high.any():
   counts=np.bincount(labels[high]);connected=high&(labels==int(np.argmax(counts)))
  else:connected=high
  front=float(s.dcoords[connected,0].max()) if connected.any() else initial_front
  extent=max(0.,front-initial_front) if not hemisphere else (float(np.ptp(s.dcoords[connected,0])) if connected.any() else 0.)
  outside=initial_d<.01 if initial_d is not None else np.zeros(N,bool);new=connected&outside
  many=int(np.sum((connected[a.dc].sum(1)>=2)&(s.dcoords[a.dc,0].mean(1)>initial_front))) if not hemisphere else int(np.sum(connected[a.dc].sum(1)>=2))
  return dict(high_damage_threshold=.95,front_m=front,extension_m=extent,multi_element_count=many,new_high_nodes_outside_initial_band=int(new.sum()),d_max_outside_initial_band=float(d[outside].max()) if outside.any() else 0.,d_increment_max=float((d-initial_d).max()) if initial_d is not None else 0.,high_nodes=int(high.sum()),high_damage_mean_volume_fraction=float(np.mean(high[a.dc])),interpretation='Nodal threshold connectivity, not a detached-fragment count')
 def accept(target,iteration,q,react,kkt,bounds,E,gates):
  nonlocal initial_d,initial_front,accepted_u,accepted_d,load
  idx=len(history);d=s.d.x.array.copy();u=s.u.x.array.copy()
  if idx==0:
   initial_d=d.copy()
   if not hemisphere:initial_front=float(s.dcoords[d>=.95,0].max())
  F=float(np.dot(react.reshape(-1,3)[:,0],a.coords[:,0])/a.R) if hemisphere else float(2*react.reshape(-1,3)[a.top,1].sum())
  displacement=target*(a.R if hemisphere else .01);E={k:efactor*v for k,v in E.items()};work=0. if not history else history[-1]['boundary_work_J']+.5*(F+history[-1]['reaction_N'])*(displacement-history[-1]['displacement_m'])
  row=dict(index=idx,load_parameter=target,wall_lambda=1+target if hemisphere else None,Delta_over_H=target if not hemisphere else None,displacement_m=displacement,reaction_N=F,energies=E,quality=q,KKT=kkt,bounds=bounds,gates=gates,crack=metrics(),inner_iterations=iteration,boundary_work_J=work,elapsed_s=time.perf_counter()-start)
  row['work_minus_energy_increment_J']=work-(E['U_total_J']-(history[0]['energies']['U_total_J'] if history else E['U_total_J']))
  fq=s.Fexpr.eval(s.domain,s.entities).reshape(-1,3,3);row['raw_material']=dict(Wiso_max_Pa=float(s.Wi.x.array.max()),W0_max_Pa=float(s.W.x.array.max()),principal_stretch_max=float(np.linalg.svd(fq,compute_uv=False).max()))
  path=out/'accepted'/f'{idx:04d}.npz'
  def save():np.savez_compressed(path,u=u,d=d,d_old=s.old.x.array,reaction_N=react,prescribed_u=pres,load_parameter=target,Wiso_cellmean_Pa=s.Wi.x.array[s.qids].mean(1),W0_cellmean_Pa=s.W.x.array[s.qids].mean(1),J_min_cell=np.linalg.det(fq).reshape(s.nc,s.nq).min(1));return None
  timed('IO',save);row['checkpoint']=dict(path=str(path),sha256=digest(path));history.append(row);dump(out/'accepted'/f'{idx:04d}.json',row);dump(out/'HISTORY.json',history)
  s.old.x.array[:]=d;accepted_u=u;accepted_d=d;load=target
  print('ACCEPT',idx,'load',target,'reaction',F,'extension',row['crack']['extension_m'],'d_increment',row['crack']['d_increment_max'],'res',q['free_residual_relative'],'KKT',kkt,'seconds',row['elapsed_s'],flush=True)
 try:
  if resume:
   history=read(out/'HISTORY.json');attempts=read(out/'ATTEMPTS.json') if (out/'ATTEMPTS.json').exists() else [];zz=np.load(history[-1]['checkpoint']['path']);initial_d=np.load(history[0]['checkpoint']['path'])['d'];initial_front=float(s.dcoords[initial_d>=.95,0].max()) if not hemisphere else 0.;accepted_u=zz['u'].copy();accepted_d=zz['d'].copy();s.u.x.array[:]=accepted_u;s.d.x.array[:]=accepted_d;s.old.x.array[:]=accepted_d;load=history[-1]['load_parameter'];pres=a.prescribe(load)
   step=read(O/'SOLVER_DIAGNOSTIC_RESTART.json')['restore_step']
   s.lower.array[:]=accepted_d;s.lower.array[a.seed]=1.;s.phase.solver.setVariableBounds(s.lower,s.upper)
   for trial in sorted((out/'attempts').glob('*')):
    if not (trial/'ATTEMPT.json').exists():
     pending=read(trial/'START.json');step=pending['target']-load;rec=dict(attempt=int(trial.name),from_index=len(history)-1,target=pending['target'],PASS=False,error='INTERRUPTED_FOR_SOLVER_DIAGNOSTICS',inner_count=0,wall_s=None);dump(trial/'ATTEMPT.json',rec);attempts.append(rec)
   dump(out/'ATTEMPTS.json',attempts);dump(out/'RESUME.json',dict(from_index=len(history)-1,step=step,reason='Add Newton diagnostics; physical configuration and accepted history unchanged',config_sha256=digest(O/'CONFIGURATION.json')))
  else:
   pres=a.prescribe(0.);s.update();timed('phase',s.phase_solve);timed('mechanics',mechanical_solve);q,react,kkt,bounds,E,gates=timed('quality',lambda:quality(pres))
   if not all(gates.values()):raise RuntimeError('Initial state gates: '+str(gates))
   accept(0.,0,q,react,kkt,bounds,E,gates)
  for attempt in range(len(attempts),cfg['maximum_attempts']):
   if len(history)>=maximum_accepted_steps+1:stop='STEP_LIMIT';break
   target=min(maxload,load+step)
   if hemisphere and load<.5<target:target=.5
   if target<=load+1e-13:stop='DECLARED_PATH_COMPLETE';break
   dest=out/'attempts'/f'{attempt:04d}';dest.mkdir();s.u.x.array[:]=accepted_u;s.d.x.array[:]=accepted_d;s.old.x.array[:]=accepted_d
   if load>0:s.u.x.array[:]*=target/load
   elif hemisphere:s.u.x.array.reshape(-1,3)[:,0]=target*a.coords[:,0]
   else:s.u.x.array.reshape(-1,3)[:,1]=target*a.coords[:,1]
   pres=a.prescribe(target);inner=[];passed=False;error=None;t0=time.perf_counter();dump(dest/'START.json',dict(from_index=len(history)-1,target=target,history_sha256=history[-1]['checkpoint']['sha256']))
   newton=[]
   def monitor(sn,iteration,norm):
    row=dict(iteration=iteration,norm=float(norm),seconds=time.perf_counter()-t0);newton.append(row);dump(dest/'NEWTON_HISTORY.json',newton)
    if iteration%5==0:print('NEWTON',attempt,iteration,'norm',norm,flush=True)
   s.mechanics.solver.cancelMonitor();s.mechanics.solver.setMonitor(monitor)
   try:
    from .predictor import predict
    predictor=timed('mechanics',lambda:predict(s,accepted_u,pres));dump(dest/'TANGENT_PREDICTOR.json',predictor)
    for it in range(cfg['alternate_max_iterations']):
     timed('mechanics',mechanical_solve);s.update();mr=int(s.mechanics.solver.getConvergedReason())
     if mr<=0:raise RuntimeError('Mechanical SNES reason '+str(mr))
     pr=timed('phase',s.phase_solve);q,react,kkt,bounds,E,gates=timed('quality',lambda:quality(pres))
     if pr<=0:raise RuntimeError('Damage VI reason '+str(pr))
     record=dict(iteration=it+1,method='STAGGERED',mechanical_reason=mr,phase_reason=pr,quality_after_damage=q,KKT=kkt,energy_J=E['U_total_J'],gates=gates);inner.append(record)
     if not all(gates.values()) and ((it==2 and q['free_residual_relative']>1e-3) or (it>=19 and (it+1)%20==0 and q['free_residual_relative']>1e-4)):
      if vi is None:vi=CoupledVI(s)
      bu=s.u.x.array.copy();bd=s.d.x.array.copy();before=E['U_total_J'];vr=timed('coupled_VI',lambda:vi.solve(pres));qv,rv,kv,bv,Ev,gv=timed('quality',lambda:quality(pres))
      keep=vr>0 and all(gv.values()) and Ev['U_total_J']<=before+max(abs(before),1e-12)*1e-8
      inner.append(dict(iteration=it+1,method='PETSC_COUPLED_VI',reason=vr,quality=qv,KKT=kv,energy_J=Ev['U_total_J'],kept=keep,gates=gv))
      if keep:q,react,kkt,bounds,E,gates=qv,rv,kv,bv,Ev,gv
      else:s.u.x.array[:]=bu;s.d.x.array[:]=bd;s.update()
     dump(dest/'INNER_ITERATIONS.json',inner)
     if all(gates.values()):passed=True;break
     if (it+1)%20==0:print('INNER',attempt,it+1,'load',target,'res',q['free_residual_relative'],'KKT',kkt,flush=True)
    if not passed:error='ALTERNATE_ITERATION_LIMIT'
   except Exception as exc:error=str(exc);(dest/'EXCEPTION.log').write_text(traceback.format_exc())
   rec=dict(attempt=attempt,from_index=len(history)-1,target=target,PASS=passed,error=error,inner_count=len(inner),wall_s=time.perf_counter()-t0);attempts.append(rec);dump(dest/'ATTEMPT.json',rec);dump(out/'ATTEMPTS.json',attempts)
   if passed:
    accept(target,len(inner),q,react,kkt,bounds,E,gates);fails=0;peak=max(history,key=lambda h:h['reaction_N']);eligible=[h for h in history if h['index']>peak['index'] and h['reaction_N']<=.8*peak['reaction_N'] and h['crack']['extension_m']>(0 if hemisphere else cfg['extension_gate_m']) and h['crack']['multi_element_count']>=3]
    if len(eligible)>=2 and eligible[-1]['reaction_N']<=eligible[-2]['reaction_N']:clear_at=len(history)-1;stop='CLEAR_FRACTURE_DEMONSTRATED';break
    step=min(maxstep,step*1.2)
   else:
    np.savez_compressed(dest/'REJECTED_TRIAL.npz',u=s.u.x.array,d=s.d.x.array,d_old=accepted_d);s.u.x.array[:]=accepted_u;s.d.x.array[:]=accepted_d;s.old.x.array[:]=accepted_d;a.prescribe(load);fails+=1;step*=.5
    print('REJECT',attempt,error,'next step',step,flush=True)
    if fails>=cfg['consecutive_bisections'] or step<cfg['minimum_step']:stop='NUMERICAL_FAILURE';break
  else:stop='ATTEMPT_LIMIT'
 except Exception as exc:stop='NUMERICAL_FAILURE';(out/'PATH_EXCEPTION.log').write_text(traceback.format_exc())
 peak=max(history,key=lambda h:h['reaction_N']) if history else None
 post=[h for h in history if h['index']>peak['index']] if peak else []
 growth=any(h['crack']['extension_m']>0 and h['crack']['multi_element_count']>=3 for h in history);damage=any(h['crack']['d_increment_max']>1e-3 for h in history)
 status='CLEAR_FRACTURE_DEMONSTRATED' if clear_at is not None else ('RESOLVED_DAMAGE_GROWTH_NO_CLEAR_CAPACITY_DROP' if growth or damage else 'NO_FRACTURE_WITHIN_DECLARED_PATH')
 if stop=='NUMERICAL_FAILURE':status='NUMERICAL_FAILURE'
 result=dict(status=status,stop_reason=stop,numerical_gate_pass=bool(history) and stop!='NUMERICAL_FAILURE' and all(all(h['gates'].values()) for h in history),accepted_states=len(history),attempts=len(attempts),clear_at=clear_at,peak_force_N=peak['reaction_N'] if peak else None,peak_index=peak['index'] if peak else None,post_peak_drop_fraction=max([1-h['reaction_N']/peak['reaction_N'] for h in post],default=0),post_peak_states=len(post),maximum_load_parameter=load,maximum_extension_m=max([h['crack']['extension_m'] for h in history],default=0),timers=timers,material_seconds=s.material_seconds,material_calls=s.material_calls,wall_s=time.perf_counter()-start,no_failed_history_commit=True,configuration_sha256=digest(O/'CONFIGURATION.json'))
 profiles={}
 result['assembly_seconds']=dict(mechanics=s.assembly_seconds,coupled_VI=vi.assembly_seconds if vi is not None else None,scope='Exclusive FE assembly within solver callbacks in this process; material evaluation separately timed. Nested inside broad solve times.')
 for name,stage in stages.items():
  profiles[name]={}
  for event in ['KSPSolve','PCSetUp','MatLUFactorNum','SNESFunctionEval','SNESJacobianEval']:
   try:profiles[name][event]=dict(PETSc.Log.Event(event).getPerfInfo(stage))
   except Exception:pass
 dump(out/'PETSC_PROFILE.json',profiles);dump(out/'PATH_SUMMARY.json',result)
 for file,columns in [('load_history.csv',['index','load_parameter','displacement_m','reaction_N','inner_iterations','elapsed_s']),('reaction_displacement.csv',['index','displacement_m','reaction_N','boundary_work_J','work_minus_energy_increment_J']),('damage_metrics.csv',['index','load_parameter','extension_m','multi_element_count','d_increment_max','high_nodes'])]:
  with (out/file).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=columns);w.writeheader()
   for h in history:w.writerow({k:{**h,**h['crack']}[k] for k in columns})
 print('PATH_SUMMARY',result,flush=True);return result
