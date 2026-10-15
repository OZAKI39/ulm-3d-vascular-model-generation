"""Read-only recovery of already saved evidence after native atom-loss exception.
Never loads LAMMPS or reconstructs a restart; never advances a physical step.
"""
from common import *
import sys,csv,re
sys.path.insert(0,str(R/'analysis/softening_feasibility_v1'))
from model import SofteningLPS
def main():
 p=O/'full';rows=list(csv.DictReader((p/'step_scalars.csv').open()));cols=list(rows[0]);arr=np.array([[float(r[c]) for c in cols] for r in rows]);tab={c:arr[:,i] for i,c in enumerate(cols)}
 tx=json.loads((p/'TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());ref=np.load(O/'provenance/reference.npz');ids=np.load(O/'provenance/341_boundary_ids.npz');b=ref['bonds'];model=SofteningLPS(ref)
 np.savez_compressed(p/'trajectory.npz',series=arr,columns=cols,reference_X=ref['X'],bonds=b,bond_ids=np.arange(1,len(b)+1),**{k:ids[k] for k in ids.files})
 paths=sorted((p/'keyframes').glob('*.npz'));last=np.load(paths[-1]);assert int(last['step'])==int(tab['step'][-1])==2755
 xpred=last['x'].copy();ii=ids['free_ids']-1;xpred[ii]+=BASE*(last['velocity'][ii]+BASE/(2*RHO*V)*last['force'][ii]);outside=np.argwhere((xpred<-.00025)|(xpred>=.00025));details=[]
 for i,axis in outside:
  details.append(dict(atom_id=int(i+1),axis='xyz'[axis],last_accepted_x_m=last['x'][i],last_accepted_v_m_s=last['velocity'][i],last_accepted_F_N=last['force'][i],predicted_next_drift_x_m=xpred[i],outside_amount_m=float(abs(xpred[i,axis])-.00025),prediction_only_not_saved_native_candidate=True))
 np.savez_compressed(p/'PREDICTED_NEXT_DRIFT_NOT_ACCEPTED.npz',x=xpred,from_time_s=float(last['time']),candidate_time_s=float(last['time'])+BASE,free_ids=ids['free_ids'])
 dump(O/'NATIVE_ATOM_LOSS_DIAGNOSIS.json',dict(status='STOPPED',last_accepted_time_s=float(last['time']),last_accepted_step=int(last['step']),attempted_time_s=float(last['time'])+BASE,attempted_step=int(last['step'])+1,native_observation='Owned atom identity not unique across ranks; subsequent write_restart reports 340 vs 341',domain_bounds_m=[[-.00025,.00025]]*3,boundary='f f f, inherited original domain',next_VV_drift_outside_domain=details,diagnosis='A free point crosses fixed nonperiodic box during drift; native exchange removes it. Supported by last exact FP64 state, analytic next drift, fixed box, and native count=340. Missing native candidate coordinates unavailable after exception.',not_a_duplicate_ID_inference='Error wording covers both duplicate and missing ID; observed count 340 indicates missing atom',last_complete_restart='Latest transactions/accepted_before_trial_*.restart, not final 687.75us; do not synthesize restart from NPZ',restart_after_atom_loss='FAILED, no valid final restart',no_continuation_or_domain_change=True))
 keys=[];maxerr={k:0. for k in ['force_N','U_J','D_increment_J','kappa','d','native_KE_J','contact_force_N','contact_U_J','force_sum_N','torque_Nm','fixed_m','loading_x_m','loading_v_m_s','decomposition_J','failed_tension_N','d_rollback','kappa_rollback']}
 # Saved full-keyframe independent verification; do not fabricate lost online maxima.
 from scipy.spatial import cKDTree
 c=15*(18*K/(np.pi*HORIZON**4))*V**2/HORIZON;codes=b[:,0]*341+b[:,1]
 for path in paths:
  s=np.load(path);step=int(s['step']);t=float(s['time']);m=model.evaluate(s['x'],s['d']);kk,dd=model.update(s['x'],s['kappa'],s_complete=SC)
  pairs=cKDTree(s['x']).query_pairs(DX,output_type='ndarray');pairs=pairs[~np.isin(pairs[:,0]*341+pairs[:,1],codes[s['d']<1])];fc=np.zeros_like(s['x']);Uc=0.
  if len(pairs):
   v=s['x'][pairs[:,1]]-s['x'][pairs[:,0]];r=np.linalg.norm(v,axis=1);forces=c*(r-DX)[:,None]*v/r[:,None];np.add.at(fc,pairs[:,0],forces);np.add.at(fc,pairs[:,1],-forces);Uc=.5*c*np.sum((r-DX)**2)
  checks=dict(force_N=abs(m['force']-s['material_force']).max(),U_J=abs(m['U']-s['scalars'][0]),kappa=abs(kk-s['kappa']).max(),d=abs(dd-s['d']).max(),native_KE_J=abs(.5*RHO*V*np.sum(s['velocity']**2)-s['K']),contact_force_N=abs(fc-(s['force']-s['material_force'])).max(),contact_U_J=abs(Uc-s['Uc']),force_sum_N=np.linalg.norm(s['force'].sum(axis=0)),torque_Nm=np.linalg.norm(np.cross(s['x'],s['force']).sum(axis=0)),failed_tension_N=np.max(abs(m['t'][s['d']==1]),initial=0))
  u=30e-6*(np.cos(np.pi*(0 if t<.00035 else min((t-.00035)/.0003,1)))-np.cos(np.pi*min(t/.0003,1)))
  vel=60e-6*np.pi/(2*.0003)*(np.sin(np.pi*min(t/.0003,1))-np.sin(np.pi*(0 if t<.00035 else min((t-.00035)/.0003,1))))
  for key,sign in [('left_ids',-1),('right_ids',1)]:
   nodes=ids[key]-1;xp=ref['X'][nodes].copy();xp[:,0]+=sign*u;vp=np.zeros_like(xp);vp[:,0]=sign*vel
   maxerr['loading_x_m']=max(maxerr['loading_x_m'],float(abs(xp-s['x'][nodes]).max()));maxerr['loading_v_m_s']=max(maxerr['loading_v_m_s'],float(abs(vp-s['velocity'][nodes]).max()))
  for k,v in checks.items():maxerr[k]=max(maxerr[k],float(v))
  reason='saved_event_or_damage_record'
  if step==0:reason='initial'
  elif step==340:reason='first_softening'
  elif step==370:reason='first_complete_failure'
  elif step==930:reason='first_active_contact'
  elif any(abs(t-a)<1e-15 for a in [.0003,.00035,.00065,.0007]):reason='physical_phase_transition'
  elif step==2755:reason='last_accepted_before_atom_loss'
  elif 2562<=step<=2569:reason='accepted_event_substep'
  keys.append(dict(step=step,time_s=t,reason=reason,path=str(path.relative_to(p))))
 maxerr['decomposition_J']=float(abs(tab['remainder_J']).max())
 # D_increment global extrema, m drift and historical monotonicity passed online,
 # but maxima were in memory when finalization failed. Mark unavailable explicitly.
 for k in ['D_increment_J','fixed_m','d_rollback','kappa_rollback']:maxerr[k]=None
 dump(p/'KEYFRAME_INDEX.json',keys)
 firstsoft=np.load(p/'keyframes/step_000340.npz');firstfail=np.load(p/'keyframes/step_000370.npz')
 events=dict(first_softening=dict(step=340,time_s=float(firstsoft['time']),bond_ids=np.flatnonzero(firstsoft['d']>0)+1),first_complete_failure=dict(step=370,time_s=float(firstfail['time']),bond_ids=np.flatnonzero(firstfail['d']==1)+1),unique_failed_ids=np.flatnonzero(last['d']==1)+1,events='Per-step counts retained in CSV; full all-ID event stream was in memory and not flushed after native exception. Saved keyframes/transactions retain exact identities there.')
 dump(p/'FIRST_SOFTENING_AND_FAILURE_EVENTS.json',events)
 audit=dict(status='FAILED_ADMISSION',executed_steps=2755,accepted_substeps=tx['accepted'],rejected_trials=tx['rejected'],minimum_accepted_h_s=float(tab['dt_s'][1:].min()),last_admitted_step=2755,final_time_s=float(last['time']),attempted_stop_time_s=float(last['time'])+BASE,planned_end_time_s=.0007,stop_reasons=['NATIVE_ATOM_LOSS_FIXED_BOX: 340 vs 341 at attempted 688us'],max_errors=maxerr,max_errors_scope='Independent saved keyframes except CSV equation/momentum/decomposition maxima; exact per-step online maxima lost during finalization exception. All accepted steps passed frozen online gates.',online_all_accepted_gates_passed=True,final_ledger=dict(zip(cols,arr[-1])),maxima={c:float(tab[c].max()) for c in cols},native_endpoint_topology_matches=None,native_topology_note='Every accepted step counted NB-failed and checked original IDs; final direct gather_bonds comparison could not run after loss. Last FP64 final d and native bond count agree.',independent_model_checked_every_step=True,full_states_saved=len(keys),contact_status='OBSERVED',missing_final_artifacts=['valid final restart','full sparse original-bond event stream','exact global max_errors','full process peak RSS'],scope_of_recovery='Read-only saved artifacts, no physical rerun')
 dump(p/'DYNAMIC_AUDIT.json',audit);dump(p/'NUMERICAL_SAFETY_STOP_RECORD.json',audit)
 start=(p/'FORMAL_CASE_STARTED.json').stat().st_mtime_ns/1e9;end=(p/'step_scalars.csv').stat().st_mtime_ns/1e9;final=(p/'ROLLBACK_STATE_CONSISTENCY.json').stat().st_mtime_ns/1e9
 dump(p/'RUN_RESOURCES.json',dict(main_wall_s=end-start,main_wall_measurement='Filesystem realtime from start lock to final flushed accepted CSV row; approximate, not lost perf_counter value',wall_to_failure_evidence_s=final-start,main_CPU_s=None,peak_RSS_bytes=None,missing_measurement_reason='Native exception followed by failed final gather/restart interrupted resource-summary persistence; no invented RSS and no rerun to replace data',MPI_ranks=1,threads=1,actual_dynamic_cases=1,actual_physical_steps=tx['accepted'],completed_rejected_trials=tx['rejected'],additional_incomplete_native_trial=1,GPU_used=False))
 print(json.dumps(dict(last_time_us=float(last['time'])*1e6,attempted_time_us=(float(last['time'])+BASE)*1e6,points=details,resources=json.loads((p/'RUN_RESOURCES.json').read_text())),default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x,indent=2))
if __name__=='__main__':main()
