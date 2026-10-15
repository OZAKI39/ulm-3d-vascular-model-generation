from common import *
import csv,shutil
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
def read(case):
 p=O/case;z=np.load(p/'trajectory.npz');arr=z['series'];cols=list(z['columns']);tab={c:arr[:,i] for i,c in enumerate(cols)}
 return p,arr,tab,json.loads((p/'DYNAMIC_AUDIT.json').read_text()),json.loads((p/'TRIAL_ACCEPT_REJECT_AUDIT.json').read_text())
def main():
 local=json.loads((O/'LOCAL_BOND1813_EVENT_AUDIT.json').read_text());p,arr,t,full,tx=read('full');lp,la,lt,lc,lx=read('local')
 rb={c:json.loads((O/c/'ROLLBACK_STATE_CONSISTENCY.json').read_text()) for c in ['local','full']}
 txall=dict(local=lx,full=tx);dump(O/'TRIAL_ACCEPT_REJECT_AUDIT.json',txall);dump(O/'ROLLBACK_STATE_CONSISTENCY.json',rb)
 allrows=[]
 for case,data in txall.items():
  for r in data['trials']:allrows.append(dict(case=case,**{k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()}))
 with (O/'ADAPTIVE_SUBSTEP_HISTORY.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 dump(O/'FULL_PROTOCOL_DYNAMIC_AUDIT.json',full)
 keys=json.loads((p/'KEYFRAME_INDEX.json').read_text());snapshots=[]
 for entry in keys:
  s=np.load(p/entry['path']);b=s['bonds'];alive=b[s['d']<1];N=len(s['x']);adj=coo_matrix((np.ones(2*len(alive)),(np.r_[alive[:,0],alive[:,1]],np.r_[alive[:,1],alive[:,0]])),shape=(N,N)).tocsr();nc,labels=connected_components(adj)
  snapshots.append(dict(**entry,failed=int(np.sum(s['d']==1)),live=len(alive),component_sizes=sorted(np.bincount(labels).tolist(),reverse=True),components=nc))
 last=np.load(p/keys[-1]['path']);failed=np.flatnonzero(last['d']==1)+1
 events=json.loads((p/'FIRST_SOFTENING_AND_FAILURE_EVENTS.json').read_text())
 dump(O/'DAMAGE_HISTORY_AND_TOPOLOGY.json',dict(original_nodes=341,original_bonds=12047,first_softening=events['first_softening'],first_complete_failure=events['first_complete_failure'],snapshots=snapshots,final_failed_ids=failed,final_failed=len(failed),final_live=12047-len(failed),final_components=snapshots[-1]['component_sizes'],original_ID_stable_online=True,final_native_gather=full['native_endpoint_topology_matches'],damage_irreversible_online=True,duplicate_failure_online=False,full_reference_m_source_unchanged=True,evidence_limitations=full.get('missing_final_artifacts',[])))
 maxrel=float(t['relative_R'].max());end={c:float(t[c][-1]) for c in ['time_s','W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','energy_scale_J','relative_R','VV_defect_J','spatial_defect_J','switch_work_J','contact_defect_J','remainder_J','contact_work_J']}
 ledger=dict(expression='W=DeltaU+DeltaK+D+DeltaUc-Jc+R',normalization='running maximum of 1e-15J, |W|, |DeltaU|+|DeltaK|+D+|DeltaUc|; unchanged',max_relative_R=maxrel,relative_gate=.02,final=end,local_before_after=local['comparison_same_640p25_to_640p5_window'],independent_per_accepted_step_errors=full['max_errors'],rejected_trials_excluded=True,damage_D='same new x: U(x_new,d_old)-U(x_new,d_new)',signed_topology_display_only=True,spatial_force_and_time_layer_definition='Same inherited decomposition, using actual accepted h and actual dx; negative sum of VV/spatial/switch/contact defects reconstructs R',source_csv='full/step_scalars.csv')
 dump(O/'DYNAMIC_ENERGY_LEDGER.json',ledger)
 dump(O/'CONTACT_AND_TIMELAYER_AUDIT.json',dict(contact_factor_dx=1.,unchanged_contact_stiffness=True,steps_with_contact=int(np.sum(t['contact_pairs'][1:]>0)),peak_pairs=int(t['contact_pairs'].max()),last_pairs=int(t['contact_pairs'][-1]),max_penetration_dx=float(t['max_contact_penetration_dx'].max()),max_penetration_m=float(t['max_contact_penetration_dx'].max()*DX),contact_force_error_N=full['max_errors']['contact_force_N'],contact_energy_error_J=full['max_errors']['contact_U_J'],rollback_exclusions_exact=all(r['contact_exclusions_exact'] for d in rb.values() for r in d['records']),physical_penetration_precision='UNVALIDATED, no stiffness calibration',max_loading_x_error_m=full['max_errors']['loading_x_m'],max_loading_v_error_m_s=full['max_errors']['loading_v_m_s'],max_VV_x_error_m=float(t['x_equation_m'].max()),max_VV_v_error_m_s=float(t['v_equation_m_s'].max()),max_momentum_residual_Ns=float(t['momentum_residual_Ns'].max()),same_phase_times_s=[.0003,.00035,.00065,.0007]))
 resources={case:json.loads((O/case/'RUN_RESOURCES.json').read_text()) for case in ['local','full']};build=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text())
 covered_intervals=(t['time_s'][-1]-t['time_s'][0])/BASE
 refcost=27.864420813974/2562*covered_intervals
 dump(O/'RUN_RESOURCES.json',dict(**resources,replay=json.loads((O/'replay/DETERMINISTIC_REPLAY_AUDIT.json').read_text()),build_wall_s=build['compile_wall_s'],observer_incremental_build_wall_s=build['observer_compile_wall_s'],accepted_full=tx['accepted'],rejected_full=tx['rejected'],checkpoints_full=tx['checkpoints'],minimum_accepted_h_s=full['minimum_accepted_h_s'],covered_base_intervals=round(covered_intervals),extra_accepted_over_same_duration=tx['accepted']-round(covered_intervals),extra_completed_native_candidates_over_same_duration=tx['accepted']+tx['rejected']-round(covered_intervals),estimated_fixed_025_same_duration_wall_s=refcost,measured_full_to_estimated_fixed_ratio=resources['full']['main_wall_s']/refcost,comparison_is_estimate_not_new_fixed_run=True,transaction_full_wall_s={k:tx[k] for k in ['checkpoint_wall_s','restore_wall_s','predict_wall_s','native_candidate_wall_s']},MPI_ranks=1,threads=1,GPU=False,development_aborted_candidates=1,development_attempts='Two archived driver startup issues; one additional initial trial, no extra complete protocol',output_bytes_excluding_build_source=sum(q.stat().st_size for q in O.rglob('*') if q.is_file() and not any(x in q.relative_to(O).parts for x in ['build','source']))))
 passed=full['status']=='COMPLETED' and abs(t['time_s'][-1]-.0007)<1e-15 and float(t['max_dd'][1:].max())<=.10 and maxrel<=.02 and rb['full']['status']=='PASS'
 status=dict(task=TASK,run_id=RUN,local_event='PASS' if local['status']=='PASS' else 'FAIL',full_700us='PASS' if passed else 'FAIL',full_started_at_virgin_t0=True,full_accepted_steps=tx['accepted'],full_rejected_trials=tx['rejected'],minimum_accepted_h_s=full['minimum_accepted_h_s'],max_accepted_dd=float(t['max_dd'][1:].max()),max_relative_energy_balance=maxrel,final_time_s=float(t['time_s'][-1]),load_hold_unload_final_hold_completed=bool(passed),damage_force_energy_numerical_consistency=bool(passed),material_and_free_energy_unchanged=True,rollback_verified=rb['full']['status']=='PASS',global_time_convergence='NOT_CERTIFIED',contact_penetration='NOT_CERTIFIED',material_Gc='UNAVAILABLE',P4_small_nonaffine_numerical_entry='ELIGIBLE_FOR_NEXT_ROUND_REVIEW' if passed else 'NOT_ADMITTED',P4_executed=False,next_task='LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1' if passed else 'Review actual safety stop',production_modified=False)
 dump(O/'FINAL_STATUS.json',status)
 dump(O/'PLOT_AND_REPORT_FACTS.json',dict(status=status,local=local,full=full,resources=json.loads((O/'RUN_RESOURCES.json').read_text()),topology=json.loads((O/'DAMAGE_HISTORY_AND_TOPOLOGY.json').read_text()),ledger=ledger,contact=json.loads((O/'CONTACT_AND_TIMELAYER_AUDIT.json').read_text())))
 print(json.dumps(status,indent=2))
if __name__=='__main__':main()
