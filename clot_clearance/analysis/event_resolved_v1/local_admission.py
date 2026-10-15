from native import *
from event_integrator import restore_comparison
audit=json.loads((O/'local/DYNAMIC_AUDIT.json').read_text());tx=json.loads((O/'local/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());rb=json.loads((O/'local/ROLLBACK_STATE_CONSISTENCY.json').read_text())
z=np.load(O/'local/trajectory.npz');cols=list(z['columns']);arr=z['series'];ci={c:i for i,c in enumerate(cols)}
rejected=[t for t in tx['trials'] if t['decision']=='REJECT'];first=rejected[0] if rejected else None
matches=[]
if first:
 candidate=np.load(O/'local/transactions'/f'rejected_trial_{first["trial"]:05d}.npz')
 previous=np.load(PREV/'keyframes/step_002562.npz')
 matches=restore_comparison(previous,candidate)
assert first and matches['pass_'] and first['max_bond_id']==1813 and first['max_dd']>.10
at=np.flatnonzero(abs(arr[:,ci['time_s']]-.0006405)<1e-15);assert len(at)==1
j=at[0];old=np.load(PREV/'trajectory.npz')['series'];newdelta=arr[j]-arr[0];olddelta=old[2562]-old[2561]
comparison={k:dict(old_single_step=float(olddelta[ci[k]]),new_two_steps=float(newdelta[ci[k]])) for k in ['W_J','U_J','K_J','D_J','Uc_J','R_J','VV_defect_J','spatial_defect_J','switch_work_J','contact_defect_J']}
statepath=next(t for t in tx['trials'] if t['decision']=='ACCEPT' and abs(t['end_s']-.0006405)<1e-15)
snew=np.load(O/'local/transactions'/f'accepted_trial_{statepath["trial"]:05d}.npz')
comparison['endpoint_state_difference_640p5']=dict(max_abs_x_m=float(abs(snew['x']-candidate['x']).max()),max_abs_v_m_s=float(abs(snew['velocity']-candidate['velocity']).max()),max_abs_force_N=float(abs(snew['force']-candidate['force']).max()),kappa1813=float(snew['kappa'][1812]),d1813=float(snew['d'][1812]),failed_old=int(np.sum(candidate['d']==1)),failed_new=int(np.sum(snew['d']==1)))
passed=audit['status']=='COMPLETED' and abs(audit['final_time_s']-.000643)<1e-15 and rb['status']=='PASS' and bool(rb['records']) and audit['maxima']['max_dd']<=.10
out=dict(status='PASS' if passed else 'FAIL',actual_window_s=[float(arr[0,1]),float(arr[-1,1])],old_violating_candidate_reproduced=matches,first_rejected_trial=first,accepted_substeps=tx['accepted'],rejected_trials=tx['rejected'],minimum_accepted_h_s=audit['minimum_accepted_h_s'],max_accepted_dd=float(arr[1:,ci['max_dd']].max()),full_rollback_verified=rb['status']=='PASS',comparison_same_640p25_to_640p5_window=comparison,ledger_origin='Old immutable prefix to 640.25us, then new local accepted substeps. Full run must begin at new t0.',no_global_time_convergence_claim=True)
dump(O/'LOCAL_BOND1813_EVENT_AUDIT.json',out);print(json.dumps(out,default=lambda x:x.item() if isinstance(x,np.generic) else x,indent=2))
assert passed,'LOCAL_FAILED_NO_FULL_RUN'
