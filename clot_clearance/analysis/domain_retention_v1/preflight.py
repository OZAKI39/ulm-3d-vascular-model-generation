from event_integrator import *
from journal import Journal
import traceback
def main():
 ids=np.load(O/'provenance/341_boundary_ids.npz');args=[ids['left_ids'],ids['right_ids'],ids['free_ids']]
 with (O/'preflight/STARTED.json').open('x') as f:json.dump(dict(purpose='one initial parity check and one checkpoint fault-output test'),f)
 old=Native('preflight/old_box_initial',dt=BASE,box=None);old.dynamics(*args);s0=old.state();top0=topology(old);m0=old.field(4);mass0=old.l.extract_atom('mass')[1];old.close()
 new=Native('preflight/new_box_initial',dt=BASE);new.dynamics(*args);s1=new.state();top1=topology(new)
 parity=restore_comparison(s0,s1);assert parity['pass_']
 assert np.array_equal(s0['atom_ids'],np.arange(1,342)) and np.array_equal(s1['atom_ids'],s0['atom_ids'])
 assert np.array_equal(new.ref['X'],s1['x']) and np.array_equal(m0,new.field(4)) and mass0==new.l.extract_atom('mass')[1]==RHO*V
 assert np.array_equal(top0[0],top1[0]) and top0[1]==top1[1] and len(top1[0])==12047
 assert s1['Uc']==0 and np.all(s1['box_lo']==-BOX) and np.all(s1['box_hi']==BOX)
 save_state(O/'preflight/initial_old_box.npz',s0);save_state(O/'preflight/initial_new_box.npz',s1)
 dump(O/'INITIAL_STATE_PARITY.json',dict(status='PASS',comparison=parity,atom_ids_exact=True,reference_bonds_exact=True,m_exact=True,r0_taper_V_reference_file_sha_equal=sha(O/'provenance/reference.npz')==sha(PREV/'provenance/reference.npz'),mass_kg=mass0,coordinates_not_scaled_or_translated=True,old_box=[s0['box_lo'],s0['box_hi']],new_box=[s1['box_lo'],s1['box_hi']],initial_Uc_J=s1['Uc'],additional_walls_or_periodicity=False,no_physical_steps=True))
 new.close()
 # ONE targeted fault test from a genuine, nontrivially damaged old checkpoint.
 cp=sorted((PREV/'full/transactions').glob('*.restart'))[-1]
 ns=Native('fault_test/native_checkpoint',restart=cp,dt=BASE);ns.dynamics(*args);s=ns.state();topbefore=topology(ns)
 journal=Journal(O/'fault_test',['step','time_s']);journal.accept(s,[s['step'],s['time']]);journal.checkpoint(ns,s,'fault_test_valid_accept')
 journal.trial_start(1,s,BASE,0,float(np.min(BOX-abs(s['x']))))
 try:raise RuntimeError('INJECTED_BINDING_EXCEPTION_BEFORE_NATIVE_ADVANCE')
 except RuntimeError as exc:record=journal.failure(exc,traceback.format_exc())
 journal.finish();ns.close()
 latest=json.loads((O/'fault_test/LATEST_NATIVE_CHECKPOINT.json').read_text());restored=Native('fault_test/readback',restart=latest['restart'],dt=BASE);restored.dynamics(*args);r=restored.state();cmp=restore_comparison(s,r);topafter=topology(restored);restored.close()
 assert cmp['pass_'] and np.array_equal(topbefore[0],topafter[0]) and topbefore[1]==topafter[1]
 events=[json.loads(v) for v in (O/'fault_test/events.jsonl').read_text().splitlines()];assert [x['kind'] for x in events]==['ACCEPTED_STATE','ACCEPTED_CHECKPOINT','FAILURE']
 assert record['last_accepted_time_s']==s['time'] and record['attempt']['candidate_time_s']>s['time'] and record['native_calls_after_exception']==0
 last=np.load(O/'fault_test/last_accepted.npz');assert np.array_equal(last['x'],s['x']) and np.array_equal(last['kappa'],s['kappa'])
 dump(O/'FAILURE_SAFE_OUTPUT_AUDIT.json',dict(status='PASS',test_count=1,existing_checkpoint=str(cp),existing_checkpoint_sha256=sha(cp),no_physical_steps_executed=True,injection='Python binding exception before native advance; no deliberate atom loss',last_accepted_readback=cmp,restart_complete_history_and_topology=True,special_exclusions_exact=True,accepted_vs_trial_distinct=True,events_flushed=True,native_calls_in_failure_handler=0,failure_report='fault_test/FAILURE_RECORD.json',scope='Handled Python/LAMMPS binding exceptions; hard process kill can retain only last periodic native checkpoint and already flushed logs. Not a full fault-injection matrix.'))
 d=json.loads((O/'DOMAIN_CHANGE_AUDIT.json').read_text());d['status']='INITIAL_PARITY_PASS_DOMAIN_FROZEN';dump(O/'DOMAIN_CHANGE_AUDIT.json',d)
 print('P0_P1_PASS',s['time'],cmp,flush=True)
if __name__=='__main__':main()
