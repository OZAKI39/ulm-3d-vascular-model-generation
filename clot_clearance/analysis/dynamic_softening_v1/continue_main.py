from dynamic_run import *
old=json.loads((O/'dynamic341_AUDIT.json').read_text());assert old['accepted_steps']==170 and old['stop_reasons']==['EXCEPTION: AssertionError()']
assert not (O/'dynamic341_continuation_AUDIT.json').exists()
z=np.load(O/'dynamic341_trajectory.npz');ids=np.load(O/'provenance/341_boundary_ids.npz')
# Check the persisted pre-trial accepted state against the restored final restart
ns=Native('resume_identity',restart=O/'dynamic341_final.restart');ns.dynamics(ids['left_ids'],ids['right_ids'],ids['free_ids']);s=ns.state()
checks={k:bool(np.array_equal(z[k][-1],s[k])) for k in ['x','velocity','kappa','d']};checks['history_counters']=bool(np.array_equal(z['scalars'][-1,:6],s['scalars'][:6]));checks['force_max_abs_N']=float(abs(z['force'][-1]-s['force']).max());assert all(checks[k] for k in ['x','velocity','kappa','d','history_counters']) and checks['force_max_abs_N']<2e-17;ns.close()
dump(O/'INSTRUMENTATION_RECOVERY.json',dict(reason='Trial x must change U; only irreversible history and D/event/commit/token are invariant. Corrected Python assertion. Native C++ unchanged.',same_frozen_case=True,continued_from_accepted_step=170,checks=checks,source_before_trial_history_retained=True))
new=run_case('dynamic341_continuation',O/'provenance',ids['left_ids'],ids['right_ids'],ids['free_ids'],1230,restart=O/'dynamic341_final.restart')
dump(O/'DYNAMIC_341_NODE_AUDIT_SEGMENTS.json',dict(segments=[old,new],reason='Single pre-registered physical trajectory continued from unchanged accepted state after instrumentation repair.'))
