"""Check mixed-DOF work ledger on the already computed last preflight step."""
from mechanics import *
ids=np.load(O/'provenance/boundary_ids.npz');g=json.loads((PREV/'DOMAIN_PROTOCOL_FREEZE.json').read_text())['gates']
ns=Native('preflight/ledger_readback',restart=O/'preflight/mixed_dof.restart',dt=BASE/2);ns.dynamics(ids['bottom_ids'],ids['top_ids'],ids['free_ids']);p=ns.state();s=dict(np.load(O/'preflight/last.npz'));mech=Mechanics(ns,ids,p,g);row,pairs,soft,fail=mech.check(p,s,BASE)
assert not len(soft) and not len(fail)
dump(O/'MIXED_DOF_LEDGER_CHECK.json',dict(status='PASS',no_new_physical_steps=True,source='Existing preflight step 2 checkpoint and saved step 3; no replay',start_time_s=p['time'],end_time_s=float(s['time']),errors=mech.maxerr,ledger=row));ns.close();print('MIXED_DOF_LEDGER_PASS',row['relative_R'],row['momentum_residual_Ns'])
