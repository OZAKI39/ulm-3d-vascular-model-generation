"""Checks only this campaign's new load, evidence and reduction paths."""
from common import *
from prescribed_traction import Load,MODES,RISE,END
import csv

def main():
 ref=load_npz(O/'provenance/reference.npz');graph=module('validation_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py');checks={};freeze=json.loads((O/'BASELINE_AND_SPECIMEN_FREEZE.json').read_text());b=ref['bonds'];cut=ref['initial_cut'];n=int(ref['clot_count'])
 assert sha(LIB)==freeze['library_sha256']
 if (O/'STARTUP_RECOVERY_FREEZE.json').exists():
  recovery=json.loads((O/'STARTUP_RECOVERY_FREEZE.json').read_text())
  if (O/'PREDICTOR_IDENTITY_FIX.json').exists():
   fix=json.loads((O/'PREDICTOR_IDENTITY_FIX.json').read_text());assert sha(O/fix['original_source_archive'])==recovery['source_sha256']==fix['original_source_sha256'];assert sha(O/fix['original_source_archive'])==fix['original_source_sha256']
   if (O/'REPLAY_TIMER_FIELD_FIX.json').exists():
    timer=json.loads((O/'REPLAY_TIMER_FIELD_FIX.json').read_text());assert sha(O/timer['original_source_archive'])==fix['updated_source_sha256']==timer['original_source_sha256'];assert sha(A/'native_case_v2.py')==timer['updated_source_sha256']
   else:assert sha(A/'native_case_v2.py')==fix['updated_source_sha256']
  else:assert sha(A/'native_case_v2.py')==recovery['source_sha256']
 for name,expected in freeze['reference_hashes'].items():assert sha(O/'provenance'/name)==expected
 for name,expected in freeze['source_hashes'].items():assert sha(R/name)==expected,('Frozen solver changed',name)
 for p in (O/'cases').glob('*/SUMMARY.json'):
  z=json.loads(p.read_text());case=p.parent;h=list(csv.DictReader((case/'history.csv').open()));assert len(h)==z['accepted_steps'];state=load_npz(case/'keys'/('final.npz' if z['status']=='COMPLETED' else 'last_accepted.npz'));initial=load_npz(case/'keys/initial.npz')
  assert np.array_equal(initial['d'],cut.astype(float)) and float(initial['time'])==0 and initial['scalars'][1]==0
  assert np.array_equal(state['atom_ids'],np.arange(1,len(ref['X'])+1))
  assert np.array_equal(state['x'][n:],ref['X'][n:]) and not state['velocity'][n:].any()
  assert len(h)==int(state['step']) and abs(float(state['time'])-float(h[-1]['time_s']))<1e-15
  rel=max(float(q['relative_R']) for q in h);assert rel<=.01 and abs(rel-z['max_relative_energy_residual'])<1e-14
  err=max(abs(float(q['W_J'])-(float(q['U_J'])-initial['scalars'][0]+float(q['K_J'])-initial['K']+float(q['D_J'])+float(q['Uc_J'])-initial['Uc']-float(q['Jc_J']))-float(q['R_J'])) for q in h);assert err<1e-20
  assert max(float(q['max_dd']) for q in h)<=.1+1e-10
  assert all(float(y['time_s'])>float(x['time_s']) for x,y in zip(h,h[1:]))
  if z['status']=='COMPLETED':
   assert abs(float(state['time'])-END)<1e-14
   expected=z['params']['peak_pa']*z['loaded_area_m2']*(END-RISE/2)
   assert abs(z['last_row']['impulse_ns']-expected)<max(1e-20,expected*1e-8)
  lab,sizes,att=graph.topology(ref,state['d']);assert int(sizes[~att].sum())==z['last_row']['detached_nodes'];assert int(sizes[~att].max(initial=0))==z['last_row']['largest_detached_nodes']
  for ev,record in z['events'].items():
   st=load_npz(case/record['evidence']);assert abs(float(st['time'])-record['time_s'])<1e-15
   if 'original_bond_ids' in record:
    ids=np.array(record['original_bond_ids'])-1;assert not cut[ids].any();assert np.array_equal(b[ids]+1,record['atom_ids']);stretch=np.linalg.norm(st['x'][b[ids,1]]-st['x'][b[ids,0]],axis=1)/ref['r0'][ids]-1;assert np.max(abs(stretch-record['stretch']))<1e-13
    if ev=='first_failure':assert np.all(st['d'][ids]==1)
    if ev=='first_softening':assert np.all(st['d'][ids]>0)
  hh=json.loads((case/'DETACHED_COMPONENT_HISTORY.json').read_text());histmax=max((max((c['nodes'] for c in row['components']),default=0) for row in hh),default=0);assert histmax==z['history_max_free_nodes']
  if histmax:
   st=load_npz(case/'keys/max_free.npz');ll,ss,aa=graph.topology(ref,st['d']);assert int(ss[~aa].max(initial=0))==histmax
  for ev in z['secondary_fracture_events']:
   parent=set(ev['parent_atom_ids']);children=[set(c) for c in ev['children_atom_ids']];assert len(children)>1 and set.union(*children)==parent and sum(map(len,children))==len(parent)
  checks[case.name]=dict(status='PASS',native_window_status=z['status'],energy_reconstruction_error_j=err,identity_events_topology_and_impulse='PASS',mapping=Load(ref,z['params']).audit(ref))
 dump(O/'DELIVERY_VALIDATION.json',dict(status='PASS',cases=checks,scope='New load mapping, ledger and accepted evidence readback; no repeated historical full benchmark'))
 print('DELIVERY_VALIDATION_PASS',len(checks))
if __name__=='__main__':main()
