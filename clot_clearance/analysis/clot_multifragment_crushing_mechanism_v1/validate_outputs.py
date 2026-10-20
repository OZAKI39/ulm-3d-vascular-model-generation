"""Read back this round's accepted evidence; never rerun historical native tests."""
from common import *
from analyze import read,histories,graph
from indenter import *
import gzip
def main():
 ref=load_npz(O/'provenance/reference.npz');freeze=read(O/'SPECIMEN_AND_MATERIAL_FREEZE.json');checks=read(O/'KEY_CONSISTENCY_CHECKS.json');assert checks['status']=='PASS';assert sha(OLDLIB)==freeze['baseline_solver_sha256'];assert sha(LIB)==read(O/'BUILD_AUDIT.json')['library_sha256']
 for name,want in freeze['reference_sha256'].items():assert sha(O/'provenance'/name)==want
 assert sha(O/'HETEROGENEITY_FIELD.npz')==read(O/'HETEROGENEITY_FREEZE.json')['field_sha256']
 result={};n=int(ref['clot_count']);cut=ref['initial_cut'];b=ref['bonds'];bt=ref['bond_types'];field=load_npz(O/'HETEROGENEITY_FIELD.npz');tool=Indenter(ref)
 for name in ['uniform','heterogeneous']:
  case=O/'cases'/name;s=read(case/'SUMMARY.json');request=read(case/'REQUEST.json');rows=histories(case);assert s['status']=='COMPLETED' and len(rows)==s['accepted_steps'];assert abs(rows[-1]['time_s']-END)<1e-14
  for filename in ['common.py','indenter.py','native.py','run_case.py']:assert sha(A/filename)==request['source_sha256'][filename],filename
  assert request['library_sha256']==sha(LIB)
  for filename,want in request['input_sha256'].items():assert sha(O/filename)==want,filename
  initial=load_npz(case/'keys/initial.npz');assert np.array_equal(initial['d'],cut.astype(float));energyerror=max(abs(q['W_J']-(q['U_J']-initial['scalars'][0]+q['K_J']-initial['K']+q['D_J']+q['Uc_J']-initial['Uc']-q['Jc_J']+q['U_indenter_J']-initial['indenter_U'])-q['R_J']) for q in rows);assert energyerror<1e-20
  assert max(q['relative_R'] for q in rows)<=.01 and max(q['max_dd'] for q in rows)<=.1+1e-10
  assert all(q['numerical_damping_work_J']==0 for q in rows)
  onset=field['s0_by_original_bond_id'] if name=='heterogeneous' else S0[bt];keys={};native_frames=0
  for p in (case/'keys').glob('*.npz'):
   if 'failed' in p.name:continue
   state=load_npz(p);assert np.array_equal(state['atom_ids'],np.arange(1,len(ref['X'])+1));assert np.array_equal(state['x'][n:],ref['X'][n:]) and not state['velocity'][n:].any();assert np.all(state['d'][cut]==1);keys.setdefault(int(state['step']),[]).append(p)
   lab,sizes,att=graph.topology(ref,state['d']);free=sizes[~att]
   if state['step']:
    row=rows[int(state['step'])-1];assert int(free.sum())==row['detached_nodes'] and int(free.max(initial=0))==row['largest_nodes'] and int(np.sum((free>=20)&(free<=200)))==row['N20_200']
   for label in np.flatnonzero(~att):
    ids=np.flatnonzero(lab==label);boundary,_=graph.boundary(ref,ids);assert not np.any(state['d'][boundary]<1)
   oracle=tool.evaluate(state['x'],float(state['time']));assert np.max(abs(oracle['force']-state['indenter_force']))<1e-15 and abs(oracle['U']-float(state['indenter_U']))<1e-17
   assert np.max(abs(state['indenter_reaction']+oracle['force'].sum(0)))<1e-12
   assert not state['indenter_force'][:n][oracle['overlap']==0].any();native_frames+=1
  # Replay only sparse damage records, with original IDs, and compare saved FP64.
  damage=cut.astype(float);pending={k:v for k,v in keys.items() if k>0};laststep=0
  with gzip.open(case/'damage_events.jsonl.gz','rt') as f:
   for line in f:
    z=json.loads(line);step=z['step'];ids=np.array(z['original_bond_ids'])-1;values=np.array(z['d']);assert step>laststep and not cut[ids].any() and np.all(values>=damage[ids]);
    for k in sorted([k for k in pending if k<step]):
     for p in pending.pop(k):assert np.array_equal(load_npz(p)['d'],damage),p
    damage[ids]=values;laststep=step
    for p in pending.pop(step,[]):assert np.array_equal(load_npz(p)['d'],damage),p
  for paths in pending.values():
   for p in paths:assert np.array_equal(load_npz(p)['d'],damage),p
  registry={int(k):v for k,v in read(case/'COMPONENT_REGISTRY.json').items()};count=0
  with gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','rt') as f:
   for line in f:
    z=json.loads(line);ids=[]
    for q in z['components']:ids.extend(registry[q['component_id']]['atom_ids'])
    assert len(ids)==len(set(ids))==z['detached_nodes'];assert z['attached_nodes']+z['detached_nodes']==n
    if z['step']:assert z['detached_nodes']==rows[z['step']-1]['detached_nodes']
    count+=1
  assert count==len(rows)+1
  for ev in s['secondary_events']:
   parent=set(registry[ev['parent_component_id']]['atom_ids']);children=[set(registry[q]['atom_ids']) for q in ev['children_component_ids']];assert set.union(*children)==parent and sum(map(len,children))==len(parent)
  for nameev in ['first_softening','first_failure']:
   ev=s['events'].get(nameev)
   if ev:
    state=load_npz(case/ev['evidence']);ids=np.array(ev['original_bond_ids'])-1;assert not cut[ids].any();assert np.all(state['d'][ids]>0) and np.all(state['kappa'][ids]>onset[ids]);assert nameev!='first_failure' or np.all(state['d'][ids]==1)
  result[name]=dict(status='PASS',accepted_records=len(rows),component_history_records=count,accepted_FP64_states=native_frames,energy_reconstruction_error_J=energyerror,original_ID_damage_replay='PASS',free_component_anchor_boundary='PASS',native_current_contact='PASS',driver_and_frozen_inputs_unchanged=True)
 dump(O/'DELIVERY_VALIDATION.json',dict(status='PASS',cases=result,scope='New interface, original-ID field, accepted evidence and bookkeeping only; no additional physics runs'));print('DELIVERY_VALIDATION_PASS',result)
if __name__=='__main__':main()
