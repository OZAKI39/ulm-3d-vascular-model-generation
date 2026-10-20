"""Explicitly audit the inherited already-precut kappa=sc sentinel semantics."""
from common import *
ref=load_npz(O/'provenance/reference.npz');cut=ref['initial_cut'];states={m:load_npz(O/'cases'/m/'keys/initial.npz') for m in CASES};u=states['uniform'];h=states['heterogeneous'];diff=u['kappa']!=h['kappa']
assert np.array_equal(u['d'],h['d']) and np.array_equal(u['d'],cut.astype(float));assert not diff[~cut].any();assert not u['kappa'][~cut].any() and not h['kappa'][~cut].any()
for key in ['x','velocity','force','material_force','external_force']:assert np.array_equal(u[key],h[key]),key
assert u['scalars'][0]==h['scalars'][0]==0. and u['scalars'][1]==h['scalars'][1]==0.
model=energy_model(ref);test=ref['X'].copy();test[:,0]+=.03*test[:,2];a=model.evaluate(test,u['d']);b=model.evaluate(test,h['d']);assert np.array_equal(a['force'],b['force']) and a['U']==b['U']
# Read-only pre-onset common saved frames: actual states must be bitwise equal.
events={m:json.loads((O/'cases'/m/'FRACTURE_EVENTS.json').read_text()) for m in CASES};on=min(e['first_softening']['time_s'] for e in events.values());bytime={}
for m in CASES:
 bytime[m]={}
 for p in (O/'cases'/m/'keys').glob('frame_*.npz'):
  s=load_npz(p)
  if float(s['time'])<on:bytime[m][float(s['time'])]=s
common=sorted(set(bytime['uniform'])&set(bytime['heterogeneous']));assert common
for t in common:
 for k in ['x','velocity','force','d','material_force','external_force']:assert np.array_equal(bytime['uniform'][t][k],bytime['heterogeneous'][t][k]),(t,k)
dump(O/'INITIAL_HISTORY_AUDIT.json',dict(status='PASS_WITH_EXPLICIT_PRECUT_KAPPA_SENTINEL_DIFFERENCE',full_initial_kappa_bitwise_equal=False,changed_kappa_original_bond_ids=np.flatnonzero(diff)+1,changed_kappa_count=int(diff.sum()),all_changed_kappa_are_original_precut=bool(np.all(cut[diff])),initial_damage_exact=True,uncut_kappa_zero_exact=True,initial_x_v_force_energy_exact=True,pre_onset_common_saved_states_bitwise_equal=len(common),first_softening_s=on,explanation='Inherited initialization uses kappa=sc to encode d=1 for original precut. H uses its frozen local sc; these kappa markers are not equal to U and are not newly imposed physical strain. U(x,d) and native material forces depend on d, not kappa directly. All changed bonds were already cut and remain d=1. This is a documented initialization representation difference, not strict preservation of the full initial kappa array.'))
print('INITIAL_HISTORY_AUDIT_PASS',int(diff.sum()),len(common),flush=True)
