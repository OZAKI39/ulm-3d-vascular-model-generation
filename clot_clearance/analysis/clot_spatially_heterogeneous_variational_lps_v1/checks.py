from common import *
from native import native,predict
from spatial_heterogeneity_field import generate
import time,ctypes
start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);model=energy_model(ref);out={}
old=load_npz(O/'checks/old_uniform.npz');new=load_npz(O/'checks/new_uniform.npz')
for key in old:
 if key.endswith('_scalars'):assert np.array_equal(old[key][:6],new[key][:6])
 else:assert np.array_equal(old[key],new[key]),key
out['uniform_fallback']=dict(status='PASS',old_new_bitwise_force_energy_history=True,states=['elastic','near_onset','softened'],original_reference_m_exact=True)
a=load_npz(O/'HETEROGENEITY_FIELD.npz');b=generate(ref,np.arange(len(ref['X']))[::-1]);assert all(np.array_equal(a[k],b[k]) for k in a)
out['reproducibility']=dict(status='PASS',reverse_ID_order_bitwise=True)
h,C=native(O/'checks/heterogeneous','heterogeneous');assert np.array_equal(h.field(14),a['bond_sc']);assert np.array_equal(h.field(13),S0[ref['bond_types']]);assert np.array_equal(h.field(4),ref['m']);assert np.array_equal(h.field(9),a['bond_ids'])
# A changed ID order on input is already used by native(): no slot-based parameters.
fdrows=[];prevk=h.field(2);prevd=h.field(3);prevD=0
for token,strain in enumerate([.01,.15001,.163],1):
 x=ref['X']*(1+strain);before=h.static(x);assert np.array_equal(h.field(2),prevk) and np.array_equal(h.field(3),prevd)
 k,d=predict(model,ref,x,prevk,h.field(13),h.field(14));priorU=model.evaluate(x,prevd)['U']
 assert h.lib.lpsvar_action(h.l.lmp,2,token,None,None)==0;s=h.static();oracle=model.evaluate(x,s['d'])
 assert np.max(abs(s['d']-d))<1e-12 and np.max(abs(s['kappa']-k))<1e-12
 assert np.max(abs(s['material_force']-oracle['force']))<1e-12 and abs(s['scalars'][0]-oracle['U'])<1e-17
 assert np.all(s['d']>=prevd-1e-14) and np.all(s['kappa']>=prevk-1e-14)
 drop=priorU-oracle['U'];assert drop>=-1e-20 and abs(s['scalars'][1]-prevD-drop)<1e-17
 # Fixed irreversible damage when differentiating stored energy, not a total derivative through d(x).
 nodes=np.unique([np.argmin(a['q'][:n]),np.argmax(a['q'][:n]),n//2]);eps=2e-10;errs=[];forces=[]
 for node in nodes:
  for c in range(3):
   xp=x.copy();xm=x.copy();xp[node,c]+=eps;xm[node,c]-=eps
   fd=-(model.evaluate(xp,s['d'])['U']-model.evaluate(xm,s['d'])['U'])/(2*eps)
   errs.append(abs(fd-s['material_force'][node,c]));forces.append(abs(s['material_force'][node,c]))
 rel=max(errs)/max(max(forces),1e-20);assert rel<1e-5,(strain,rel)
 assert np.linalg.norm(s['material_force'].sum(0))<1e-15 and np.linalg.norm(np.cross(x,s['material_force']).sum(0))<1e-18
 fdrows.append(dict(strain=strain,selected_original_atom_ids=nodes+1,max_gradient_error_N=max(errs),relative_gradient_error=rel,oracle_force_error_N=float(np.max(abs(s['material_force']-oracle['force']))),oracle_U_error_J=float(abs(s['scalars'][0]-oracle['U'])),damage_drop_J=drop))
 prevk=s['kappa'];prevd=s['d'];prevD=s['scalars'][1]
# Rejected predictor / trial force evaluation cannot commit history.
xtrial=ref['X']*1.3;trial_k,trial_d=predict(model,ref,xtrial,prevk,h.field(13),h.field(14));h.static(xtrial)
assert np.array_equal(h.field(2),prevk) and np.array_equal(h.field(3),prevd) and h.field(8)[1]==prevD
h.static(x);before_tags=h.l.numpy.extract_atom('id')[:h.N].copy();h.l.commands_string('atom_modify sort 1 0.0001\nrun 1');ss=h.state();after_tags=h.l.numpy.extract_atom('id')[:h.N].copy()
assert np.array_equal(ss['d'],prevd) and np.array_equal(h.field(14),a['bond_sc'])
restart=O/'checks/softening.restart';h.l.command(f'write_restart {restart}');h.close()
r,C=native(O/'checks/restarted','heterogeneous',restart);sr=r.static()
for key in ['x','d','kappa']:assert np.array_equal(ss[key],sr[key]),key
assert np.array_equal(r.field(14),a['bond_sc']) and np.array_equal(r.field(4),ref['m']) and sr['scalars'][1]==prevD
r.static(ref['X']);assert r.lib.lpsvar_action(r.l.lmp,2,4,None,None)==0;un=r.static();assert np.array_equal(un['d'],prevd) and np.array_equal(un['kappa'],prevk) and un['scalars'][1]==prevD
r.close();out['variational']=fdrows;out['history']=dict(status='PASS',sorting_changed_slots=int(np.sum(before_tags!=after_tags)),restart_V4_sc_and_original_ID_history_exact=True,unloading_irreversible=True,rejected_trial_no_commit=True,nonnegative_damage_dissipation=True)
# Exactly one new-path native VV step per material mode with real nonzero drag.
rows=[]
for mode in CASES:
 ns,C=native(O/'checks'/('dynamic_'+mode),mode);s=ns.dynamics();ns.l.command(f'reset_timestep 0 time {HOLD:.17g}')
 vel=np.zeros_like(ref['X']);vel[:n]=[.02,.003,0];ids=ns.l.numpy.extract_atom('id')[:ns.N]-1;ns.l.numpy.extract_atom('v')[:ns.N]=vel[ids];ns.l.command('run 0');p=ns.state();h=MIN_H
 xp=p['x'].copy();xp[:n]+=h*(p['velocity'][:n]+h*p['force'][:n]/(2*RHO*V));k,d=predict(model,ref,xp,p['kappa'],ns.field(13),ns.field(14));s=ns.step(h)
 uc,fc,_=C(ref).evaluate(s['x'],s['d']);z=ns.tool.evaluate(s['x'],s['velocity'],s['d'],s['time']);w=ns.tool.work(p,s)
 res=w.sum()-(s['scalars'][0]-p['scalars'][0]+s['K']-p['K']+s['scalars'][1]-p['scalars'][1]+s['Uc']-p['Uc']);rel=abs(res)/max(abs(w.sum()),1e-30)
 assert np.max(abs(xp-s['x']))<2e-18 and np.max(abs(d-s['d']))<1e-12 and np.max(abs(s['force']-s['material_force']-fc-z['force']))<1e-12 and w[1]<0 and rel<.01
 rows.append(dict(mode=mode,steps=1,relative_energy_error=rel,work_drive_J=w[0],work_drag_J=w[1]));ns.close()
out['dynamic_energy']=rows
# Direct evidence from fixed d: sc does not enter U except through native irreversible d.
out['energy_interpretation']='Same positive endpoint-star energy U(x,d); pair sc enters irreversible damage update, D is the exact material-energy drop at fixed x; no force multiplier. Gradient holds d fixed.'
dump(O/'VARIATIONAL_ENERGY_AUDIT.json',dict(status='PASS',checks=out,wall_s=time.perf_counter()-start));dump(O/'KEY_CHECKS.json',dict(status='PASS',new_dynamic_microsteps=2,formal_cases_started=0,details='VARIATIONAL_ENERGY_AUDIT.json'));print('CHECKS_PASS',out,flush=True)
