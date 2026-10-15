from native import *
from typed_model import TypedLPS
import resource
ref=dict(np.load(O/'provenance/reference.npz'));model=TypedLPS(ref);contact=Contact(ref);legacy=np.load(O/'preflight/legacy_force_parity.npz');b=ref['bonds'];bt=ref['bond_types'];N=len(ref['X']);n=int(ref['clot_count']);X=ref['X'];nb=len(b)
assert not (O/'MINIMAL_INTEGRATION_CHECK.json').exists()
started=time.perf_counter();ns=Native('preflight/parity');parity=[]
for label,x in [('initial',X),('perturbed',legacy['perturbed'])]:
 s=ns.static(x);pred=model.evaluate(x);uc,fc,pairs=contact.evaluate(x,s['d']);f=s['force'];old=legacy[label+'_force'];scale=float(np.max(abs(old)))
 err=float(abs(f-old).max());errref=float(abs(s['material_force']-pred['force']).max());errm=float(np.max(abs(ns.field(4)-legacy[label+'_m'])/ref['m']))
 assert err<2e-17 and errref<2e-17 and errm<2e-12,(label,err,errref,errm)
 assert np.array_equal(ns.field(9),np.arange(1,nb+1)) and np.array_equal(ns.field(12),bt)
 assert np.array_equal(ns.field(13),S0[bt]) and np.array_equal(ns.field(14),SC[bt])
 assert abs(s['scalars'][0]-pred['U'])<1e-20
 rec=dict(state=label,max_node_force_delta_N=err,max_material_node_delta_N=float(abs(f[:n]-old[:n]).max()),max_anchor_reaction_delta_N=float(abs(f[n:]-old[n:]).max()),independent_force_delta_N=errref,force_scale_N=scale,relative_force_delta=None if scale<1e-20 else err/scale,m_max_relative=errm,net_force_delta_N=np.sum(f-old,axis=0),torque_delta_Nm=np.sum(np.cross(x-x.mean(0),f-old),axis=0),native_U_J=s['scalars'][0],independent_U_J=pred['U'],native_pair_contact_N=float(abs(fc).max()))
 parity.append(rec)
 save_state(O/'preflight'/f'{label}.npz',s,old_force=old,independent_force=pred['force'])
ns.close();dump(O/'INITIAL_FORCE_PARITY.json',dict(status='PASS',physical_steps=0,states=parity,all_original_atoms_preserved=True,all_three_types_and_thresholds_exact=True))
# One localized real-bond commit/unload/reload check. These are static diagnostic
# coordinates only, excluded from the sole formal trajectory.
ib=int(np.flatnonzero(bt==1)[0]);jb=int(np.flatnonzero(bt==2)[-1]);a,z=b[ib];c,d=b[jb];assert c<n<=d
assert len(set([int(a),int(z),int(c),int(d)]))==4
ns=Native('preflight/mixed');s=ns.static();lastk=s['kappa'];lastd=s['d'];events=0;stages=[]
for token,stretch in enumerate([.16,0.,.205,0.],1):
 x=X.copy();x[z]=X[z]+stretch*(X[z]-X[a]);x[c]=X[c]+stretch*(X[c]-X[d]);ns.static(x)
 assert ns.lib.lpsvar_action(ns.l.lmp,2,token,None,None)==0
 s=ns.static();kk,dd=model.update(x,lastk);pred=model.evaluate(x,s['d'])
 assert np.max(abs(kk-s['kappa']))<2e-13 and np.max(abs(dd-s['d']))<2e-13
 assert not np.any(s['d']<lastd-2e-13) and not np.any(s['kappa']<lastk-2e-13)
 assert np.array_equal(s['x'][n:],X[n:]) and not s['velocity'][n:].any()
 assert s['failed_tension_max_N']<2e-17 and abs(pred['force']-s['material_force']).max()<2e-17
 events+=int(np.sum((s['d']==1)&(lastd<1)));assert s['scalars'][3]==events
 if token==1:assert 0<s['d'][ib]<s['d'][jb]<1
 if token==3:assert s['d'][ib]==s['d'][jb]==1
 if token in [2,4]:assert np.array_equal(s['d'],lastd)
 stages.append(dict(token=token,local_imposed_stretch=stretch,bulk_bond_id=ib+1,interface_bond_id=jb+1,bulk_kappa=s['kappa'][ib],interface_kappa=s['kappa'][jb],bulk_d=s['d'][ib],interface_d=s['d'][jb],total_events=events,native_bonds=s['native_bonds'],failed_tension_N=s['failed_tension_max_N']))
 lastk=s['kappa'].copy();lastd=s['d'].copy()
save_state(O/'preflight/mixed_static_terminal.npz',s);ns.close()
# Native full failure + contact exclusion + immediate detached drag test.
# Imported damage is explicit diagnostic preparation, never a natural event.
ns=Native('preflight/detached');ids=np.flatnonzero(np.bincount(b[bt==2].ravel(),minlength=N)[:n]>0);node=int(ids[len(ids)//2]);cut=np.any(b==node,axis=1);dk=cut.astype(float);kk=SC[bt]*dk
assert ns.lib.lpsvar_action(ns.l.lmp,1,0,ns.ptr(kk),ns.ptr(dk))==0
s=ns.static();assert s['native_bonds']==nb-int(cut.sum()) and s['failed_tension_max_N']==0
ns.dynamics();s=ns.state();assert s['detached'][node] and np.count_nonzero(s['detached'])==1
assert s['valid_flow'][node] and np.linalg.norm(s['drag'][node])>0 and not s['numeric'][node].any()
assert np.all(s['drag'][~s['detached']]==0) and len(s['drag'])==n
assert np.array_equal(s['x'][n:],X[n:]) and not s['velocity'][n:].any()
save_state(O/'preflight/imported_detachment.npz',s)
detcheck=dict(status='PASS',explicitly_imported_diagnostic_failure=True,formal_data=False,isolated_material_id=node+1,cut_original_bond_ids=np.flatnonzero(cut)+1,drag_active_same_native_state=True,drag_N=s['drag'][node],all_anchor_velocities_zero=True,anchor_drag_not_allocated=True)
ns.close()
# Short intact V3 dynamics and exact restart are part of this one integration suite.
ns=Native('preflight/dynamic');initial=ns.dynamics();p=initial;times=[];max_x=max_v=max_force=max_energy=0.
for k in range(8):
 ts=time.perf_counter();ns.l.command('run 1 pre no post no');s=ns.state();dx=s['x'][:n]-p['x'][:n];dv=s['velocity'][:n]-p['velocity'][:n]
 ex=float(abs(dx-BASE*p['velocity'][:n]-BASE**2/(2*MASS)*p['force'][:n]).max());ev=float(abs(dv-BASE/(2*MASS)*(p['force'][:n]+s['force'][:n])).max());max_x=max(max_x,ex);max_v=max(max_v,ev)
 pred=model.evaluate(s['x'],s['d']);max_force=max(max_force,float(abs(pred['force']-s['material_force']).max()));max_energy=max(max_energy,float(abs(pred['U']-s['scalars'][0])))
 assert abs(float(s['external_time'])-s['time'])<1e-15 and ex<2e-18 and ev<2e-12
 assert np.array_equal(s['x'][n:],X[n:]) and not s['velocity'][n:].any()
 ns.ext.accept(s);p=s;times.append(time.perf_counter()-ts)
cp=O/'preflight/dynamic.restart';ns.l.command(f'write_restart {cp}');save_state(O/'preflight/dynamic.state.npz',s);query=ns.ext.query;ns.close()
rs=Native('preflight/restart',restart=cp,restart_state=s,query=query);rr=rs.dynamics()
exact={key:float(np.max(abs(s[key]-rr[key]))) for key in ['x','velocity','kappa','d','force','material_force','scalars'] if key!='scalars'}
assert all(exact[k]==0 for k in ['x','velocity','kappa','d']) and exact['force']<2e-17 and exact['material_force']<2e-17,exact
assert np.max(abs(s['scalars'][:6]-rr['scalars'][:6]))<1e-20
assert np.array_equal(rs.field(12),bt) and np.array_equal(rs.field(13),S0[bt]) and np.array_equal(rs.field(14),SC[bt])
rs.close();q=resource.getrusage(resource.RUSAGE_SELF)
dump(O/'MINIMAL_INTEGRATION_CHECK.json',dict(status='PASS',mixed_static_stages=stages,detachment=detcheck,dynamic_steps=8,actual_dynamic_end_s=float(s['time']),max_VV_x_error_m=max_x,max_VV_v_error_m_s=max_v,max_independent_force_error_N=max_force,max_independent_energy_error_J=max_energy,typed_restart_max_errors=exact,anchors_fixed=True,median_short_step_wall_s=float(np.median(times)),nominal_protocol_steps=round(END/BASE),estimated_native_and_one_reference_wall_s=float(np.median(times)*END/BASE),peak_RSS_bytes=q.ru_maxrss*1024,wall_s=time.perf_counter()-started,formal_cases_run=0))
print('MINIMAL_V3_INTEGRATION_PASS',np.median(times),q.ru_maxrss*1024,flush=True)
