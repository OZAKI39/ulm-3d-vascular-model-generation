from native import *
import resource,importlib.util
spec=importlib.util.spec_from_file_location('frozen',R/'analysis/softening_feasibility_v1/model.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
ref=np.load(PREV/'provenance/reference.npz');m=mod.SofteningLPS(ref);X=m.X;start=time.perf_counter();records=[];freeze=json.loads((O/'NATIVE_SOFTENING_DESIGN_FREEZE.json').read_text())

def compare(native,off):
 fscale=float(np.linalg.norm(off['force']));ferr=float(np.linalg.norm(native['force']-off['force']));ua=abs(native['scalars'][0]-off['U']);fm=float(abs(native['force']-off['force']).max())
 return dict(force_max_abs_N=fm,force_relative_L2=ferr/fscale if fscale else None,force_scale_L2_N=fscale,energy_abs_J=ua,energy_relative=ua/abs(off['U']) if off['U'] else None,theta_max_abs=float(abs(native['theta']-off['theta']).max()),endpoint_v_max_abs_m=float(abs(native['v']-off['v']).max()),net_force_N=float(np.linalg.norm(native['force'].sum(0))),torque_Nm=float(np.linalg.norm(np.cross(native['x'],native['force']).sum(0))),native_tally_energy_abs_J=abs(native['ebond_J']-off['U']),pass_force=fm<freeze['force_parity_absolute_N'] and (not fscale or ferr/fscale<freeze['force_parity_relative']),pass_energy=ua<freeze['energy_parity_absolute_J'] and (not off['U'] or ua/abs(off['U'])<freeze['energy_parity_relative']))
ns=Native('undamaged_new');old=Native('undamaged_legacy',legacy=True)
for name,x in [('reference',X),('small_affine',np.load(PREV/'baseline_uniform_small.npz')['x']),('nonaffine',np.load(PREV/'baseline_shear_nonaffine.npz')['x'])]:
 st=ns.eval(x);lg=old.eval(x);off=m.evaluate(x);z=compare(st,off);z.update(state=name,legacy_force_max_abs_N=float(abs(st['force']-lg['force']).max()),legacy_force_relative_L2=float(np.linalg.norm(st['force']-lg['force'])/np.linalg.norm(lg['force'])) if np.linalg.norm(lg['force']) else None,legacy_energy_tally_not_reference=True)
 records.append(z);np.savez_compressed(O/(name+'_parity.npz'),**st,legacy_force_N=lg['force'],offline_force_N=off['force'],offline_U_J=off['U'])
ns.close();old.close();passed=all(z['pass_force'] and z['pass_energy'] and z['legacy_force_max_abs_N']<2e-17 for z in records);dump(O/'UNDAMAGED_NATIVE_PARITY.json',dict(status='PASS' if passed else 'FAIL',records=records));assert passed
# Historical binary topology comparison: candidate and old LPS have DIFFERENT relaxed damaged energies.
hist=np.load(PREV/'provenance/postbreak_same_coordinates.npz');d=(~hist['alive']).astype(float);ns=Native('binary_five_break');ns.import_state(d*.24,d);st=ns.eval(hist['x']);off=m.evaluate(hist['x'],d);diff=compare(st,off);diff.update(broken=5,old_vs_new_force_difference_N=float(np.linalg.norm(st['force']-hist['force_N'])),old_binary_definition_preserved_by_offline_baseline=float(np.linalg.norm(m.baseline(hist['x'],hist['alive'])['force']-hist['force_N'])),interpretation='Expected constitutive binary-limit difference: new endpoint-star free internal variables relax; original binary LPS deletes extensions. New candidate must match frozen softening, not be silently redefined to old deletion.')
np.savez_compressed(O/'binary_five_break.npz',**st,legacy_force_N=hist['force_N'],offline_force_N=off['force']);dump(O/'BINARY_LIMIT_DIFFERENCE.json',diff);ns.close();assert diff['pass_force'] and diff['pass_energy']
records=[];partial_states={}
for name,sc in [('soft_sc018',.18),('soft_sc024',.24)]:
 z=np.load(PREV/(name+'.npz'));it=int(np.flatnonzero((z['phase']==2)&(z['series'][:,0]>.16))[-1]);ns=Native(name+'_fixed',sc);ns.import_state(z['kappa'][it],z['d'][it]);st=ns.eval(z['x'][it]);off=m.evaluate(st['x'],st['d']);v=compare(st,off);v.update(case=name,frame=it,kappa_max_abs=float(abs(st['kappa']-z['kappa'][it]).max()),d_max_abs=float(abs(st['d']-z['d'][it]).max()));records.append(v);np.savez_compressed(O/(name+'_partial_parity.npz'),**st,offline_force_N=off['force'],offline_U_J=off['U']);partial_states[name]=st;ns.close()
passed=all(v['pass_force'] and v['pass_energy'] for v in records);dump(O/'PARTIAL_DAMAGE_FORCE_ENERGY_PARITY.json',dict(status='PASS' if passed else 'FAIL',records=records));assert passed
# Fixed nonuniform diagnostic; imported history uses inverse law only to remain admissible.
z=np.load(PREV/'fixed_damage_diagnostic.npz');d=z['d'];kappa=.24*.15/(.24-d*(.24-.15));ns=Native('fixed_gradient');ns.import_state(kappa,d);st=ns.eval(z['x']);off=m.evaluate(st['x'],d);check=compare(st,off);assert check['pass_force'] and check['pass_energy']
# Six complete endpoint-star dense solves using exact frozen matrices.
solves=[]
for a in [0,34,85,170,255,340]:
 qq=np.flatnonzero(m.center==a);b=m.edge[qq];li=m.l[qq];H=np.diag(m.Ddir[qq])+m.A*np.outer(li,li);pen=m.krefdir[qq]*(1-d[b])/d[b];vd=np.linalg.solve(H+np.diag(pen),pen*off['e'][b]);err=float(np.max(abs(st['v'][qq]-vd)));solves.append(dict(atom_id=a+1,degree=len(qq),dense_vs_native_v_abs_m=err))
dump(O/'LOCAL_SOLVER_NUMERICAL_AUDIT.json',dict(status='PASS' if max(v['dense_vs_native_v_abs_m'] for v in solves)<1e-16 else 'FAIL',records=solves,method='Exact rank-one elimination vs independent dense FP64 endpoint solves'))
assert max(v['dense_vs_native_v_abs_m'] for v in solves)<1e-16
rows=[];scale=float(abs(st['force']).max());dofs=np.r_[np.argsort(abs(st['force']).ravel())[-4:],[0,111,444,900]]
for h in freeze['gradient_steps_m']:
 for dof in dofs:
  xp=z['x'].copy();xm=z['x'].copy();xp.flat[dof]+=h;xm.flat[dof]-=h
  sp=ns.eval(xp);sm=ns.eval(xm);fd=-(sp['scalars'][0]-sm['scalars'][0])/(2*h);err=abs(fd-st['force'].flat[dof]);rows.append(dict(atom_id=int(dof//3+1),axis=int(dof%3),h_m=h,force_N=float(st['force'].flat[dof]),FD_N=fd,abs_N=err,relative_global=err/scale))
  assert np.array_equal(sp['d'],d) and np.array_equal(sm['d'],d) and sp['scalars'][4]==0 and sm['scalars'][4]==0
passed=max(v['relative_global'] for v in rows)<freeze['gradient_relative_global'] and max(v['abs_N'] for v in rows)<freeze['gradient_absolute_N']
dump(O/'FIXED_DAMAGE_GRADIENT_AUDIT.json',dict(status='PASS' if passed else 'FAIL',force_scale_N=scale,rows=rows,max_error=max(rows,key=lambda a:a['abs_N']),nonuniform_parity=check,damage_commits_during_trial=0));np.savez_compressed(O/'nonuniform_fixed_damage.npz',**st,offline_force_N=off['force']);ns.close();assert passed
resources=dict(wall_s=time.perf_counter()-start,CPU_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,native_instances=6,fixed_FD_evaluations=len(rows)*2)
dump(O/'CHECK_RESOURCES.json',resources);print('MECHANICAL_GATES_PASS',resources,flush=True)
