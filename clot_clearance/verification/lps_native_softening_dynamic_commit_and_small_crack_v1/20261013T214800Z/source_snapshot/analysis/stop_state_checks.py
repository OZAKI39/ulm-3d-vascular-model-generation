from native import *
ref=np.load(O/'provenance/reference.npz');ids=np.load(O/'provenance/341_boundary_ids.npz');records=[]
for name,rst in [('accepted170','dynamic341_final.restart'),('stopped244','dynamic341_continuation_final.restart'),('half_window278','dt_half_window_final.restart')]:
 ns=Native('readonly_'+name,restart=O/rst,dt=.25e-6 if 'half' in name else .5e-6);ns.dynamics(ids['left_ids'],ids['right_ids'],ids['free_ids']);s=ns.state();nb,raw=ns.l.gather_bonds();arr=np.ctypeslib.as_array(raw,shape=(3*nb,)).reshape(-1,3).copy();actual=sorted(tuple(sorted(row[1:])) for row in arr);expected=sorted(tuple(row+1) for row in ref['bonds'][s['d']<1]);same=actual==expected;assert same and np.all(arr[:,0]==1)
 merror=float(abs(ns.field(4)-ref['m']).max());assert merror==0;check=dict(name=name,native_unique_bonds=nb,all_original_live_endpoint_pairs_match=same,fixed_m_abs=merror,accepted_commit_count=int(s['scalars'][4]),physical_step=s['step'])
 if name=='stopped244':
  x=s['x'].copy();components=[int(np.argmax(abs(s['material_force']).ravel())),3*(int(ids['free_ids'][len(ids['free_ids'])//2])-1)+1,3*(int(ids['free_ids'][0])-1)+2];eps=1e-10;rr=[]
  for j in components:
   energies=[]
   for sign in [-1,1]:
    xx=x.copy();xx.ravel()[j]+=sign*eps;ns.l.numpy.extract_atom('x')[:ns.N]=xx;ns.l.command('run 0');st=ns.state();assert np.array_equal(s['d'],st['d']) and np.array_equal(s['kappa'],st['kappa']) and np.array_equal(s['scalars'][1:6],st['scalars'][1:6]);energies.append(st['scalars'][0]+st['Uc'])
   gradient=-(energies[1]-energies[0])/(2*eps);err=abs(gradient-s['force'].ravel()[j]);rr.append(dict(node_id=j//3+1,axis=j%3,minus_energy_gradient_N=gradient,native_force_N=float(s['force'].ravel()[j]),abs_error_N=err,relative_global_force=err/abs(s['force']).max()))
  ns.l.numpy.extract_atom('x')[:ns.N]=x;ns.l.command('run 0');restored=ns.state();assert np.array_equal(restored['d'],s['d']) and abs(restored['force']-s['force']).max()<2e-17;check['fixed_history_gradient']=rr;check['gradient_max_relative']=max(r['relative_global_force'] for r in rr);assert check['gradient_max_relative']<2e-7
 records.append(check);ns.close()
dump(O/'STOP_STATE_TOPOLOGY_GRADIENT_AUDIT.json',dict(status='PASS_FROZEN_READONLY_CHECKS',records=records,no_dynamic_steps=True))
print('STOP_STATE_CHECKS',records,flush=True)
