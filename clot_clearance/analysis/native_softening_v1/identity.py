from native import *
import resource
assert json.loads((O/'FULL_DAMAGE_COMPONENT_INDEPENDENCE.json').read_text())['status']=='PASS'
start=time.perf_counter();mode=sys.argv[1];z=np.load(PREV/'planar_partial_sc024.npz');idx=int(np.argmax(z['series'][:,-1]));x=z['x'][idx];d=z['d'][idx];k=z['kappa'][idx];ref=np.load(PREV/'provenance/reference.npz');edges=ref['bonds']
ns=Native(mode+'_identity');ns.import_state(k,d);before=ns.eval(x);owners_before=COMM.allgather(ns.owned());ghost_before=COMM.allgather(ns.l.extract_global('nghost'))
if mode=='sort':
 assert COMM.size==1
 ns.l.command(f'atom_modify sort 1 {DX*.65:.17g}');after=ns.eval();owners_after=COMM.allgather(ns.owned());changed=int(np.sum(owners_before[0]!=owners_after[0]));triggered=changed>0
else:
 assert COMM.size==2
 ns.l.command('balance 1.0 x 0.42');after=ns.eval();owners_after=COMM.allgather(ns.owned());ob=np.zeros(341,int);oa=ob.copy()
 for rankid,ids in enumerate(owners_before):ob[ids-1]=rankid
 for rankid,ids in enumerate(owners_after):oa[ids-1]=rankid
 changed=int(np.sum(ob!=oa));triggered=changed>0
# Native explicit commit across ranks at a real additional admissible opening.
y=x.copy();y[:,0]+=np.where(ref['X'][:,0]<-30e-6,-1.,1.)*.15*DX
ns.eval(y);accepted=ns.commit(0);repeat=ns.commit(0)
# One-rank frozen offline reference, evaluated directly without rerunning a path.
import importlib.util
sp=importlib.util.spec_from_file_location('frozen',R/'analysis/softening_feasibility_v1/model.py');mm=importlib.util.module_from_spec(sp);sp.loader.exec_module(mm);m=mm.SofteningLPS(ref);kn,dn=m.update(y,k,.15,.24);off=m.evaluate(y,dn);Dref=m.evaluate(y,d)['U']-off['U']
allk=COMM.allgather(accepted['kappa']);alld=COMM.allgather(accepted['d']);allstats=COMM.allgather(accepted['scalars'][:6]);safety=max(float(abs(t-allk[0]).max()) for t in allk)==0 and max(float(abs(t-alld[0]).max()) for t in alld)==0 and all(np.array_equal(t,allstats[0]) for t in allstats)
result=dict(status='PASS' if triggered else ('MPI_MIGRATION_NOT_EXERCISED' if mode=='mpi' else 'SORT_NOT_EXERCISED'),source_offline_frame=idx,original_bond_ids_unchanged=bool(np.array_equal(before['bond_ids'],after['bond_ids'])),changed_owned_slots_or_atoms=changed,actual_owners_before=owners_before,actual_owners_after=owners_after,ghost_counts_before=ghost_before,ghost_counts_after=COMM.allgather(ns.l.extract_global('nghost')),damage_sort_or_migration_error=float(abs(after['d']-before['d']).max()),kappa_error=float(abs(after['kappa']-before['kappa']).max()),reference_m_error=float(abs(after['m']-before['m']).max()),force_error_N=float(abs(after['force']-before['force']).max()),energy_error_J=abs(after['scalars'][0]-before['scalars'][0]),commit_all_ranks_identical=safety,commit_native_offline_force_max_N=float(abs(accepted['force']-off['force']).max()),commit_native_offline_damage_max=float(abs(accepted['d']-dn).max()),commit_native_offline_D_abs_J=abs(accepted['scalars'][1]-Dref),accepted_global_commits=int(accepted['scalars'][4]),actual_new_failure_events=int(accepted['scalars'][3]),expected_new_failure_events=int(np.sum((dn==1)&(d<1))),repeat_did_not_increment=bool(np.array_equal(accepted['scalars'][:6],repeat['scalars'][:6])),native_bonds=accepted['bonds'],storage='replicated keyed history; owned-coordinate collective, no ghost history slots',wall_s=time.perf_counter()-start,RSS_KiB_per_rank=COMM.allgather(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
if mode=='mpi':
 across=oa[edges[:,0]]!=oa[edges[:,1]];result.update(actual_crossrank_bonds=int(across.sum()),crossrank_partial_before=int(np.sum(across&(d>0)&(d<1))),crossrank_failed_after=int(np.sum(across&(dn==1))),owner_before=ob,owner_after=oa)
 # Compare single-rank actual run of this identical imported state and commit.
 single=np.load(O/'identity_single_reference.npz');result['single_vs_MPI_force_max_N']=float(abs(single['force']-accepted['force']).max());result['single_vs_MPI_U_abs_J']=abs(float(single['U_J'])-accepted['scalars'][0])
else:np.savez_compressed(O/'identity_single_reference.npz',force=accepted['force'],U_J=accepted['scalars'][0],d=accepted['d'],kappa=accepted['kappa'])
passed=triggered and safety and result['original_bond_ids_unchanged'] and result['damage_sort_or_migration_error']==0 and result['reference_m_error']==0 and result['force_error_N']<2e-17 and result['commit_native_offline_force_max_N']<2e-17 and result['commit_native_offline_damage_max']<2e-13 and result['commit_native_offline_D_abs_J']<1e-20 and result['actual_new_failure_events']==result['expected_new_failure_events'] and result['repeat_did_not_increment']
if mode=='mpi':passed=passed and result['crossrank_partial_before']>0 and result['single_vs_MPI_force_max_N']<2e-17
result['status']='PASS' if passed else result['status'] if not triggered else 'FAIL'
if rank==0:
 dump(O/('ATOM_SORT_DAMAGE_IDENTITY.json' if mode=='sort' else 'MPI_DAMAGE_HISTORY_AUDIT.json'),result);np.savez_compressed(O/(mode+'_identity.npz'),x=x,bonds=edges,d_before=before['d'],d_after=after['d'],kappa_before=before['kappa'],kappa_after=after['kappa'],force_before=before['force'],force_after=after['force'],force_accepted=accepted['force'],d_accepted=accepted['d'],bond_ids=accepted['bond_ids']);print(mode,result['status'],changed,result['actual_new_failure_events'],flush=True)
ns.close();assert passed,result
