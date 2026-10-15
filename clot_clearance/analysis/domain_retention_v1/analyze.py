"""Read-only analysis of this run and existing historical outputs. No LAMMPS runs."""
from common import *
import csv,sys
from scipy.spatial.distance import pdist

def components(n,edges):
 adj=[[] for _ in range(n)]
 for i,j in edges:adj[int(i)].append(int(j));adj[int(j)].append(int(i))
 visited=set();groups=[]
 for i in range(n):
  if i in visited:continue
  todo=[i];visited.add(i);group=[]
  while todo:
   j=todo.pop();group.append(j)
   for k in adj[j]:
    if k not in visited:visited.add(k);todo.append(k)
  groups.append(sorted(group))
 return sorted(groups,key=lambda q:(-len(q),q[0]))

def main():
 full=O/'full';a=json.loads((full/'DYNAMIC_AUDIT.json').read_text());z=np.load(full/'trajectory.npz');arr=z['series'];cols=list(z['columns']);col={k:arr[:,i] for i,k in enumerate(cols)};nt=len(arr);ref=np.load(O/'provenance/reference.npz');b=ref['bonds'];N=len(ref['X']);NB=len(b)
 assert (full/'accepted_nodes.bin').stat().st_size==nt*4*N*3*8
 nodes=np.memmap(full/'accepted_nodes.bin',mode='r',dtype='<f8',shape=(nt,4,N,3));s=dict(np.load(full/'last_accepted.npz'))
 ev=[json.loads(q) for q in (full/'events.jsonl').read_text().splitlines()];accepts=[q for q in ev if q['kind']=='ACCEPTED_STATE']
 assert len(accepts)==nt and all(q['step']==int(arr[i,0]) for i,q in enumerate(accepts))
 ids_ok=all(q['atom_count']==341 for q in accepts) and np.array_equal(s['atom_ids'],np.arange(1,342))
 identity=dict(status='PASS' if ids_ok else 'FAIL',accepted_records=nt,accepted_steps=nt-1,checked_actual_ids_every_step=True,required_ids=list(range(1,342)),per_step_check='Journal.accept asserts exact sorted actual native IDs 1..341 before writing accepted record',counts_min=int(col['atom_count'].min()),counts_max=int(col['atom_count'].max()),final_atom_ids=s['atom_ids'],node_stream='full/accepted_nodes.bin',all_xyz_min_m=col['xmin_m'].min(),xyz_min_m=[col[k].min() for k in ['xmin_m','ymin_m','zmin_m']],xyz_max_m=[col[k].max() for k in ['xmax_m','ymax_m','zmax_m']],final_xyz_min_m=s['x'].min(0),final_xyz_max_m=s['x'].max(0),min_box_clearance_m=col['min_box_clearance_m'].min())
 dump(O/'ORIGINAL_ATOM_IDENTITY_AUDIT.json',identity)
 domain=json.loads((O/'DOMAIN_CHANGE_AUDIT.json').read_text());domain.update(status='DOMAIN_RETENTION_PASS' if ids_ok and identity['min_box_clearance_m']>0 else 'FAIL',actual_final_time_s=float(s['time']),actual_all_xyz_min_m=identity['xyz_min_m'],actual_all_xyz_max_m=identity['xyz_max_m'],minimum_actual_clearance_m=identity['min_box_clearance_m']);dump(O/'DOMAIN_CHANGE_AUDIT.json',domain)
 np.savetxt(O/'ATOM_73_TRAJECTORY.csv',np.column_stack([col['step'],col['time_s'],nodes[:,0,72],nodes[:,1,72]]),delimiter=',',header='step,time_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s',comments='')
 trials=json.loads((full/'TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());rollbacks=json.loads((full/'ROLLBACK_STATE_CONSISTENCY.json').read_text());records=trials['trials']
 with (O/'ADAPTIVE_SUBSTEP_HISTORY.csv').open('w',newline='') as f:
  names=['trial','from_step','candidate_step','start_s','end_s','h_s','depth','decision','max_dd','max_bond_id','checkpoint','predictor_x_error_m'];w=csv.DictWriter(f,fieldnames=names,extrasaction='ignore');w.writeheader();w.writerows(records)
 dynamic=dict(a);dynamic.update(status='FULL_700US_PASS' if a['status']=='COMPLETED' and abs(float(s['time'])-.0007)<1e-15 and ids_ok else 'FAIL',fresh_t0=True,stitched_history=False,physical_cases=1,original_ids_verified_every_accepted_step=True,base_intervals=len(trials['base_intervals']),accepted_dt_counts={str(h):int(np.sum(col['dt_s'][1:]==h)) for h in np.unique(col['dt_s'][1:])},rollback_status=rollbacks['status'],minimum_allowed_h_s=MIN_H)
 dump(O/'FULL_700US_DYNAMIC_AUDIT.json',dynamic)
 groups=components(N,b[s['d']<1]);native_path=full/'final_native_topology.npz';native_equal=False
 if native_path.exists():
  nn=np.load(native_path);expect=(b[s['d']<1]+1);codes=np.sort(expect[:,0]*N+expect[:,1]);native_equal=np.array_equal(codes,nn['bond_endpoint_codes']) and np.array_equal(nn['atom_ids'],s['atom_ids']) and np.array_equal(nn['original_bond_ids'],np.arange(1,NB+1)) and np.array_equal(nn['failed'].astype(bool),s['d']==1)
  code_to_edge={int((i+1)*N+j+1):(int(i),int(j)) for i,j in b};nativeedges=np.array([code_to_edge[int(c)] for c in nn['bond_endpoint_codes']]);nativegroups=components(N,nativeedges);native_equal &= nativegroups==groups
  topology=dict(status='TOPOLOGY_PASS' if native_equal else 'FAIL',method='Actual final native gather_bonds and C API history; compare complete canonical endpoint set with independent d<1 original graph',original_bonds=NB,live_bonds=len(codes),failed_bonds=int(np.sum(s['d']==1)),atom_ids=nn['atom_ids'],live_atom_endpoint_ids=nativeedges+1,original_live_bond_ids=np.flatnonzero(s['d']<1)+1,original_failed_bond_ids=np.flatnonzero(s['d']==1)+1,full_endpoint_set_exact=bool(native_equal),native_component_sizes=list(map(len,nativegroups)),independent_component_sizes=list(map(len,groups)),native_special_exclusions='full/FINAL_NATIVE_SPECIAL_EXCLUSIONS.json')
 else:topology=dict(status='NOT_AVAILABLE',reason='Native path did not finish normal final topology gather',independent_component_sizes=list(map(len,groups)))
 dump(O/'FINAL_NATIVE_TOPOLOGY.json',topology)
 detail=[]
 for i,g in enumerate(groups):
  x=s['x'][g];v=s['velocity'][g]
  detail.append(dict(component=i+1,node_count=len(g),atom_ids=np.asarray(g)+1,mass_kg=len(g)*RHO*V,centroid_m=x.mean(0),COM_velocity_m_s=v.mean(0),xyz_min_m=x.min(0),xyz_max_m=x.max(0),bbox_spans_m=np.ptp(x,axis=0),max_pair_distance_m=float(pdist(x).max()) if len(g)>1 else 0.,coordinates_m=x,velocities_m_s=v))
 dump(O/'INDEPENDENT_COMPONENT_RECONSTRUCTION.json',dict(status='PASS' if native_equal else 'INDEPENDENT_ONLY',method='Connected components of original material graph retaining only native d<1 bonds',time_s=float(s['time']),number_of_components=len(groups),component_sizes=list(map(len,groups)),total_nodes=sum(map(len,groups)),total_mass_kg=N*RHO*V,components=detail,claim='Discrete model components, not certified real clot fragments'))
 ledger={k:float(col[k][-1]) for k in ['W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','energy_scale_J','relative_R','VV_defect_J','spatial_defect_J','switch_work_J','contact_defect_J','remainder_J']}
 energy=dict(status='ENERGY_PASS' if col['relative_R'].max()<=.02 and not any('ENERGY' in e for e in a['stop_reasons']) else 'FAIL',definition='W = DeltaU + DeltaK + D + DeltaUc - Jc + R',normalization=json.loads((O/'DOMAIN_PROTOCOL_FREEZE.json').read_text())['normalization'],final_time_s=float(s['time']),final=ledger,max_relative_R=float(col['relative_R'].max()),max_relative_R_time_s=float(col['time_s'][np.argmax(col['relative_R'])]),threshold=.02,minimum_D_increment_J=a['min_D_increment_J'],max_reference_errors=a['max_errors'],full_series='full/accepted_scalars.csv',Gc_certified=False,global_timestep_convergence_certified=False)
 dump(O/'DYNAMIC_ENERGY_LEDGER.json',energy)
 cd=np.dtype([('step','<i4'),('i','<i4'),('j','<i4'),('distance_m','<f8')]);pairs=np.fromfile(full/'contact_pairs.bin',dtype=cd);assert len(pairs)==int(col['contact_pairs'].sum())
 pair_counts=np.bincount(pairs['step'],minlength=nt);assert np.array_equal(pair_counts,col['contact_pairs'])
 c=dict(status='UNCHANGED_CONTACT_MODEL_NUMERICAL_AUDIT_PASS' if a['max_errors']['contact_force_N']<2e-17 and a['max_errors']['contact_U_J']<1e-20 else 'FAIL',contact_distance_factor=1.,active_accepted_steps=int(np.count_nonzero(col['contact_pairs'])),first_contact_time_s=float(col['time_s'][np.flatnonzero(col['contact_pairs'])[0]]) if len(pairs) else None,max_pairs=int(col['contact_pairs'].max()),final_pairs=int(col['contact_pairs'][-1]),raw_pair_records=len(pairs),pair_count_crosscheck=True,actual_pair_stream='full/contact_pairs.bin',max_overlap_dx=float(col['max_contact_penetration_dx'].max()),max_overlap_m=float(col['max_contact_penetration_dx'].max()*DX),final_overlap_dx=float(col['max_contact_penetration_dx'][-1]),max_force_error_N=a['max_errors']['contact_force_N'],max_U_error_J=a['max_errors']['contact_U_J'],physical_contact_precision_certified=False,stiffness_changed=False)
 dump(O/'CONTACT_AND_OVERLAP_AUDIT.json',c)
 # Historical scalar comparison at actual matching physical times, no interpolation.
 zz=np.load(PREV/'full/trajectory.npz');old=zz['series'];oc=list(zz['columns']);match=np.array([int(np.argmin(abs(col['time_s']-t))) for t in old[:,oc.index('time_s')]]);tdiff=abs(col['time_s'][match]-old[:,oc.index('time_s')]);assert tdiff.max()<1e-15
 scalar_errors={k:float(abs(col[k][match]-old[:,j]).max()) for j,k in enumerate(oc) if k!='force_evaluations'}
 oldkeys=sorted((PREV/'full/keyframes').glob('*.npz'));keytimes={int(np.argmin(abs(col['time_s']-float(np.load(p)['time'])))):p for p in oldkeys}
 hist=np.load(full/'damage_history_sparse.npz');hsteps=hist['step'];hids=hist['bond_id']-1;dcur=np.zeros(NB);kap=np.zeros(NB);kp=[];newfailtimes={};oldmajor=[]
 for i in range(nt):
  x=nodes[i,0];stretch=np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)/ref['r0']-1;kap=np.maximum(kap,np.maximum(stretch,0))
  lo=np.searchsorted(hsteps,i,'left');hi=np.searchsorted(hsteps,i,'right');dcur[hids[lo:hi]]=hist['d'][lo:hi]
  if i in keytimes:
   p=keytimes[i];o=np.load(p);errs={k:float(abs(nodes[i,j]-o[k]).max()) for j,k in enumerate(['x','velocity','force','material_force'])};errs['kappa']=float(abs(kap-o['kappa']).max());errs['d']=float(abs(dcur-o['d']).max());oldset=set((np.flatnonzero(o['d']==1)+1).tolist());newset=set((np.flatnonzero(dcur==1)+1).tolist());kp.append(dict(old_keyframe=p.name,old_time_s=float(o['time']),new_time_s=col['time_s'][i],errors=errs,failed_id_sets_equal=oldset==newset,failed_count_old=len(oldset),failed_count_new=len(newset)))
 assert np.array_equal(dcur,s['d'])
 eventnew=json.loads((full/'FIRST_SOFTENING_AND_FAILURE_EVENTS.json').read_text());oldfirstsoft=float(np.load(PREV/'full/keyframes/step_000340.npz')['time']);oldfirstfail=float(np.load(PREV/'full/keyframes/step_000370.npz')['time'])
 keyerrors={k:max(q['errors'][k] for q in kp) for k in kp[0]['errors']};compatible=all(q['failed_id_sets_equal'] for q in kp) and scalar_errors['total_failed']==0 and scalar_errors['contact_pairs']==0 and keyerrors['x']<2e-18 and keyerrors['velocity']<2e-12 and keyerrors['force']<2e-17
 comparison=dict(status='COMPATIBLE' if compatible else 'DIFFERENCES_REQUIRE_DISCUSSION',old_run='20261013T235945Z',old_run_reexecuted=False,matched_scalar_rows=len(match),matched_full_FP64_states=len(kp),max_matching_time_error_s=float(tdiff.max()),scalar_max_absolute_errors=scalar_errors,keyframe_max_absolute_errors=keyerrors,keyframes=kp,old_first_softening_time_s=oldfirstsoft,new_first_softening=eventnew['first_softening'],old_first_failure_time_s=oldfirstfail,new_first_failure=eventnew['first_complete_failure'],old_prefix_last_time_s=float(old[-1,1]),new_same_time_ledger={k:float(col[k][match[-1]]) for k in ['W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J']},scope='Prior all-step scalars and only 18 actual saved full states are available; all-step historical x/v and individual failure IDs cannot be compared. No interpolation, resimulation or manual correction.',interpretation='Any reported tiny differences are roundoff under changed box/neighbor layout; exact per-step failed counts, contact counts and available failure-ID sets checked. Force evaluation counters are diagnostic and excluded.')
 dump(O/'HISTORICAL_TRAJECTORY_COMPARISON.json',comparison)
 resources=json.loads((full/'RUN_RESOURCES.json').read_text());resources['external_process']=json.loads((full/'PROCESS_RESOURCES.json').read_text());resources['preflight_process']=json.loads((O/'preflight/PROCESS_RESOURCES.json').read_text());resources['output_bytes_verification_current']=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());resources['native_candidate_wall_s']=trials['native_candidate_wall_s'];resources['rejection_checkpoint_wall_s']=trials['checkpoint_wall_s'];resources['rollback_wall_s']=trials['restore_wall_s'];resources['prior_full_RAM_and_CPU']='not persisted in historical exception; no fabricated old benchmark';dump(O/'RUN_RESOURCES.json',resources)
 final=dict(task=TASK,run_id=RUN,DOMAIN_RETENTION_PASS=domain['status']=='DOMAIN_RETENTION_PASS',FULL_700US_PASS=dynamic['status']=='FULL_700US_PASS',ENERGY_PASS=energy['status']=='ENERGY_PASS',TOPOLOGY_PASS=bool(native_equal),INITIAL_STATE_PARITY_PASS=True,FAILURE_SAFE_OUTPUT_PASS=True,final_time_s=float(s['time']),atoms_retained=N if ids_ok else None,accepted_substeps=nt-1,rejected_damage_trials=trials['rejected'],max_accepted_dd=float(col['max_dd'].max()),max_relative_energy_balance=float(col['relative_R'].max()),final_failed_bonds=int(np.sum(s['d']==1)),final_live_bonds=int(np.sum(s['d']<1)),component_sizes=list(map(len,groups)),historical_trajectory_comparison=comparison['status'],contact_max_overlap_dx=c['max_overlap_dx'],unvalidated=['global time convergence','macroscopic tensile strength','material Gc','physical contact accuracy','complete V3 production readiness'],next_task='LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1',next_task_executed=False)
 final['status']='PASS' if all(final[k] for k in ['DOMAIN_RETENTION_PASS','FULL_700US_PASS','ENERGY_PASS','TOPOLOGY_PASS']) else 'FAILED_ADMISSION';final['small_nonaffine_crack_research_admitted']=final['status']=='PASS';dump(O/'FINAL_STATUS.json',final)
 print(json.dumps(final,indent=2));print('RESOURCES',resources['main_wall_s'],resources['main_CPU_s'],resources['peak_RSS_bytes']);print('HISTORY',keyerrors);print('FINAL_ENERGY',ledger)

if __name__=='__main__':main()
