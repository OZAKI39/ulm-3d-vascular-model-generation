"""Read-only analysis of the single stopped physical run; never instantiate LAMMPS."""
from common import *
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist
import shutil,csv,difflib,gzip

ref=dict(np.load(O/'provenance/reference.npz'));s=dict(np.load(O/'full/last_accepted.npz'))
hist=np.genfromtxt(O/'NATIVE_DYNAMIC_HISTORY.csv',delimiter=',',names=True)
audit=json.loads((O/'full/DYNAMIC_AUDIT.json').read_text());pre=json.loads((O/'MINIMAL_IMPLEMENTATION_CHECK.json').read_text())
resources=json.loads((O/'RUN_RESOURCES.json').read_text());freeze=json.loads((O/'SMALL_CRACK_DESIGN_FREEZE.json').read_text())
X=ref['X'];b=ref['bonds'];cut=ref['initial_cut'];natural=~cut;live=s['d']<1;N=len(X)
assert len(hist)==201 and int(s['step'])==200 and np.array_equal(s['d'],cut.astype(float))
assert np.array_equal(s['atom_ids'],ref['atom_ids']) and np.all(hist['atom_count']==N)
assert np.all(hist['new_fail']==0) and np.all(hist['new_soft']==0) and np.all(hist['D_J']==0)
assert np.all(s['kappa'][cut]==SC) and np.all(s['kappa'][natural]<S0)
assert np.array_equal(b,np.unique(b,axis=0)) and np.allclose(ref['r0'],np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1),rtol=0,atol=0)
frameidx=np.genfromtxt(O/'full/native_frames.csv',delimiter=',',names=True)
frames=np.memmap(O/'full/native_frames.bin',dtype='<f8',mode='r',shape=(len(frameidx),4,N,3))
assert len(frameidx)==11 and np.array_equal(frames[-1,0],s['x'])
ci=ref['cross_ids'];cells=ref['cross_cells'];numcells=COUNTS[0]*COUNTS[1]
cross_count=np.bincount(cells,minlength=numcells)
alive_count=np.bincount(cells[live[ci]],minlength=numcells)
cut_count=np.bincount(cells[cut[ci]],minlength=numcells)
opening=s['x'][ref['face_upper'],2]-s['x'][ref['face_lower'],2]-DX
centers=X[ref['face_upper'],:2];new_cells=centers[:,0]>TIP
candidate=new_cells&(alive_count==0)&(opening>=.1*DX)
assert not candidate.any()
np.savez_compressed(O/'PLANE_CELL_AUDIT.npz',cell_centers_xy_m=centers,original_count=cross_count,initial_cut_count=cut_count,live_count=alive_count,relative_opening_m=opening,new_area_eligible_cell=new_cells,qualifying_cell=candidate)
np.savez_compressed(O/'ORIGINAL_TOPOLOGY_SETS.npz',initial_cut_bond_ids=ref['bond_ids'][cut],natural_failed_bond_ids=ref['bond_ids'][natural&~live],live_cross_bond_ids=ref['bond_ids'][ref['cross']&live])
bb=b[live];g=coo_matrix((np.ones(len(bb)),(bb[:,0],bb[:,1])),shape=(N,N)).tocsr();nc,labels=connected_components(g,directed=False)
components=[]
for k in range(nc):
 ids=np.flatnonzero(labels==k);q=s['x'][ids]
 components.append(dict(component=k,node_count=len(ids),mass_kg=len(ids)*RHO*V,atom_ids=ids+1,center_of_mass_m=q.mean(0),actual_xyz_min_m=q.min(0),actual_xyz_max_m=q.max(0),actual_center_span_m=np.ptp(q,axis=0),max_pair_distance_m=float(pdist(q).max()) if len(q)>1 else 0.))
np.savez_compressed(O/'COMPONENT_LABELS.npz',atom_ids=ref['atom_ids'],labels=labels)
dump(O/'NATIVE_COMPONENT_TOPOLOGY.json',dict(status='INDEPENDENT_REFERENCE_RECONSTRUCTION_LAST_ACCEPTED',time_s=float(s['time']),components=components,initial_cut_bonds=int(cut.sum()),natural_failed_bonds=0,live_bonds=int(live.sum()),native_bond_count_matches_every_accepted_step=bool(np.all(hist['live_bonds']==live.sum())),native_initial_and_P1_restart_endpoint_match=pre['restart_native_topology_and_special_exact'],native_final_endpoint_gather='NOT_PERFORMED_AFTER_RESOURCE_STOP',raw_native_final_topology_matches_false_means='Not evaluated; no observed endpoint mismatch',all_original_atom_ids_preserved=True,partly_softened_bonds_remain_connected=True))
natids=np.flatnonzero(natural);order=natids[np.argsort(s['kappa'][natids])[-12:][::-1]]
stretch=np.linalg.norm(s['x'][b[:,1]]-s['x'][b[:,0]],axis=1)/ref['r0']-1
maxrec=[dict(bond_id=int(j+1),atom_ids=b[j]+1,reference_midpoint_m=X[b[j]].mean(0),actual_midpoint_m=s['x'][b[j]].mean(0),region=REGIONS[int(ref['bond_regions'][j])],kappa=float(s['kappa'][j]),instantaneous_stretch=float(stretch[j]),d=float(s['d'][j]),label='UNFAILED_OBSERVATION_NOT_FIRST_EVENT') for j in order]
dump(O/'FIRST_SOFTENING_AND_FAILURE.json',dict(status='NOT_OBSERVED_IN_COMPLETED_WINDOW',window_s=[0,float(s['time'])],first_softening=None,first_complete_failure=None,first_localized_damage_cloud=None,first_continuous_crack_advance=None,maximum_natural_kappa=float(s['kappa'][natural].max()),maximum_natural_instantaneous_stretch=float(stretch[natural].max()),initial_cut_kappa_SC_is_not_a_natural_event=True,highest_unfailed_kappa_bonds_at_last_state=maxrec,planned_protocol_completed=False))
regional=[]
dn=np.load(O/'full/accepted_node_D.npz')['node_D_J']
for k,name in enumerate(REGIONS):
 m=(ref['bond_regions']==k)&natural
 regional.append(dict(region=name,natural_reference_bonds=int(m.sum()),naturally_softened=0,naturally_failed=0,max_natural_kappa=float(s['kappa'][m].max()) if m.any() else None,star_center_D_J=float(dn[ref['node_regions']==k].sum())))
dump(O/'CRACK_LOCALIZATION_AUDIT.json',dict(classification='NO_NATURAL_DAMAGE_IN_EARLY_WINDOW',first_region_order='UNDETERMINED',notch_before_grip=None,natural_crack_propagation=False,damage_band_outside_precut=False,regions=regional,classification_basis=freeze['scientific_criteria']['region_definition'],no_full_protocol_inference=True))
dump(O/'CROSS_PLANE_BRIDGE_AUDIT.json',dict(status='BRIDGING_REMAINS',plane='original z=0',initial_reference_cross_plane_bonds=len(ci),initial_precut_cross_bonds=int(cut[ci].sum()),initial_live_cross_bonds=int((~cut[ci]).sum()),last_live_cross_bonds=int(live[ci].sum()),natural_failed_cross_bonds=0,all_original_lengths_included=True,long_live_cross_bonds_r0_greater_than_dx=int(np.sum(ref['cross']&live&(ref['r0']>DX*(1+1e-12)))),last_relative_opening_min_m=float(opening.min()),last_relative_opening_max_m=float(opening.max()),opening_threshold_m=.1*DX,candidate_new_cells=int(candidate.sum()),cell_data='PLANE_CELL_AUDIT.npz',id_sets='ORIGINAL_TOPOLOGY_SETS.npz',front_advance_observed=False))
dump(O/'CRACK_AREA_IDENTIFIABILITY.json',dict(status='CRACK_AREA_UNRESOLVED',PROJECTED_SEPARATION_AREA_m2=None,initial_precut_projected_area_m2=(TIP+SIZE[0]/2)*SIZE[1],initial_area_excluded_from_new_area=True,qualifying_new_separation_cells=0,qualifying_grid_proxy_area_m2=0.,proxy_zero_does_not_establish_a_measured_unique_physical_area=True,reason='No natural damage through 50 us; all 4576 initially live cross-plane bonds remain alive. Observation ended at resource admission stop.',finite_horizon_m=HORIZON,cell_size_m=DX,criterion_frozen=True,unique_surface_from_damage_cloud=False,G_SEPARATION_PROTOCOL='UNAVAILABLE',GC_EFFECTIVE='UNAVAILABLE'))
last={k:float(hist[k][-1]) for k in hist.dtype.names}
uvalid=hist['U_J']>=.01*hist['U_J'].max()
energy=dict(status='PASS_FOR_ACCEPTED_0_TO_50_US_ONLY',equation='W = DeltaU + DeltaK + D + DeltaUc - Jc + R',initial_U_K_Uc_J=[0.,0.,0.],last={k:last[k] for k in ['W_J','U_J','K_J','Kfree_J','D_J','Uc_J','Jc_J','R_J']},max_absolute_R_J=float(abs(hist['R_J']).max()),max_relative_R=float(hist['relative_R'].max()),max_relative_R_percent=float(hist['relative_R'].max()*100),relative_R_gate=.02,last_K_over_U=last['K_J']/last['U_J'],max_K_over_U_where_U_ge_1pct_peak=float(np.max(hist['K_J'][uvalid]/hist['U_J'][uvalid])),inertia_classification='DYNAMIC_EARLY_TRANSIENT_NOT_QUASISTATIC',W_definition='Prescribed z DOFs only: -0.5*(F_old+F_new)*du plus change of kinetic energy of prescribed z components. Includes actuator acceleration.',regional_D_method='Exact declared endpoint-star center partition, D_i = U_i(x_new,d_old)-U_i(x_new,d_new); state coupling makes attribution convention-dependent, not unique bond/surface work.',regional_D=regional,initial_cut_D_J=0.,G_SEPARATION_PROTOCOL='UNAVAILABLE',GC_EFFECTIVE='UNAVAILABLE')
dump(O/'DYNAMIC_ENERGY_LEDGER.json',energy)
dump(O/'NUMERICAL_VALIDITY.json',dict(status='PARTIAL_MECHANICS_PASS_RESOURCE_ADMISSION_FAIL',planned_end_s=END,actual_end_s=float(s['time']),protocol_completed=False,accepted_steps=200,rejections=0,actual_dt_range_s=[float(hist['dt_s'][1:].min()),float(hist['dt_s'][1:].max())],max_dd=float(hist['max_dd'].max()),damage_event_subdivision='NOT_EXERCISED_NO_DAMAGE',max_relative_R=energy['max_relative_R'],max_absolute_R_J=energy['max_absolute_R_J'],max_momentum_residual_Ns=float(hist['momentum_residual_Ns'].max()),min_box_clearance_m=float(hist['min_box_clearance_m'].min()),max_speed_m_s=float(hist['max_speed_m_s'].max()),max_displacement_m=float(hist['max_displacement_m'].max()),last_nonaffine_RMS_m=last['nonaffine_RMS_m'],all_3080_atom_ids_retained=True,initial_tangent_limit_s=pre['initial_tangent_Verlet_limit_s'],max_independent_errors=audit['max_accepted_errors'],max_contact_pairs=int(hist['contact_pairs'].max()),max_contact_overlap_dx=float(hist['max_contact_penetration_dx'].max()),irreversibility_pass=True,original_identity_pass=True,native_final_endpoint_check='UNAVAILABLE',preflight_resource_admission_missed=True,full_fracture_or_timestep_convergence_certified=False))
# Resource diagnosis uses existing source and logs only. No native calls or rerun.
src=NATIVE_PREV/'source/src/special.cpp';shutil.copy2(src,O/'source_identity/resource_special.cpp')
lines=src.read_text().splitlines();evidence=[]
for lo,hi in [(115,140),(505,545)]:evidence.append('\n'.join(f'{i+1}: {lines[i]}' for i in range(lo-1,min(hi,len(lines)))))
(O/'RESOURCE_SOURCE_EXCERPTS.txt').write_text('\n\n'.join(evidence)+'\n')
log=(O/'full/native.log').read_text();logstart=log[:log.index('lattice sc')]
(O/'RESOURCE_LOG_EXCERPT.txt').write_text(logstart)
diag=dict(status='RESOURCE_LIMIT_STOP_CONFIRMED',main_peak_RSS_bytes=resources['peak_RSS_bytes'],main_peak_RSS_GiB=resources['peak_RSS_bytes']/1024**3,budget_GiB=4,preflight_peak_RSS_bytes=pre['peak_RSS_bytes'],resource_guard_stage='AFTER_STEP_200_ACCEPTED; first periodic check',preflight_admission_gap='P1 saved >4 GiB peak RSS, but freeze/main checked only numerical PASS. Resource failure should have blocked formal launch.',main_observed='read_data used default 0/0/0 special weights before formal special_bonds 0/1/1; max raw 1-2/1-3/1-4 counts 122/14718/1730928; special builder CPU 10.951 s. Later correct weights skip 1-3/1-4 (0.010 s).',source_evidence='Special::build returns early when 1-3 and 1-4 weights are 1 for both LJ and Coulomb. onefour_build creates nlocal x maxall tagint BEFORE dedup.',raw_onefour_nominal_allocation_elements=3080*1730928,raw_onefour_nominal_bytes_assuming_32bit_tagint=3080*1730928*4,allocation_scope='Nominal virtual allocation; not equal to resident RSS. Only visited portions are touched. Process high-water RSS is measured, attribution to an exact allocation was not separately instrumented.',primary_memory_cause='Strong source/log evidence for unnecessary initial three-hop special-neighbor construction; no solver force-state mismatch observed.',proposed_next_task='Move the SAME final special_bonds command before read_data in a new isolated revision; prove identical 1-2 exclusions, precut state, initial forces and energies with a short check; gate actual peak RSS before any newly authorized full run of the SAME geometry/loading.',change_or_rerun_performed_after_stop=False,immutable_source_sha256=sha(src))
dump(O/'RESOURCE_STOP_DIAGNOSIS.json',diag)
dump(O/'FAILURE_RECORD_INTERPRETATION.json',dict(raw_record='full/FAILURE_RECORD.json',raw_preserved=True,actual_stop_stage='Periodic resource guard AFTER accepted step 200',last_accepted_step=200,last_accepted_time_s=float(s['time']),unaccepted_candidate_at_stop=False,raw_candidate_not_admitted_true_is_generic_handler_field_and_not_applicable_here=True,raw_attempt_is_last_accepted_attempt_not_a_rejected_trial=True,native_final_topology_matches_false_means='Not performed after stop, not a discovered mismatch',latest_native_restart_time_s=0.,last_accepted_FP64_time_s=float(s['time']),FP64_is_not_a_native_restart=True,exit_zero_means='Exception archived and process closed; NOT scientific PASS'))
# Preserve raw event streams, normalize the legacy literal backslash-n separator for consumers.
streams={}
for name in ['events','trials']:
 p=O/'full'/f'{name}.jsonl';raw=p.read_text();decoder=json.JSONDecoder();pos=0;records=[]
 while pos<len(raw):
  while pos<len(raw) and raw[pos].isspace():pos+=1
  if raw[pos:pos+2]=='\\n':pos+=2;continue
  if pos==len(raw):break
  obj,end=decoder.raw_decode(raw,pos);records.append(obj);pos=end
 output=O/f'{name.upper()}_NORMALIZED.jsonl';output.write_text(''.join(json.dumps(z,ensure_ascii=False)+'\n' for z in records))
 streams[name]=dict(records=len(records),raw_sha256=sha(p),normalized_sha256=sha(output),literal_backslash_n_separator_count=raw.count('\\n'),raw_preserved=True)
dump(O/'EVENT_STREAM_FORMAT_AUDIT.json',streams)
status=dict(task=TASK,run_id=RUN,status='STOPPED_RESOURCE_BUDGET_SCIENTIFIC_VALIDATION_NOT_COMPLETED',delivery_status='ANALYSIS_IN_PROGRESS',formal_cases=1,accepted_steps=200,actual_end_s=float(s['time']),planned_end_s=END,completed_time_fraction=float(s['time'])/END,geometry_nodes=N,original_reference_bonds=len(b),initial_cut_bonds=int(cut.sum()),natural_softened_bonds=0,natural_failed_bonds=0,live_cross_plane_bonds=int(live[ci].sum()),component_sizes=[c['node_count'] for c in components],natural_crack_propagation_observed=False,CRACK_AREA='UNRESOLVED',G_SEPARATION_PROTOCOL='UNAVAILABLE',GC_EFFECTIVE='UNAVAILABLE',partial_mechanics='PASS_0_TO_50_US_ONLY',resource_admission='FAIL',preflight_resource_admission_gap=True,prior_341_700us_baseline='PRESERVED_PASS',macro_strength='NOT_VALIDATED',production_V3='NOT_ADMITTED',rerun_after_stop=False,next_task='Same-physics special-neighbor initialization memory fix plus early resource admission; no automatic rerun')
dump(O/'FINAL_STATUS.json',status)
resources.update(preflight=dict(wall_s=pre['wall_s'],CPU_s=pre['CPU_s'],peak_RSS_bytes=pre['peak_RSS_bytes'],physical_steps=3),parent_process=json.loads((O/'full/PROCESS_RESOURCES.json').read_text()),resource_admission='FAIL',resource_diagnosis='RESOURCE_STOP_DIAGNOSIS.json',postprocessing_contains_no_native_steps=True)
dump(O/'RUN_RESOURCES.json',resources)
print(json.dumps(dict(status=status['status'],max_relative_R_percent=energy['max_relative_R_percent'],max_natural_kappa=float(s['kappa'][natural].max()),opening_max_m=float(opening.max()),components=nc,event_streams=streams),indent=2))
