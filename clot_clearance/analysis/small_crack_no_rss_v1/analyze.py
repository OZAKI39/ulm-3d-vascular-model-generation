"""Analyze the explicitly authorized continuation; no LAMMPS calls."""
from common import *
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
import gzip,csv

def analyze():
 ref=dict(np.load(O/'provenance/reference.npz'));s=dict(np.load(O/'full/last_accepted.npz'));h=np.genfromtxt(O/'NATIVE_DYNAMIC_HISTORY.csv',delimiter=',',names=True)
 audit=json.loads((O/'full/DYNAMIC_AUDIT.json').read_text());freeze=json.loads((O/'SMALL_CRACK_DESIGN_FREEZE.json').read_text());events=json.loads((O/'full/FIRST_EVENTS.json').read_text());resources=json.loads((O/'RUN_RESOURCES.json').read_text())
 resource_stop=any('WALL_OR_OUTPUT_BUDGET' in z for z in audit['stop_reasons']);stop_category=('WALL_CLOCK_BUDGET' if resources['main_wall_s']>freeze['budget']['main_wall_s'] else 'OUTPUT_BUDGET') if resource_stop else ('NONE' if not audit['stop_reasons'] else 'NUMERICAL_OR_RUNTIME_GATE')
 X=ref['X'];b=ref['bonds'];cut=ref['initial_cut'];nat=~cut;N=len(X);M=len(b);complete=audit['status']=='COMPLETED' and abs(float(s['time'])-END)<1e-15
 assert len(h)==int(s['step'])+1 and np.all(h['atom_count']==N) and np.array_equal(s['atom_ids'],ref['atom_ids'])
 assert np.all(s['d'][cut]==1) and np.all(s['kappa'][cut]>=SC)
 assert np.sum(s['d']<1)==s['native_bonds']==h['live_bonds'][-1]
 assert np.all(h['live_bonds']==M-cut.sum()-h['total_natural_failed'])
 assert np.all(np.diff(h['D_J'])>=-1e-20)
 # Exact chronological d updates reconstruct all accepted topology changes.
 dtype=np.dtype([('step','<i4'),('bond_id','<i4'),('kappa','<f8'),('d','<f8')]);raw=gzip.open(O/'full/damage_updates.bin.gz','rb').read();updates=np.frombuffer(raw,dtype=dtype)
 if len(updates):assert np.all(np.diff(updates['step'])>=0) and not np.any(cut[updates['bond_id']-1])
 starts=np.r_[0,np.flatnonzero(np.diff(updates['step']))+1,len(updates)] if len(updates) else np.array([0]);group=0;d=cut.astype(float)
 fi=np.genfromtxt(O/'full/native_frames.csv',delimiter=',',names=True);frames=np.memmap(O/'full/native_frames.bin',dtype='<f8',mode='r',shape=(len(fi),4,N,3));assert np.array_equal(frames[-1,0],s['x'])
 ci=ref['cross_ids'];cells=ref['cross_cells'];centers=X[ref['face_upper'],:2];newcell=(centers[:,0]>TIP).reshape(COUNTS[:2]);tipcol=int(np.flatnonzero(newcell[:,0])[0]);interior=abs(centers.reshape(*COUNTS[:2],2)[0,:,1])<SIZE[1]/2-HORIZON
 priorarea=[];plane_arrays=[];exact_component_sizes=[];firstprop=None;previous_front=None;previous_time=None
 for f,row in enumerate(fi):
  step=int(row['step'])
  while group+1<len(starts) and int(updates['step'][starts[group]])<=step:
   z=updates[starts[group]:starts[group+1]];ids=z['bond_id']-1;assert np.all(z['d']>=d[ids]);d[ids]=z['d'];group+=1
  assert np.sum((d==1)&nat)==h['total_natural_failed'][step]
  alive=d[ci]<1;counts=np.bincount(cells[alive],minlength=len(centers)).reshape(COUNTS[:2]);x=frames[f,0];gap=(x[ref['face_upper'],2]-x[ref['face_lower'],2]-DX).reshape(COUNTS[:2]);eligible=newcell&(counts==0)&(gap>=.1*DX);connected=np.zeros_like(eligible)
  # Frozen x-contiguous prefix from original front in each thickness strip.
  for y in range(COUNTS[1]):
   for ix in range(tipcol,COUNTS[0]):
    if not eligible[ix,y]:break
    connected[ix,y]=True
  if np.any(alive):distance=cKDTree(ref['cross_intersections'][alive,:2]).query(centers)[0].reshape(COUNTS[:2])
  else:distance=np.full(COUNTS[:2],np.inf)
  strict=connected&(distance>HORIZON+DX/np.sqrt(2))
  front=TIP+connected.sum(axis=0)*DX
  if previous_front is not None and firstprop is None:
   persistent=(front-TIP>=2*DX-1e-15)&(previous_front-TIP>=2*DX-1e-15)&interior
   if persistent.any():
    firstprop=dict(step=step,time_s=float(row['time_s']),previous_saved_time_s=previous_time,thickness_strip_ids=np.flatnonzero(persistent),candidate_front_x_m=front[persistent],scope='Two saved times, interior strip >=2dx, original-plane no-bridge/opening/notch-prefix rule; not unique 3D crack surface')
    save_state(O/'front_candidate_visual_state.npz',dict(x=np.asarray(frames[f,0]),velocity=np.asarray(frames[f,1]),force=np.asarray(frames[f,2]),material_force=np.asarray(frames[f,3]),d=d.copy(),step=step,time=float(row['time_s']),scope='Actual saved nodal frame plus exact chronological d reconstruction; no complete kappa history, not a native restart'))
  previous_front=front;previous_time=float(row['time_s'])
  bb=b[d<1];g=coo_matrix((np.ones(len(bb)),(bb[:,0],bb[:,1])),shape=(N,N)).tocsr();nc,labels=connected_components(g,directed=False);assert nc==int(h['components'][step]);sizes=np.sort(np.bincount(labels))[::-1];exact_component_sizes.append(sizes.tolist())
  record=dict(step=step,time_s=float(row['time_s']),live_cross=int(alive.sum()),naturally_failed_cross=int(np.sum((d[ci]==1)&~cut[ci])),candidate_open_cells=int(eligible.sum()),notch_connected_cells=int(connected.sum()),conservative_cells=int(strict.sum()),notch_connected_proxy_area_m2=float(connected.sum()*DX**2),conservative_proxy_area_m2=float(strict.sum()*DX**2),max_front_advance_m=float((front-TIP).max()),min_gap_m=float(gap.min()),max_gap_m=float(gap.max()),components=int(nc))
  priorarea.append(record);plane_arrays.append((counts.copy(),gap.copy(),connected.copy(),strict.copy(),front.copy()))
 assert group+1==len(starts) and np.array_equal(d,s['d'])
 with (O/'CRACK_PLANE_HISTORY.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(priorarea[0]));w.writeheader();w.writerows(priorarea)
 np.savez_compressed(O/'PLANE_HISTORY.npz',time_s=fi['time_s'],steps=fi['step'],cell_centers_xy_m=centers,live_counts=np.array([p[0] for p in plane_arrays]),opening_m=np.array([p[1] for p in plane_arrays]),notch_connected=np.array([p[2] for p in plane_arrays]),conservative=np.array([p[3] for p in plane_arrays]),front_x_m=np.array([p[4] for p in plane_arrays]),initial_cross_counts=np.bincount(cells,minlength=len(centers)).reshape(COUNTS[:2]),initial_cut_counts=np.bincount(cells[cut[ci]],minlength=len(centers)).reshape(COUNTS[:2]))
 # Last accepted full reference graph and actual component geometry.
 components=[]
 for k in range(nc):
  ids=np.flatnonzero(labels==k);q=s['x'][ids];components.append(dict(label=k,node_count=len(ids),atom_ids=ids+1,mass_kg=len(ids)*RHO*V,center_of_mass_m=q.mean(0),xyz_min_m=q.min(0),xyz_max_m=q.max(0),center_span_m=np.ptp(q,axis=0),max_pair_distance_m=float(pdist(q).max()) if len(ids)>1 else 0.))
 components.sort(key=lambda q:q['node_count'],reverse=True);np.savez_compressed(O/'COMPONENT_LABELS.npz',labels=labels,atom_ids=ref['atom_ids'])
 native_final=None
 if (O/'full/final_native_topology.npz').exists():
  z=np.load(O/'full/final_native_topology.npz');bb=b[s['d']<1]+1;native_final=bool(np.array_equal(z['endpoint_codes'],np.sort(bb[:,0]*N+bb[:,1])) and np.array_equal(z['original_bond_ids'],ref['bond_ids']) and np.array_equal(z['reference_m'],ref['m']));assert native_final
 dump(O/'NATIVE_COMPONENT_TOPOLOGY.json',dict(status='PASS_NATIVE_AND_INDEPENDENT' if native_final else 'INDEPENDENT_LAST_ACCEPTED_NATIVE_FINAL_NOT_AVAILABLE',native_final_endpoint_match=native_final,components=components,component_sizes=[q['node_count'] for q in components],atom_count=N,original_reference_bonds=M,initial_cut=int(cut.sum()),natural_failed=int(np.sum((s['d']==1)&nat)),partly_softened_still_connected=True,live_bonds=int(np.sum(s['d']<1)),sampled_frame_component_sizes=exact_component_sizes))
 np.savez_compressed(O/'ORIGINAL_TOPOLOGY_SETS.npz',initial_cut_ids=ref['bond_ids'][cut],natural_softened_ids=ref['bond_ids'][nat&(s['d']>0)],natural_failed_ids=ref['bond_ids'][nat&(s['d']==1)],live_cross_ids=ref['bond_ids'][ref['cross']&(s['d']<1)])
 firstsoft=events['global_events'].get('softening');firstfail=events['global_events'].get('failure');conc=events['first_observed_concentration'];regional=[];dn=np.load(O/'full/accepted_node_D.npz')['node_D_J']
 for k,name in enumerate(REGIONS):
  m=(ref['bond_regions']==k)&nat;regional.append(dict(region=name,naturally_softened=int(np.sum(m&(s['d']>0))),naturally_failed=int(np.sum(m&(s['d']==1))),maximum_kappa=float(s['kappa'][m].max()) if m.any() else None,star_center_D_J=float(dn[ref['node_regions']==k].sum()),first_softening=events['first_region_softening'].get(name),first_failure=events['first_region_failure'].get(name)))
 if firstfail and 'GRIP' in firstfail['regions']:classification='GRIP_DOMINATED_FAILURE' if set(firstfail['regions'])=={'GRIP'} else 'MIXED_FIRST_FAILURE_INCLUDING_GRIP'
 elif conc and firstsoft and any(z.startswith('NOTCH_TIP') for z in firstsoft['regions']):classification='NOTCH_LOCALIZED_DAMAGE'
 elif np.any(nat&(s['d']>0)):classification='DIFFUSE_OR_MIXED_DAMAGE'
 else:classification='NO_NATURAL_DAMAGE'
 dump(O/'FIRST_SOFTENING_AND_FAILURE.json',dict(first_softening=firstsoft,first_complete_failure=firstfail,first_localized_cloud=conc,first_candidate_front_advance=firstprop,initial_cut_excluded=True,region_first_softening=events['first_region_softening'],region_first_failure=events['first_region_failure'],time_window_s=[0,float(s['time'])]))
 gripsoft=events['first_region_softening'].get('GRIP');tipsoft=[z['time_s'] for name,z in events['first_region_softening'].items() if name.startswith('NOTCH_TIP')]
 dump(O/'CRACK_LOCALIZATION_AUDIT.json',dict(classification=classification,classification_scope='First natural onset and observed localized cloud; later regional competition is reported separately',regions=regional,first_softening_regions=None if not firstsoft else firstsoft['regions'],first_failure_regions=None if not firstfail else firstfail['regions'],notch_softening_before_grip=bool(tipsoft and (gripsoft is None or min(tipsoft)<gripsoft['time_s'])),notch_softening_lead_over_grip_s=None if not tipsoft or gripsoft is None else gripsoft['time_s']-min(tipsoft),localized_cloud_observed=conc is not None,candidate_plane_propagation_observed=firstprop is not None,first_candidate_propagation=firstprop,region_rule=freeze['scientific_criteria']['region_definition'],off_plane_natural_failed=int(np.sum(nat&(s['d']==1)&~ref['cross']))))
 lastplane=priorarea[-1];live=s['d']<1
 dump(O/'CROSS_PLANE_BRIDGE_AUDIT.json',dict(initial_cross_bonds=len(ci),initial_cut_cross=int(cut[ci].sum()),initial_live_cross=int((~cut[ci]).sum()),last_live_cross=int(live[ci].sum()),natural_failed_cross=int(np.sum(nat[ci]&~live[ci])),last_long_live_cross=int(np.sum(ref['cross']&live&(ref['r0']>DX*(1+1e-12)))),all_reference_lengths_included=True,last_plane=lastplane,front_history='PLANE_HISTORY.npz'))
 # No unique surface or Gc follows from the grid proxy alone.
 areaok=firstprop is not None and lastplane['conservative_cells']>0 and classification=='NOTCH_LOCALIZED_DAMAGE'
 dump(O/'CRACK_AREA_IDENTIFIABILITY.json',dict(status='PROJECTED_CANDIDATE_REQUIRES_GEOMETRIC_REVIEW' if areaok else 'CRACK_AREA_UNRESOLVED',PROJECTED_SEPARATION_AREA_m2=None,notch_connected_grid_proxy_m2=lastplane['notch_connected_proxy_area_m2'],conservative_grid_proxy_m2=lastplane['conservative_proxy_area_m2'],initial_cut_area_m2=(TIP+SIZE[0]/2)*SIZE[1],initial_area_excluded=True,propagation_candidate=firstprop,strict_rule='No alive crossing within delta + half-cell diagonal of cell center, opening>=0.1dx, contiguous from notch',finite_horizon_m=HORIZON,grid_m=DX,grid_ties='Original frozen FP64 floor assignment retained; exact boundary ties can alter adjacent bin counts. No new positive area certification without boundary-tie sensitivity and 3D surface review.',unique_continuous_surface_certified=False,G_SEPARATION_PROTOCOL='UNAVAILABLE',GC_EFFECTIVE='UNAVAILABLE'))
 last={k:float(h[k][-1]) for k in h.dtype.names};valid=h['U_J']>=max(.01*h['U_J'].max(),1e-30);ku=h['K_J'][valid]/h['U_J'][valid]
 energy=dict(status='PASS_ACCEPTED_STATES_ONLY',equation='W = DeltaU + DeltaK + D + DeltaUc - Jc + R',last={k:last[k] for k in ['W_J','U_J','K_J','Kfree_J','D_J','Uc_J','Jc_J','R_J']},max_absolute_R_J=float(abs(h['R_J']).max()),max_relative_R=float(h['relative_R'].max()),max_relative_R_percent=float(h['relative_R'].max()*100),gate=.02,last_K_over_U=last['K_J']/last['U_J'] if last['U_J']>0 else None,max_K_over_U_where_U_ge_1pct_peak=float(ku.max()) if len(ku) else None,regional_D=regional,regional_D_sum_J=float(dn.sum()),native_D_J=last['D_J'],regional_sum_error_J=float(dn.sum()-last['D_J']),partition='Declared endpoint-star center partition, not unique bond/surface energy',G_SEPARATION_PROTOCOL='UNAVAILABLE',GC_EFFECTIVE='UNAVAILABLE',quasistatic_certified=False)
 assert abs(energy['regional_sum_error_J'])<=freeze['gates']['energy_abs_J'];dump(O/'DYNAMIC_ENERGY_LEDGER.json',energy)
 cp=json.loads((O/'full/LATEST_NATIVE_CHECKPOINT.json').read_text());trial=json.loads((O/'full/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text())
 numerical=dict(status='FULL_PROTOCOL_NUMERICAL_PASS' if complete else 'STOPPED_AT_EXISTING_ADMISSION_GATE',complete=complete,stop_reasons=audit['stop_reasons'],accepted_steps=int(s['step']),rejected_steps=trial['rejected'],actual_time_s=float(s['time']),planned_end_s=END,max_dd=float(h['max_dd'].max()),dt_min_s=float(h['dt_s'][1:].min()),dt_max_s=float(h['dt_s'][1:].max()),max_relative_R=energy['max_relative_R'],max_momentum_residual_Ns=float(h['momentum_residual_Ns'].max()),max_speed_m_s=float(h['max_speed_m_s'].max()),max_displacement_m=float(h['max_displacement_m'].max()),min_domain_clearance_m=float(h['min_box_clearance_m'].min()),all_original_ids_preserved=True,max_independent_errors=audit['max_accepted_errors'],latest_native_checkpoint_time_s=cp['time_s'],native_final_endpoint_match=native_final,global_timestep_convergence_certified=False,memory_limit_GiB=None)
 dump(O/'NUMERICAL_VALIDITY.json',numerical)
 numerical.update(stop_category=stop_category,resource_stop_is_not_a_failed_accepted_numerical_check=resource_stop);dump(O/'NUMERICAL_VALIDITY.json',numerical)
 cdtype=np.dtype([('step','<i4'),('i','<i4'),('j','<i4'),('distance_m','<f8')]);contact=np.frombuffer(gzip.open(O/'full/contact_pairs.bin.gz','rb').read(),dtype=cdtype);active=contact['distance_m']<DX;active_steps,activecounts=np.unique(contact['step'][active],return_counts=True)
 dump(O/'CONTACT_STREAM_AUDIT.json',dict(geometric_candidate_records=len(contact),force_bearing_records=int(active.sum()),zero_force_touching_records=int(np.sum(contact['distance_m']==DX)),first_force_bearing_step=int(active_steps[0]) if len(active_steps) else None,max_force_bearing_pairs=int(activecounts.max()) if len(activecounts) else 0,max_geometric_candidates=int(h['contact_pairs'].max()),max_overlap_dx=float(h['max_contact_penetration_dx'].max()),raw_first_contact_is_geometric_candidate=True,contact_model_changed=False))
 if (O/'full/FAILURE_RECORD.json').exists():
  raw=json.loads((O/'full/FAILURE_RECORD.json').read_text());attempt=raw['attempt'];unaccepted=attempt['candidate_step']>int(s['step']);dump(O/'FAILURE_RECORD_INTERPRETATION.json',dict(raw_preserved=True,last_accepted_step=int(s['step']),last_accepted_time_s=float(s['time']),attempt=attempt,unaccepted_candidate_at_stop=unaccepted,raw_candidate_not_admitted_field_applicable=unaccepted,latest_native_checkpoint_time_s=cp['time_s'],FP64_is_not_a_native_restart=True,exit_zero_is_not_scientific_PASS=True))
 status=dict(task=TASK,run_id=RUN,status='FULL_PROTOCOL_NUMERICAL_PASS' if complete else 'STOPPED_EXISTING_ADMISSION_GATE',scientific_classification=classification,full_protocol_completed=complete,accepted_steps=int(s['step']),rejected_steps=trial['rejected'],actual_end_s=float(s['time']),planned_end_s=END,memory_cap_GiB=None,natural_softened=int(np.sum(nat&(s['d']>0))),natural_failed=int(np.sum(nat&~live)),initial_cut=int(cut.sum()),live_cross_bonds=int(live[ci].sum()),component_sizes=[q['node_count'] for q in components],max_relative_R_percent=energy['max_relative_R_percent'],first_softening_time_s=None if not firstsoft else firstsoft['time_s'],first_failure_time_s=None if not firstfail else firstfail['time_s'],crack_propagation_candidate=firstprop is not None,unique_surface_certified=False,CRACK_AREA='UNRESOLVED',G_SEPARATION_PROTOCOL='UNAVAILABLE',GC_EFFECTIVE='UNAVAILABLE',new_formal_cases_this_authorization=1,source_run=str(PRIOR_SMALL),prefix_50us='EXACT_MATCH',production_and_previous_baselines_unchanged=True,delivery_status='POSTPROCESSING')
 status['stop_category']=stop_category
 if resource_stop:status['status']='STOPPED_'+stop_category
 dump(O/'FINAL_STATUS.json',status)
 resources.update(parent_process=json.loads((O/'full/PROCESS_RESOURCES.json').read_text()),RSS_ceiling_disabled_by_user=True,no_new_CPP_build=True,no_new_preflight_physical_run=True);dump(O/'RUN_RESOURCES.json',resources)
 print(json.dumps(status,indent=2))

if __name__=='__main__':analyze()
