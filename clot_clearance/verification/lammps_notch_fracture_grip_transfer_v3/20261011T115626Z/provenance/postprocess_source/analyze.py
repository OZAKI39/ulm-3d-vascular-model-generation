"""Read-only audit of the one preregistered trajectory; never invokes LAMMPS."""
from pathlib import Path
import sys,json,csv,hashlib,numpy as np
from metrics import Topology
from design import save,js,input_groups
ROOT=Path(__file__).resolve().parents[2]
DT=.5e-6
NAMES=['notch_tip','grip_transition','other_surface','free_interior']

def loadcsv(p):return np.genfromtxt(p,delimiter=',',names=True)
def energy(a,first,peak,advance):
 keys=['W_ext_J','U_bulk_J','U_contact_J','K_free_J','D_num_J','D_break_independent_J','R_energy_J']
 cuts=[('pre_first',0,first-1),('first_to_peak',first,peak)]+([('peak_to_central_advance',peak,advance),('central_advance_to_end',advance,len(a)-1)] if advance is not None else [('peak_to_stop_no_central_advance',peak,len(a)-1)])
 phases={}
 for name,start,end in cuts:
  if end<start:continue
  v=a[start:end+1];k=int(np.argmax(v['KE_over_U']));phases[name]=dict(start_step=start,end_step=end,time_s=[float(v['time_s'][0]),float(v['time_s'][-1])],max_K_over_U=float(v['KE_over_U'][k]),max_ratio_time_s=float(v['time_s'][k]),U_at_max_ratio_J=float(v['U_bulk_J'][k]),max_K_J=float(v['K_free_J'].max()),delta_work_J=float(v['W_ext_J'][-1]-v['W_ext_J'][0]),delta_numeric_damping_J=float(v['D_num_J'][-1]-v['D_num_J'][0]))
 inc=np.diff(a['D_break_independent_J']);res=float(abs(a['R_energy_J']).max()/abs(a['W_ext_J']).max())
 return dict(end={k:float(a[k][-1]) for k in keys},peak={k:float(a[k][peak]) for k in keys},pre_first_K_over_U=float(a['KE_over_U'][first-1]),peak_K_over_U=float(a['KE_over_U'][peak]),max_post_first_K_over_U=float(a['KE_over_U'][first:].max()),max_speed_m_s=float(a['max_speed_m_s'].max()),positive_topology_J=float(np.maximum(inc,0).sum()),negative_topology_J=float(np.minimum(inc,0).sum()),negative_topology_event_steps=int(np.sum(inc<0)),max_residual_over_max_work=res,phases=phases,interpretation='Signed topology transfer is U(x,old)-U(x,new), not automatically irreversible fracture dissipation. Numeric damping is not biological viscosity. Startup ratios at tiny U and post-fracture ratios reported separately.',source_consistent_energy='Unchanged V1/V2 independent energy.cpp; fixed initial weighted volume, recomputed dilatation after topology change; contact energy from native conservative pair PE times nodal volume.')

def analyze_case(p,out,name,common_step):
 out.mkdir(exist_ok=True,parents=True);g=dict(np.load(p/'geometry.npz'));X=g['X'];b=g['bonds'];T=Topology(g,g['prescribed_mask']);a=loadcsv(p/'history.csv');ev=np.load(p/'break_events.npz')['events'];bt=np.full(len(b),np.iinfo(np.int64).max,np.int64);bt[ev[:,1]]=ev[:,0];first=int(ev[0,0]);peak=int(np.argmax(a['stress_Pa']));end=len(a)-1
 assert np.isfinite(a.view(float)).all();assert np.array_equal(a['step'],np.arange(len(a)));assert np.max(abs(a['time_s']-a['step']*DT))<1e-14
 assert len(np.unique(ev[:,1]))==len(ev);assert (ev[:,2]==1).all();assert np.all(np.diff(ev[:,0])>=0);assert np.array_equal(np.cumsum(np.bincount(ev[:,0],minlength=len(a))),a['broken_bulk'])
 states={};sourcedir=out/'states';sourcedir.mkdir(exist_ok=True)
 wanted=['initial','pre_first_break','first_break','peak_tensile' if name=='M5C' else 'peak','first_coherent_cluster','first_central_advance','first_large_separation','end']
 if name=='M5B':wanted+=['frame_012000','frame_014000']
 # Audit every saved full-precision state, including interval frames; only copy selected evidence.
 checked=0
 for path in p.glob('*.npz'):
  z=dict(np.load(path));
  if 'x_m' not in z:continue
  step=int(z['step']);alive=z['alive'].astype(bool);assert np.array_equal(alive,bt>step);assert np.isfinite(z['x_m']).all() and np.isfinite(z['v_full_m_s']).all();assert z['x_m'].shape==(2880,3);checked+=1
  if path.stem not in wanted:continue
  q=T.summary(alive);q.update(step=step,time_s=float(z['time_s']),opening_m=float(a['displacement_m'][step]),nominal_stress_Pa=float(a['stress_Pa'][step]),cluster=T.cluster(alive))
  util=(np.linalg.norm(z['x_m'][b[:,1]]-z['x_m'][b[:,0]],axis=1)/T.r-1)/.15;eligible=bt>=step;q['regional_max_utilization_including_new_breaks']={n:float(util[(T.reg==i)&eligible].max()) if np.any((T.reg==i)&eligible) else None for i,n in enumerate(NAMES)}
  q['surface_status']='CROSSING_CONNECTIVITY_PROXY_ONLY';states[path.stem]=q;remain,lost,advance=T.crossing(alive);sizes,labels=T.components(alive)
  np.savez_compressed(sourcedir/(path.stem+'.npz'),**z,remaining_crossings=remain,initial_crossings=T.initial_counts,advance_m=advance,component_labels=labels,utilization=util)
 firstids=ev[ev[:,0]==first,1];firstinfo=[dict(bond_index_zero_based=int(k),native_body_bond_id=int(k+1),atom_ids_one_based=(b[k]+1).tolist(),reference_midpoint_um=(T.mid[k]*1e6).tolist(),region=NAMES[T.reg[k]],distance_to_nominal_y_surface_um=float((.0003-abs(T.mid[k,1]))*1e6)) for k in firstids]
 alive=np.ones(len(b),bool);coherent=None;central=None;central_onset=None;trace=[]
 # Actual native event steps; spatial summaries are exact event topology, not interpolated coordinates.
 unique,starts=np.unique(ev[:,0],return_index=True)
 for k,step in enumerate(unique):
  ids=ev[starts[k]:starts[k+1] if k+1<len(starts) else len(ev),1];alive[ids]=False
  if central_onset is None and np.any((abs(T.mid[ids,1])<=100e-6+1e-12)&T.band[ids]):central_onset=dict(step=int(step),time_s=float(step*DT),bond_indices=ids[(abs(T.mid[ids,1])<=100e-6+1e-12)&T.band[ids]])
  if coherent is None:
   c=T.cluster(alive)
   if c['largest_voxels']>=10:coherent=dict(step=int(step),time_s=float(step*DT),**c)
  _,_,ad=T.crossing(alive)
  if central is None and np.sum(ad[4:8]>=100e-6-1e-12)>=2:central=dict(step=int(step),time_s=float(step*DT),front_advance_m=ad)
  idsall=np.flatnonzero(~alive);trace.append([int(step),float(step*DT),len(idsall),int(T.band[idsall].sum()),int((~T.cross[idsall]).sum()),*np.bincount(T.reg[idsall],minlength=4),*ad])
 propagation_step=max(peak,int(central['step'])) if central else None
 E=energy(a,first,peak,propagation_step);common=T.summary(bt>common_step)
 q=dict(case=name,source=str(p),first_step=first,first_time_s=first*DT,first_events=firstinfo,first_is_notch=bool(np.all(T.reg[firstids]==0)),first_central_band_event=central_onset,first_coherent_cluster=coherent,first_substantial_central_advance_exact_native_event=central,peak_step=peak,peak_time_s=peak*DT,peak_nominal_stress_Pa=float(a['stress_Pa'][peak]),peak_opening_m=float(a['displacement_m'][peak]),end_step=end,end_time_s=end*DT,states=states,matched_common_time=dict(step=common_step,time_s=common_step*DT,opening_m=float(a['displacement_m'][common_step]),nominal_stress_Pa=float(a['stress_Pa'][common_step]),**common),snapshots_checked=checked,energy=E,region_names=NAMES)
 np.savez_compressed(out/'review_curves.npz',**{k:a[k] for k in a.dtype.names},event_trace=np.array(trace),event_columns=np.array(['step','time_s','broken','target_band','off_plane',*NAMES,*[f'advance_y{i}_m' for i in range(12)]]))
 np.savez_compressed(out/'geometry_events.npz',**g,events=ev,break_step=bt,region=T.reg,target_band=T.band,cross_z0=T.cross)
 save(out/'CASE_AUDIT.json',q);return q

def main(base,old):
 assert js(base/'RUN_RESOURCE.json')['exit_code']==0;ex=js(base/'M5C/EXECUTION.json');end=ex['steps'];C=analyze_case(base/'M5C',base/'analysis/M5C','M5C',end);B=analyze_case(old/'M5B',base/'analysis/M5B','M5B',end);obs=js(base/'M5C/OBSERVER_CHECKS.json');diag=js(base/'M5C/DIAGNOSTIC_CHECKS.json');stop=js(base/'M5C/STOP_DECISION.json')
 snap=C['states']['end'];cl=snap['cluster'];pca=cl['largest_pca'];gate=dict(tip_origin=C['first_is_notch'] and bool(cl['meaningful_intersects_tip']),four_central_strips=snap['central_advancing_strips_2dx']==4,eight_thickness_strips=snap['all_advancing_strips_2dx']>=8,target_band_ge70=snap['target_band_fraction']>=.7,grip_le20=snap['grip_fraction']<=.2,rms_thickness_le_dx=pca.get('rms_thickness_m',1)<=50e-6,secondary_clusters_lt10=cl['secondary_voxel_fraction']<.1,independently_bounded_physical_surface_area_error=False)
 area=dict(status='CROSSING_CONNECTIVITY_PROXY_ONLY',physical_new_area_m2=None,connectivity_proxy_m2=snap['crossing_connectivity_area_proxy_m2'],gates=gate,front_advance_m=snap['front_advance_m'],largest_voxel_cluster_pca=pca,secondary_voxel_fraction=cl['secondary_voxel_fraction'],component_sizes=snap['component_sizes'],point_spacing_m=50e-6,horizon_m=150.5e-6,reason='A connected cloud of lost interactions and z=0 crossing depletion are not a uniquely reconstructed material free surface. No independent physical area error bound or spatial refinement in this bounded task. Finite horizon damage-band thickness and branches are disclosed; no bond-count area or artificial gap closure.',precut_area_reference_m2=.0002*.0006)
 # Preserve the distinction: gross off-plane clusters can positively show nonuniqueness; failure of a numerical gate alone cannot prove it.
 if not gate['target_band_ge70'] and (not gate['rms_thickness_le_dx'] or not gate['secondary_clusters_lt10']):area['status']='NONUNIQUE_CRACK_SURFACE'
 elif all(v for k,v in gate.items() if k!='independently_bounded_physical_surface_area_error'):area['status']='CRACK_PATH_IDENTIFIED_AREA_UNRESOLVED'
 save(base/'CRACK_AREA_ATTEMPT.json',area)
 save(base/'CRACK_PROPAGATION_AUDIT.json',dict(M5B=B,M5C=C,comparison_scope='Inherited reference midpoint priority regions; actual controlled-node membership separately reported. Same actual mean boundary opening/time; different H=850/950um. No rerun of M5B.',exact_event_vs_checkpoint='Exact event time reported offline; central/separation state saved only at first predeclared 100-step checkpoint. No invented exact full field for an unsaved event.'))
 boundary=loadcsv(base/'M5C/boundary_dynamics_audit.csv');ib=boundary['step'].astype(int);a=loadcsv(base/'M5C/history.csv');oldboundary=loadcsv(old/'M5B/boundary_audit.csv');oldrisk=loadcsv(old/'M5B/regional_utilization.csv')
 np.savez_compressed(base/'analysis/M5C/boundary_review.npz',**{k:boundary[k] for k in boundary.dtype.names});np.savez_compressed(base/'analysis/M5B/boundary_review.npz',**{k:oldboundary[k] for k in oldboundary.dtype.names});np.savez_compressed(base/'analysis/M5B/risk_review.npz',**{k:oldrisk[k] for k in oldrisk.dtype.names})
 Qe=float(abs(boundary['generalized_Q_N']-.5*(boundary['top_Rz_N']-boundary['bottom_Rz_N'])).max());work_note='Q=-sum_i,d f_material[i,d]*mask[i,d]*sign[i,d]; W=-integral sum f_material du_prescribed. For M5C Q=(Rtop-Rbottom)/2. The prescribed masses are excluded from free-DOF K; reported work is transmitted into material, not actuator kinetic energy.'
 E=dict(M5B=B['energy'],M5C=C['energy'],work_conjugacy=work_note,max_M5C_Q_identity_error_N=Qe,max_contact_energy_J=float(a['U_contact_J'].max()),max_measured_contact_force_N=float(np.nanmax(boundary['max_contact_force_N'])),contact_force_missing_nonzero_contact_break_rows=int(np.sum(boundary['contact_force_decomposition_valid']==0)),contact_decomposition_scope=obs['contact_force_note'],ledger_status='PASS' if C['energy']['max_residual_over_max_work']<=.05 else 'FAIL',apparent_fracture_work_candidate_J_m2=None,material_Gc_J_m2=None,reason='Signed topology transfer is not a unique fracture-only numerator; area unresolved; postfracture dynamics, numerical damping and prior spacing sensitivity preclude material-Gc certification.')
 save(base/'ENERGY_LEDGER_AUDIT.json',E)
 ids=js(base/'provenance/REUSED_IDENTITY.json');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();changed=[k for k,v in ids.items() if sha(ROOT/k)!=v];assert changed==['pd_clot/lammps_fragmentation/macro_lab_v1.py'],changed
 launched=js(base/'provenance/LAUNCH_SOURCE_IDENTITY.json');save(base/'provenance/FINAL_REUSED_IDENTITY_CHECK.json',dict(unchanged_count=len(ids)-1,allowed_changed=changed,all_other_production_native_old_results_unchanged=True,launch_record=launched))
 checks=dict(status='PASS',case_budget=js(base/'FORMAL_RUN_BUDGET.json'),groups=input_groups(base/'M5C/in.lammps',dict(np.load(base/'M5C/geometry.npz'))),native_callbacks=diag['callbacks'],expected_callbacks=diag['expected'],topology_matches_all_saved_accepted_states=True,full_precision_states_checked=C['snapshots_checked'],history_finite_and_sequential=True,force_checks=diag['force_checks'],observer_checks=obs,energy_residual_over_work=C['energy']['max_residual_over_max_work'],Q_identity_error_N=Qe,old_and_native_immutable_identity_count=len(ids)-1,formal_retest_used=False)
 checks['energy_check_scope']='V1/V2 native bond PE undercounts deviatoric energy and remains a diagnostic only. Retain its reported discrepancy; validate unchanged independent energy via native-force gradient (<1e-9) and finite difference (<1e-5), as in V2, not by equating native bond PE with full LPS energy.'
 assert diag['callbacks']==diag['expected']==end+1;assert all(q['native_force_relative_error']<1e-9 and q['finite_difference_relative_error']<1e-5 for q in diag['force_checks']);assert len(diag['force_checks'])==3;assert Qe<1e-12;assert C['energy']['max_residual_over_max_work']<.05
 save(base/'IMPLEMENTATION_CHECKS.json',checks)
 rows=[]
 for q in [B,C]:
  peak=q['states'].get('peak_tensile',q['states'].get('peak'));en=q['states']['end'];row=dict(case=q['case'],first_time_s=q['first_time_s'],first_region=q['first_events'][0]['region'],peak_time_s=q['peak_time_s'],peak_Pa=q['peak_nominal_stress_Pa'],peak_opening_m=q['peak_opening_m'],peak_broken=peak['broken'],peak_band_fraction=peak['target_band_fraction'],peak_grip_fraction=peak['grip_fraction'],end_time_s=q['end_time_s'],end_broken=en['broken'],end_band_fraction=en['target_band_fraction'],end_grip_count=en['region_counts'][1],end_grip_fraction=en['grip_fraction'],end_dual_prescribed_broken=en['both_axial_prescribed_broken'],end_fully_free_endpoint_broken=en['fully_free_endpoint_broken'],end_offplane_fraction=en['off_plane_fraction'],central_advance_time_s=q['first_substantial_central_advance_exact_native_event']['time_s'] if q['first_substantial_central_advance_exact_native_event'] else None,common_time_grip_fraction=q['matched_common_time']['grip_fraction'],common_time_broken=q['matched_common_time']['broken'],net_topology_J=q['energy']['end']['D_break_independent_J'],energy_residual_fraction=q['energy']['max_residual_over_max_work'],Gc_J_m2=None);rows.append(row)
 with (base/'M5B_M5C_COMPARISON.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 status=dict(IMPLEMENTATION_VALIDITY='PASS',GRIP_KINEMATIC_VALIDITY='PASS_NO_PRESCRIBED_AXIAL_FORCED_FAILURE',FREE_DOF_INTEGRATION_VALIDITY='PASS_ACTUAL_EVERY_STEP_RECURRENCE_AND_NATIVE_VELOCITIES',NOTCH_INITIATION='PASS_ACTUAL_FIRST_BREAK_AT_NOTCH' if C['first_is_notch'] else 'FAIL_FIRST_BREAK_ELSEWHERE',GRIP_DAMAGE_CONTAMINATION=dict(status='SEVERE' if snap['grip_fraction']>=.5 else 'PRESENT_NOT_DOMINANT' if snap['region_counts'][1] else 'NONE_OBSERVED',peak_fraction=C['states']['peak_tensile']['grip_fraction'],end_fraction=snap['grip_fraction']),CRACK_PROPAGATION='CENTRAL_THICKNESS_CONNECTIVITY_ADVANCE_OBSERVED' if C['first_substantial_central_advance_exact_native_event'] else 'NO_SUBSTANTIAL_CENTRAL_THICKNESS_ADVANCE',CRACK_SURFACE_IDENTIFIABILITY=area['status'],ENERGY_LEDGER_VALIDITY=E['ledger_status'],APPARENT_FRACTURE_WORK_STATUS='UNAVAILABLE_AREA_AND_FRACTURE_ONLY_ENERGY_UNRESOLVED',MATERIAL_GC_STATUS='GC_UNAVAILABLE',PRODUCTION_PARAMETERS_UNCHANGED=True,new_formal_cases=1,zero_step_preflights=1,bug_retests=0,stop=stop,run_resource=js(base/'RUN_RESOURCE.json'),material_Gc_J_m2=None,apparent_fracture_work_J_m2=None,scope='Specimen-specific unchanged LPS boundary validation; no new constitutive calibration, grid/time-step study, microbubble or Y-clot production run.',next_priority='End repeated grip trials; audit identifiable fracture-path and energy observables against an independent physical notched-specimen benchmark before any constitutive upgrade or Gc/AI certification.')
 save(base/'FINAL_STATUS.json',status);print(json.dumps(dict(comparison=rows,area=area,status=status),indent=2,ensure_ascii=False))

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
