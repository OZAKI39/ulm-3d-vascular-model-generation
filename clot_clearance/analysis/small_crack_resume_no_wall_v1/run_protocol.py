from integrator import *
from journal import Journal
from mechanics import Mechanics
from scipy import ndimage
import resource,traceback

def concentration(ref,d):
 damaged=(d>0)&~ref['initial_cut'];count=int(damaged.sum());tipfrac=float(np.mean(ref['bond_regions'][damaged]<2)) if count else 0.
 if count<20 or tipfrac<.6:return None
 mid=.5*(ref['X'][ref['bonds'][damaged,0]]+ref['X'][ref['bonds'][damaged,1]]);ij=np.clip(np.floor((mid+SIZE/2)/DX).astype(int),0,np.array(COUNTS)-1);grid=np.zeros(COUNTS,bool);grid[tuple(ij.T)]=True;labels,n=ndimage.label(grid,structure=np.ones((3,3,3)))
 for k in range(1,n+1):
  pts=np.argwhere(labels==k);x=(pts+.5)*DX-SIZE/2
  if len(pts)>=6 and np.ptp(x,axis=0).max()>=2*DX-1e-15 and np.min(np.hypot(x[:,0]-TIP,x[:,2]))<=DX:
   return dict(damaged_bonds=count,notch_fraction=tipfrac,connected_voxels=len(pts),reference_voxel_centers_m=x,interpretation='Localized damage cloud, not a proved continuous crack surface')
 return None

def main():
 freeze=json.loads((O/'SMALL_CRACK_DESIGN_FREEZE.json').read_text());assert freeze['budget']['main_wall_s'] is None and freeze['budget']['RSS_stop_GiB'] is None
 for n,h in freeze['driver_sha256'].items():assert sha(A/n)==h,n
 assert sha(LIB)==freeze['library_sha256']
 with (O/'full/FORMAL_CASE_STARTED.json').open('x') as f:json.dump(dict(pid=os.getpid(),fresh_t0=False,continuation_of=str(PARENT),start_step=9600),f)
 start=time.perf_counter();cpu0=time.process_time();cp=freeze['continuation']['checkpoint']
 for key in ['restart','state']:assert sha(cp[key])==cp[key+'_sha256']
 original=dict(np.load(PARENT/'full/keyframes/step_000000.npz'));saved=dict(np.load(cp['state']));ref=np.load(O/'provenance/reference.npz');ids=np.load(O/'provenance/boundary_ids.npz')
 import csv
 with (PARENT/'NATIVE_DYNAMIC_HISTORY.csv').open() as f:
  reader=csv.DictReader(f);columns=reader.fieldnames;rows=[{k:float(v) for k,v in row.items()} for row in reader]
 lastrow=rows[-1];assert len(rows)==9601 and int(lastrow['step'])==9600
 ns=Native('full/native_resume',restart=cp['restart'],dt=lastrow['dt_s']);ns.dynamics(ids['bottom_ids'],ids['top_ids'],ids['free_ids']);s=ns.state()
 restored=compare(saved,s);codes,special=topology(ns);expected=ref['bonds'][s['d']<1]+1
 restored.update(original_bond_ids_exact=bool(np.array_equal(ns.field(9),ref['bond_ids'])),reference_m_exact=bool(np.array_equal(ns.field(4),ref['m'])),native_endpoint_graph_exact=bool(np.array_equal(codes,np.sort(expected[:,0]*len(ref['X'])+expected[:,1]))),special_exclusions_sha256=special,physical_steps_advanced=0,source_checkpoint=cp)
 dump(O/'NATIVE_RESTART_RESTORE_AUDIT.json',restored)
 assert restored['status']=='PASS' and restored['original_bond_ids_exact'] and restored['reference_m_exact'] and restored['native_endpoint_graph_exact'],restored
 np.savez_compressed(O/'RESTORED_NATIVE_STATE.npz',**s)
 mech=Mechanics(ns,ids,original,freeze['gates']);mech.W=lastrow['W_J'];mech.Jc=lastrow['Jc_J'];mech.cwork=lastrow['contact_work_J'];mech.scale=lastrow['energy_scale_J'];mech.defects=np.array([lastrow[k] for k in ['VV_defect_J','spatial_defect_J','switch_work_J','contact_defect_J']]);mech.Dnode=np.load(PARENT/'full/accepted_node_D.npz')['node_D_J'].copy();mech.maxerr=json.loads((PARENT/'full/DYNAMIC_AUDIT.json').read_text())['max_accepted_errors'];mech.last=mech.model.evaluate(s['x'],s['d']);mech.prevfc=mech.contact.evaluate(s['x'],s['d'])[1];mech.update_components(s)
 restored_R=mech.W-(s['scalars'][0]-original['scalars'][0]+s['K']-original['K']+s['scalars'][1]+s['Uc']-original['Uc']-mech.Jc)
 assert abs(restored_R-lastrow['R_J'])<=freeze['gates']['energy_abs_J'] and abs(mech.Dnode.sum()-s['scalars'][1])<=freeze['gates']['energy_abs_J']
 dump(O/'LEDGER_RESTORE_AUDIT.json',dict(status='PASS',global_t0_energy_origin=True,cumulative_W_J=mech.W,cumulative_D_J=float(s['scalars'][1]),cumulative_Jc_J=mech.Jc,cumulative_contact_work_J=mech.cwork,energy_scale_J=mech.scale,defects_J=mech.defects,restored_R_J=float(restored_R),source_R_J=lastrow['R_J'],R_jump_J=float(restored_R-lastrow['R_J']),node_D_sum_error_J=float(mech.Dnode.sum()-s['scalars'][1]),previous_max_errors_retained=True))
 journal=Journal(columns,ref,resume=True);journal.last=s;journal.last_row=lastrow;journal.records=len(rows);journal.checkpoint(ns,s,'verified_resume_2p4ms');journal.event('NATIVE_RESTART_RESUMED',source_run=str(PARENT),step=s['step'],time_s=s['time'],physical_steps_advanced=0)
 engine=Integrator(ns,ids,journal);oldtrial=json.loads((PARENT/'full/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());engine.trials=oldtrial['trials'];engine.rollbacks=oldtrial['rollbacks'];engine.serial=max(z['trial'] for z in engine.trials);engine.accepts=oldtrial['accepted_damage_gate'];engine.rejects=oldtrial['rejected'];engine.base_intervals=oldtrial['base_intervals']
 prior_events=json.loads((PARENT/'full/FIRST_EVENTS.json').read_text());events=prior_events['global_events'];first_region_soft=prior_events['first_region_softening'];first_region_fail=prior_events['first_region_failure'];first_concentration=prior_events['first_observed_concentration'];first_open=prior_events['first_open_no_bridge_cell'];first_contact=prior_events['first_contact']
 prev=s;last_Dnode=mech.Dnode.copy();accepted_maxerr=dict(mech.maxerr);stop=[];native_topology=False
 print('RESUME_VERIFIED',s['step'],s['time'],'ledger_R_jump',restored_R-lastrow['R_J'],flush=True)
 try:
  while prev['time']<END-1e-15:
   s,h=engine.next(prev);ns=engine.ns;mech.ns=ns;row,pairs,newsoft,newfail=mech.check(prev,s,h)
   journal.accept(s,row,prev=prev,pairs=pairs,newsoft=newsoft,newfail=newfail);rows.append(row);prev=s;last_Dnode=mech.Dnode.copy();accepted_maxerr=dict(mech.maxerr)
   for kind,bondids,seen in [('softening',newsoft,first_region_soft),('failure',newfail,first_region_fail)]:
    if len(bondids):
     if kind not in events:
      events[kind]=dict(step=s['step'],time_s=s['time'],bond_ids=bondids+1,reference_midpoints_m=.5*(ref['X'][ref['bonds'][bondids,0]]+ref['X'][ref['bonds'][bondids,1]]),actual_midpoints_m=.5*(s['x'][ref['bonds'][bondids,0]]+s['x'][ref['bonds'][bondids,1]]),regions=[REGIONS[int(j)] for j in ref['bond_regions'][bondids]],kappa=s['kappa'][bondids],d=s['d'][bondids],stretch=np.linalg.norm(s['x'][ref['bonds'][bondids,1]]-s['x'][ref['bonds'][bondids,0]],axis=1)/ref['r0'][bondids]-1);journal.key(s,'first_'+kind)
     for reg in np.unique(ref['bond_regions'][bondids]):
      name=REGIONS[int(reg)]
      if name not in seen:seen[name]=dict(time_s=s['time'],step=s['step'],bond_ids=bondids[ref['bond_regions'][bondids]==reg]+1);journal.key(s,'first_'+kind+'_'+name)
   if first_contact is None and len(pairs):first_contact=dict(step=s['step'],time_s=s['time'],pairs=pairs+1);journal.key(s,'first_contact')
   if s['step']%20==0:
    if first_concentration is None:
     conc=concentration(ref,s['d'])
     if conc is not None:first_concentration=dict(step=s['step'],time_s=s['time'],**conc);journal.key(s,'first_observed_concentration')
    if first_open is None:
     counts=np.bincount(ref['cross_cells'][s['d'][ref['cross_ids']]<1],minlength=COUNTS[0]*COUNTS[1]);gap=s['x'][ref['face_upper'],2]-s['x'][ref['face_lower'],2]-DX;xx=ref['X'][ref['face_upper'],0];candidate=(xx>TIP)&(counts==0)&(gap>=.1*DX)
     if np.any(candidate):first_open=dict(step=s['step'],time_s=s['time'],candidate_cells=np.flatnonzero(candidate),label='FIRST_OPEN_NO_BRIDGE_CELL_NOT_YET_AREA_CERTIFIED');journal.key(s,'first_open_no_bridge_cell')
   if any(abs(s['time']-phase)<1e-15 for phase in [.0005,.001,.0015,.002,.0025]):journal.key(s,'fixed_observation_time')
   if s['step']%400==0:journal.checkpoint(ns,s,'periodic_400_accepted_steps')
   if s['step']%200==0:
    elapsed=time.perf_counter()-start;rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024;bytes_now=sum(p.stat().st_size for p in (O/'full').rglob('*') if p.is_file());print('STEP',s['step'],'t_ms',s['time']*1e3,'new_failed_total',row['total_natural_failed'],'components',row['components'],'max_dd',row['max_dd'],'relR',row['relative_R'],'wall_s',elapsed,flush=True)
    assert freeze['budget']['RSS_stop_GiB'] is None and freeze['budget']['main_wall_s'] is None,'USER_AUTHORIZED_NO_MEMORY_OR_WALL_CAP'
    assert bytes_now<=freeze['budget']['output_limit_GiB']*1024**3,'OUTPUT_BUDGET'
 except Exception as exc:stop.append(repr(exc));journal.failure(exc,traceback.format_exc())
 finally:
  engine.save();ns=engine.ns;s=journal.last
  if not stop:
   try:
    codes,special_hash=topology(ns);b=ref['bonds'][s['d']<1]+1;native_topology=np.array_equal(codes,np.sort(b[:,0]*ns.N+b[:,1]));assert native_topology
    assert np.array_equal(ns.field(4),ref['m']);np.savez_compressed(O/'full/final_native_topology.npz',atom_ids=ns.gather('id'),endpoint_codes=codes,original_bond_ids=ns.field(9),failed=ns.field(10),reference_m=ns.field(4));dump(O/'full/FINAL_SPECIAL_EXCLUSIONS_HASH.json',dict(sha256=special_hash));journal.checkpoint(ns,s,'final_accepted');journal.key(s,'final_accepted')
   except Exception as exc:stop.append('FINALIZATION '+repr(exc));journal.failure(exc,traceback.format_exc())
  else:journal.key(s,'last_accepted_after_stop')
  journal.finish();dump(O/'full/FIRST_EVENTS.json',dict(global_events=events,first_region_softening=first_region_soft,first_region_failure=first_region_fail,first_observed_concentration=first_concentration,first_open_no_bridge_cell=first_open,first_contact=first_contact))
  np.savez_compressed(O/'full/accepted_node_D.npz',node_D_J=last_Dnode,node_regions=ref['node_regions'])
  dump(O/'RUN_RESOURCES.json',dict(main_wall_s=time.perf_counter()-start,main_CPU_s=time.process_time()-cpu0,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,formal_dynamic_cases=0,continuation_segments=1,new_steps=len(rows)-1-9600,accepted_steps=len(rows)-1,damage_rejections=engine.rejects,checkpoints_written=journal.checkpoints,checkpoint_wall_s=journal.checkpoint_wall_s,native_candidate_wall_s=engine.native_wall,rollback_wall_s=engine.restore_wall,MPI_ranks=1,threads=1,GPU_used=False,output_bytes_at_main_end=sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),budget=freeze['budget']))
  maps=sorted(set(q.split()[-1] for q in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in q));dump(O/'full/DYNAMIC_AUDIT.json',dict(status='COMPLETED' if not stop else 'FAILED_ADMISSION',stop_reasons=stop,from_original_precut_t0=False,continuous_original_protocol=True,source_run=str(PARENT),resume_step=9600,last_accepted_step=s['step'],last_accepted_time_s=s['time'],planned_end_time_s=END,accepted_steps=len(rows)-1,damage_rejections=engine.rejects,last_accepted_ledger=journal.last_row,max_accepted_errors=accepted_maxerr,native_final_topology_matches=bool(native_topology),original_atom_ids_checked_every_step=True,independent_mechanics_every_step=True,actual_loaded_libraries={p:sha(p) for p in maps},next_case_executed=False))
  try:ns.close()
  except Exception as exc:dump(O/'full/CLOSE_WARNING.json',dict(exception=repr(exc)))
 print('FORMAL_DONE',json.dumps(dict(status='COMPLETED' if not stop else 'FAILED_ADMISSION',time_s=float(s['time']),accepted=len(rows)-1,rejected=engine.rejects,stop_reasons=stop)),flush=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:
  dump(O/'SETUP_OR_FINALIZATION_FAILURE.json',dict(exception=repr(exc),traceback=traceback.format_exc()));raise
