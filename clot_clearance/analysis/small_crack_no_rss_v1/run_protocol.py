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
 freeze=json.loads((O/'SMALL_CRACK_DESIGN_FREEZE.json').read_text());assert json.loads((O/'MINIMAL_IMPLEMENTATION_CHECK.json').read_text())['status']=='PASS'
 for n,h in freeze['driver_sha256'].items():assert sha(A/n)==h,n
 assert sha(LIB)==freeze['library_sha256']
 with (O/'full/FORMAL_CASE_STARTED.json').open('x') as f:json.dump(dict(pid=os.getpid(),fresh_t0=True,formal_case=1),f)
 start=time.perf_counter();cpu0=time.process_time();ns=Native('full/native',dt=BASE);ref=ns.ref;ids=np.load(O/'provenance/boundary_ids.npz');ns.dynamics(ids['bottom_ids'],ids['top_ids'],ids['free_ids']);s=ns.state();assert s['step']==0 and s['time']==0 and s['scalars'][1]==s['scalars'][3]==0 and np.array_equal(s['d'],ref['initial_cut'].astype(float))
 mech=Mechanics(ns,ids,s,freeze['gates']);columns=['step','time_s','dt_s','opening_m','reaction_top_N','reaction_bottom_N','generalized_material_reaction_N','actuator_top_N','W_J','U_J','K_J','Kfree_J','D_J','Uc_J','Jc_J','R_J','relative_R','energy_scale_J','max_speed_m_s','max_displacement_m','nonaffine_RMS_m','step_disp_dx','relative_endpoint_step_dx','max_dd','max_kappa','new_soft','new_fail','total_natural_failed','initial_cut_count','live_bonds','components','contact_pairs','max_contact_penetration_dx','relative_contact_step_dx','atom_count','Px_Ns','Py_Ns','Pz_Ns','momentum_residual_Ns','x_equation_m','v_equation_m_s','contact_work_J','VV_defect_J','spatial_defect_J','switch_work_J','contact_defect_J','remainder_J','min_box_clearance_m','commits','token']
 for name in REGIONS:columns+=['D_'+name+'_J','failed_'+name]
 for axis in 'xyz':columns+=[axis+'min_m',axis+'max_m']
 row={k:0. for k in columns};row.update(step=0,time_s=0,dt_s=BASE,U_J=float(s['scalars'][0]),K_J=s['K'],Uc_J=s['Uc'],energy_scale_J=1e-15,initial_cut_count=int(ref['initial_cut'].sum()),live_bonds=s['native_bonds'],components=mech.components,atom_count=ns.N,max_kappa=SC,min_box_clearance_m=float(np.min(BOX-abs(s['x']))),commits=0,token=-1)
 for j,axis in enumerate('xyz'):row[axis+'min_m']=float(s['x'][:,j].min());row[axis+'max_m']=float(s['x'][:,j].max())
 journal=Journal(columns,ref);journal.accept(s,row);journal.key(s,'initial');journal.checkpoint(ns,s,'initial_reference_with_precut');engine=Integrator(ns,ids,journal);prev=s;rows=[row];events={};first_region_soft={};first_region_fail={};last_Dnode=mech.Dnode.copy();accepted_maxerr={};stop=[];native_topology=False;first_concentration=None;first_open=None;first_contact=None
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
    assert freeze['budget']['RSS_stop_GiB'] is None,'USER_AUTHORIZED_NO_RSS_LIMIT'
    assert elapsed<=freeze['budget']['main_wall_s'] and bytes_now<=freeze['budget']['output_limit_GiB']*1024**3,'WALL_OR_OUTPUT_BUDGET'
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
  dump(O/'RUN_RESOURCES.json',dict(main_wall_s=time.perf_counter()-start,main_CPU_s=time.process_time()-cpu0,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,formal_dynamic_cases=1,accepted_steps=len(rows)-1,damage_rejections=engine.rejects,checkpoints_written=journal.checkpoints,checkpoint_wall_s=journal.checkpoint_wall_s,native_candidate_wall_s=engine.native_wall,rollback_wall_s=engine.restore_wall,MPI_ranks=1,threads=1,GPU_used=False,output_bytes_at_main_end=sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),budget=freeze['budget']))
  maps=sorted(set(q.split()[-1] for q in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in q));dump(O/'full/DYNAMIC_AUDIT.json',dict(status='COMPLETED' if not stop else 'FAILED_ADMISSION',stop_reasons=stop,from_original_precut_t0=True,last_accepted_step=s['step'],last_accepted_time_s=s['time'],planned_end_time_s=END,accepted_steps=len(rows)-1,damage_rejections=engine.rejects,last_accepted_ledger=journal.last_row,max_accepted_errors=accepted_maxerr,native_final_topology_matches=bool(native_topology),original_atom_ids_checked_every_step=True,independent_mechanics_every_step=True,actual_loaded_libraries={p:sha(p) for p in maps},next_case_executed=False))
  try:ns.close()
  except Exception as exc:dump(O/'full/CLOSE_WARNING.json',dict(exception=repr(exc)))
 print('FORMAL_DONE',json.dumps(dict(status='COMPLETED' if not stop else 'FAILED_ADMISSION',time_s=float(s['time']),accepted=len(rows)-1,rejected=engine.rejects,stop_reasons=stop)),flush=True)
if __name__=='__main__':main()
