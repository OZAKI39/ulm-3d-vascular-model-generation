from common import *
from native import native,predict
from current_surface import CurrentSurface
from multifragment_shear_objective import scalar_metrics,shape,selection_key,SPEC
import time,resource,traceback,csv,gzip,argparse,ctypes
graph=module('crushing_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
def save(case,name,s):
 p=case/'keys'/name;tmp=p.with_suffix('.tmp.npz');np.savez_compressed(tmp,**s);os.replace(tmp,p)
def run(name):
 start=time.perf_counter();cpu=time.process_time();case=O/'cases'/name
 assert not (case/'history.csv').exists(),'Preserve interrupted attempts; use a new output name for a justified recovery'
 ns,C=native(case,name);ref=ns.ref;n=ns.n;N=ns.N;b=ref['bonds'];bt=ref['bond_types'];cut=ref['initial_cut'];model=energy_model(ref);contact=C(ref);on=ns.field(13);off=ns.field(14)
 assert json.loads((O/'KEY_CHECKS.json').read_text())['status']=='PASS'
 freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
 for path,want in freeze['implementation_sha256'].items():assert sha(R/path)==want,path
 for path,want in freeze['input_sha256'].items():assert sha(O/path)==want,path
 dump(case/'REQUEST.json',dict(name=name,uniform_material=(name=='uniform'),label=ns.tool.label,library_sha256=sha(LIB),source_sha256={p.name:sha(p) for p in A.glob('*.py')},input_sha256={name:sha(O/name) for name in ['SOURCE_FREEZE.json','FROZEN_SHEAR_DRIVE.json','BIPHASIC_MATERIAL_SPEC.json','CURRENT_EXPOSED_SURFACE_SPEC.json']}))
 s=ns.dynamics();ns.tool.commit(ns.tool.evaluate(s['x'],s['velocity'],s['d'],s['time']));initial=s;prev=s;save(case,'initial.npz',s);W=0.;Wdrive=0.;Wdrag=0.;Jc=0.;impulse=np.zeros(3);rows=[];events={};err={};stop=[];substeps=0;rejects=0;scale=1e-24;secondary=[];maxscore=(-1,);snap_next=25e-6;absolute_impulse=0.;rotation_impulse=np.zeros(3);geometry_first=None;surface_events={};stage_work=0.;new_surface_work=0.;facet_impulse=0.
 cached=graph.topology(ref,s['d']);registry={};lookup={};current=[];previous_free={}
 history=gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','wt');changes=gzip.open(case/'damage_events.jsonl.gz','wt');stream=(case/'history.csv').open('w');writer=None
 def error(key,value):err[key]=max(err.get(key,0.),float(value))
 def event(name,ids,state):
  if name in events or not len(ids):return
  ends=b[ids];mid=ref['X'][ends].mean(1)
  events[name]=dict(time_s=state['time'],step=state['step'],original_bond_ids=ids+1,atom_ids=ends+1,bond_types=bt[ids],reference_midpoints_m=mid,actual_midpoints_m=state['x'][ends].mean(1),stretch=np.linalg.norm(state['x'][ends[:,1]]-state['x'][ends[:,0]],axis=1)/ref['r0'][ids]-1,kappa=state['kappa'][ids],damage=state['d'][ids],s0=on[ids],near_anchor=(mid[:,2]-ref['X'][:n,2].min())<=HORIZON,evidence='keys/'+name+'.npz');save(case,name+'.npz',state)
 def component_row(state,topology_changed):
  nonlocal current,previous_free
  lab,sizes,attached=cached;cs=[]
  if topology_changed:
   current=[]
   for label in np.flatnonzero(~attached):
    ids=np.flatnonzero(lab==label);key=tuple((ids+1).tolist())
    if key not in lookup:
     cid=len(registry)+1;lookup[key]=cid;registry[cid]=dict(atom_ids=ids+1,nodes=len(ids),reference_volume_m3=len(ids)*V,equivalent_diameter_m=(6*len(ids)*V/np.pi)**(1/3),birth_time_s=state['time'],birth_step=state['step'])
    cid=lookup[key];current.append((label,cid,ids))
   for old,oldset in previous_free.items():
    children=[cid for _,cid,ids in current if oldset.intersection((ids+1).tolist())]
    if len(children)>1:
     secondary.append(dict(time_s=state['time'],step=state['step'],parent_component_id=old,children_component_ids=children))
     if len(secondary)==1:save(case,'first_secondary.npz',state)
   previous_free={cid:set((ids+1).tolist()) for _,cid,ids in current}
  bulk=(bt==1)&~cut;bb=b[bulk];same=lab[bb[:,0]]==lab[bb[:,1]];ilab=lab[bb[same,0]];dam=state['d'][bulk][same];counts=np.bincount(ilab,minlength=len(sizes));alive=np.bincount(ilab,weights=dam<1,minlength=len(sizes));dsum=np.bincount(ilab,weights=dam,minlength=len(sizes))
  for label,cid,ids in current:
   v=state['velocity'][ids];vc=v.mean(0);relativeKE=.5*RHO*V*np.sum((v-vc)**2)
   cs.append(dict(component_id=cid,internal_initial_bonds=int(counts[label]),internal_alive_bonds=int(alive[label]),internal_survival=float(alive[label]/counts[label]) if counts[label] else None,mean_internal_damage=float(dsum[label]/counts[label]) if counts[label] else None,centroid_m=state['x'][ids].mean(0),centroid_velocity_m_s=vc,internal_relative_kinetic_J=float(relativeKE),remaining_live_external_connections=0,anchor_reachable=False,**shape(state['x'][ids]),internal_alive_edge_density_m3=float(alive[label]/(len(ids)*V)),internal_graph_density=float(2*alive[label]/(len(ids)*(len(ids)-1))) if len(ids)>1 else None))
  history.write(json.dumps(dict(time_s=state['time'],step=state['step'],attached_nodes=int(sizes[attached].sum()),detached_nodes=int(sizes[~attached].sum()),components=cs),default=cv)+'\n')
  return cs
 component_row(s,True)
 try:
  while prev['time']<END-1e-14:
   boundary=min([t for t in ns.tool.boundaries if t>prev['time']+1e-14]);startup=BASE/8 if prev['time']<2e-6-1e-15 else BASE/4 if prev['time']<4e-6-1e-15 else BASE/2 if prev['time']<8e-6-1e-15 else BASE
   h=min(startup,boundary-prev['time']);depth=0
   while True:
    xp=prev['x'].copy();xp[:n]+=h*(prev['velocity'][:n]+(.5*h/(RHO*V))*prev['force'][:n]);kk,dp=predict(model,ref,xp,prev['kappa'],on,off)
    if np.max(dp-prev['d'])<=.1:break
    assert h>MIN_H*(1+1e-12),'MINIMUM_EVENT_SUBSTEP_EXHAUSTED'
    h/=2;depth+=1;rejects+=1
   substeps+=bool(depth)
   ctrl=dict(nominal_amplitude_Pa=ns.tool.field.amplitude(prev['time']+h),nominal_shear_Pa=1000*float(ns.tool.field.envelope(prev['time']+h)),drive_mode='FROZEN_OPEN_LOOP')
   s=ns.step(h);p=prev;dd=s['d']-p['d'];sc=s['scalars'];nf=np.flatnonzero((s['d']==1)&(p['d']<1)&~cut)
   error('predicted_x_m',np.max(abs(xp-s['x'])));error('predicted_damage',np.max(abs(dp-s['d'])));error('predicted_kappa',np.max(abs(kk-s['kappa'])));error('irreversibility',max(np.max(-dd),np.max(p['kappa']-s['kappa'])));error('velocity_VV_m_s',np.max(abs(s['velocity'][:n]-p['velocity'][:n]-h/(2*RHO*V)*(p['force'][:n]+s['force'][:n]))))
   assert err['predicted_x_m']<2e-18 and err['velocity_VV_m_s']<1e-10 and max(err['predicted_damage'],err['predicted_kappa'],err['irreversibility'])<1e-10
   error('slip_end_velocity_m_s',np.max(abs(s['force_stage_velocity']-s['velocity'])));assert err['slip_end_velocity_m_s']<1e-10
   assert np.max(dd)<=.1+1e-10 and np.array_equal(s['x'][n:],ref['X'][n:]) and not s['velocity'][n:].any()
   assert s['native_bonds']==int(np.sum(s['d']<1)) and not s['d'][bt==3].any() and np.isfinite(s['force']).all()
   work=ns.tool.work(p,s);winc=float(work.sum());wnew=W+winc;jnew=Jc
   # Signed passive endpoint work is audited; tiny quadrature effects are retained.
   assert work[1]<=1e-17
   new_surface_work+=.5*float(np.sum((p['surface_force']*p['surface_accessible'][:,None]+s['surface_force']*s['surface_accessible'][:,None])*(s['surface_centroid']-p['surface_centroid'])))
   facet_impulse+=.5*h*float(np.linalg.norm(p['surface_force'],axis=1).sum()+np.linalg.norm(s['surface_force'],axis=1).sum())
   stage_work+=.5*h*(float(np.sum(p['external_force']*p['velocity']))+float(np.sum(s['external_force']*s['velocity'])))
   if len(nf):
    rr=np.linalg.norm(s['x'][b[nf,1]]-s['x'][b[nf,0]],axis=1);jnew+=float(.5*contact.c*np.sum(np.minimum(rr-DX,0)**2))
   res=wnew-(sc[0]-initial['scalars'][0]+s['K']-initial['K']+sc[1]+s['Uc']-initial['Uc']-jnew);newscale=max(scale,abs(wnew),abs(sc[0]-initial['scalars'][0])+s['K']+sc[1]+s['Uc']);rel=abs(res)/newscale
   assert rel<=.01,('ENERGY_RESIDUAL',rel,res,newscale)
   if s['step']%100==0 or np.any(dd):
    new=model.evaluate(s['x'],s['d']);old=model.evaluate(s['x'],p['d']) if np.any(dd) else new;uc,fc,pairs=contact.evaluate(s['x'],s['d']);error('native_force_N',np.max(abs(new['force']-s['material_force'])));error('force_port_N',np.max(abs(s['force']-s['material_force']-fc-s['external_force'])));error('native_U_J',abs(new['U']-sc[0]));error('dissipation_increment_J',abs(old['U']-new['U']-(sc[1]-p['scalars'][1])));error('pair_contact_U_J',abs(uc-s['Uc']))
    assert max(err['native_force_N'],err['force_port_N'])<1e-12 and max(err['native_U_J'],err['dissipation_increment_J'],err['pair_contact_U_J'])<1e-17
   W=wnew;Wdrive+=work[0];Wdrag+=work[1];Jc=jnew;scale=newscale
   meanforce=.5*(p['external_force']+s['external_force']);impulse+=h*meanforce.sum(0);absolute_impulse+=h*.5*(np.linalg.norm(p['external_force'],axis=1).sum()+np.linalg.norm(s['external_force'],axis=1).sum());rotation_impulse+=.5*h*(np.cross(p['x'],p['external_force']).sum(0)+np.cross(s['x'],s['external_force']).sum(0))
   surface=ns.tool.evaluate(s['x'],s['velocity'],s['d'],s['time']);geometry=ns.tool.metrics(surface,s['x']);ns.tool.commit(surface)
   error('current_surface_force_replay_N',np.max(abs(surface['force']-s['external_force'])));assert err['current_surface_force_replay_N']<1e-15
   if geometry_first is None and (geometry['max_loaded_displacement_m']>.1*H_EFF or geometry['max_normal_rotation_deg']>15):
    geometry_first=dict(time_s=s['time'],step=s['step'],**geometry);save(case,'first_layer_geometry_diagnostic.npz',s)
   for key,condition in [('first_surface_support_loss',geometry['inactive_broken_edge_faces']>0),('first_area_validity_diagnostic',geometry['max_active_area_ratio']>2 or (geometry['min_active_area_ratio'] is not None and geometry['min_active_area_ratio']<.5)),('first_degenerate_live_face',geometry['degenerate_live_faces']>0)]:
    if condition and key not in surface_events:surface_events[key]=dict(time_s=s['time'],step=s['step'],**geometry);save(case,key+'.npz',s)
   event('first_kappa_exceeds_s0',np.flatnonzero((s['kappa']>on)&~cut),s);event('first_softening',np.flatnonzero((s['d']>0)&~cut),s);event('first_failure',nf,s)
   if len(nf):cached=graph.topology(ref,s['d'])
   lab,sizes,attached=cached;free=sizes[~attached];large=int(free.max(initial=0));medium=(free>=20)&(free<=200);debris=free[free<=2].sum();cs=component_row(s,bool(len(nf)))
   if large>=2 and 'first_multinode' not in events:
    events['first_multinode']=dict(time_s=s['time'],step=s['step'],component_ids=[q['component_id'] for q in cs],evidence='keys/first_multinode.npz');save(case,'first_multinode.npz',s)
   score=selection_key(free,n)
   if score>maxscore:maxscore=score;save(case,'max_multifragment.npz',s)
   for t,key in [(RISE,'rise_end'),(HOLD,'peak_load'),(OFF,'upper_drive_off'),(COHORT,'common_cohort')]:
    if abs(s['time']-t)<1e-14:save(case,key+'.npz',s)
   if int(np.sum((s['d']==1)&~cut)) and s['failed_tension_max_N']>1e-20:raise AssertionError('FAILED_BOND_RESIDUAL_FORCE')
   if s['time']>=snap_next-1e-14:save(case,f'frame_{s["step"]:06d}.npz',s);snap_next+=25e-6
   diameter=(6*free*V/np.pi)**(1/3);dfrac=float(debris/free.sum()) if free.sum() else None;instant=scalar_metrics(free,n)['instant_target']
   if instant and 'first_instant_target' not in events:
    events['first_instant_target']=dict(time_s=s['time'],step=s['step'],evidence='keys/first_instant_target.npz');save(case,'first_instant_target.npz',s)
   row=dict(step=s['step'],time_s=s['time'],dt_s=h,phase=ns.tool.phase(s['time']),shear_peak_Pa=geometry['actual_peak_slip_shear_Pa'],net_force_x_N=float(s['external_force'][:,0].sum()),net_force_y_N=float(s['external_force'][:,1].sum()),net_force_z_N=float(s['external_force'][:,2].sum()),sum_force_norm_N=float(np.linalg.norm(s['external_force'],axis=1).sum()),net_moment_z_Nm=float(np.cross(s['x'],s['external_force'])[:,2].sum()),absolute_impulse_Ns=absolute_impulse,impulse_x_Ns=impulse[0],impulse_y_Ns=impulse[1],impulse_z_Ns=impulse[2],moment_impulse_z_Nms=rotation_impulse[2],bulk_failed=int(np.sum((s['d']==1)&~cut&(bt==1))),interface_failed=int(np.sum((s['d']==1)&~cut&(bt==2))),softened=int(np.sum((s['d']>0)&~cut)),attached_nodes=int(sizes[attached].sum()),detached_nodes=int(free.sum()),components=len(free),largest_nodes=large,largest_fraction=large/n,N20_200=int(medium.sum()),mass_fraction20_200=float(free[medium].sum()/n),N_ge34=int(np.sum(free>=34)),N_ge114=int(np.sum(free>=114)),N_ge200um=int(np.sum(diameter>=200e-6)),N_ge300um=int(np.sum(diameter>=300e-6)),debris12_nodes=int(debris),debris12_detached_fraction=dfrac,instant_target=instant,W_J=W,U_J=float(sc[0]),K_J=s['K'],D_J=float(sc[1]),Uc_J=s['Uc'],Jc_J=Jc,numerical_damping_work_J=0.,R_J=res,relative_R=rel,max_dd=float(dd.max()),max_speed_m_s=float(np.linalg.norm(s['velocity'],axis=1).max()),max_kappa=float(s['kappa'][~cut].max()))
   row.update(scalar_metrics(free,n));row.update(ctrl)
   row.update(cumulative_unique_loaded_area_m2=ns.tool.cumulative_new_area,ever_loaded_face_count=int(ns.tool.ever_loaded.sum()),initial_loaded_area_m2=ns.tool.initial_area,W_accessible_new_faces_J=new_surface_work,facet_shear_impulse_Ns=facet_impulse)
   row.update(geometry)
   active=surface['active'];free_faces=active&np.r_[~attached[lab],np.zeros(N-n,dtype=bool)][ns.tool.tri[:,0]];supported=np.unique(ns.tool.tri[active]);free_ids=np.flatnonzero(~attached[lab]);
   row.update(W_drive_J=float(Wdrive),W_drag_J=float(Wdrag),W_endpoint_power_J=stage_work,work_path_vs_endpoint_J=W-stage_work,drive_power_W=surface['drive_power'],drag_power_W=surface['drag_power'],active_free_faces=int(free_faces.sum()),free_nodes_with_retained_surface=int(np.isin(free_ids,supported).sum()),free_nodes_without_resolved_surface=int((~np.isin(free_ids,supported)).sum()),slip_iterations=int(s['slip_iterations']),Re_h_outer_max=float(RHO*H_EFF/MU*np.linalg.norm(surface['outer_velocity'][active],axis=1).max(initial=0)),Re_h_slip_max=float(RHO*H_EFF/MU*np.linalg.norm(surface['traction'][active]/BETA,axis=1).max(initial=0)))
   if writer is None:writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader()
   writer.writerow(row);rows.append(row);changed=np.flatnonzero(dd>0)
   if len(changed):changes.write(json.dumps(dict(step=s['step'],time_s=s['time'],original_bond_ids=changed+1,d=s['d'][changed]),default=cv)+'\n')
   if s['step']%250==0:stream.flush();history.flush();changes.flush();print('STEP',name,s['step'],s['time'],'tau',row['shear_peak_Pa'],'soft',row['softened'],'N20',row['N20_200'],'largest',large,'R',rel,'wall',time.perf_counter()-start,flush=True)
   prev=s
  save(case,'final.npz',s)
 except BaseException as exc:
  stop.append(repr(exc));(case/'FAILURE_TRACEBACK.txt').write_text(traceback.format_exc());save(case,'failed_unaccepted.npz',s);s=prev;save(case,'last_accepted.npz',s)
 finally:stream.close();history.close();changes.close()
 count,raw=ns.l.gather_bonds();actual=np.ctypeslib.as_array(raw,shape=(3*count,)).reshape(-1,3)[:,1:].copy();actual.sort(1);wanted=b[ns.field(3)<1]+1;native_topology=bool(np.array_equal(np.sort(actual[:,0]*N+actual[:,1]),np.sort(wanted[:,0]*N+wanted[:,1])));assert native_topology
 assert np.array_equal(ns.field(13),on) and np.array_equal(ns.field(14),off);ns.l.command(f"write_restart {case / 'final.restart'}");ns.close()
 summary=dict(status='COMPLETED' if not stop else 'NUMERICALLY_UNRESOLVED',name=name,uniform_material=(name=='uniform'),accepted_steps=len(rows),time_s=float(s['time']),stop_reasons=stop,events=events,last_row=rows[-1] if rows else None,max_independent_errors=err,max_relative_energy_residual=max((q['relative_R'] for q in rows),default=0),max_medium_components=max((q['N20_200'] for q in rows),default=0),history_max_free_nodes=max((q['largest_nodes'] for q in rows),default=0),secondary_events=secondary,predictive_substeps=substeps,predictive_rejections=rejects,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,library_sha256=sha(LIB),native_topology_matches_current_native=True,original_onsets_unchanged=True,saved_state_count=len(list((case/'keys').glob('*.npz'))),label=ns.tool.label,mechanical_status=LABEL,layer_geometry_first_diagnostic=geometry_first,surface_events=surface_events,field_parameters=ns.tool.field.theta,field_applicability=ns.tool.field.applicability(),new_crack_wetting='QUASISTATIC_FLUID_ACCESSIBILITY_PROXY',surface_force_scope='Supported current material tetra boundaries; exterior-connected dilated gap graph; current solid occlusion probes; no attachment-based force switch',work_definition='Actual accepted endpoint forces averaged and dotted with actual native drift; signed drive and drag retained separately; no residual added to D')
 np.savez_compressed(case/'accepted_surface_history.npz',ever_loaded_faces=ns.tool.ever_loaded,tri_atom_ids=ns.tool.tri+1)
 dump(case/'COMPONENT_REGISTRY.json',registry);dump(case/'FRACTURE_EVENTS.json',events);dump(case/'SUMMARY.json',summary);print('DONE',name,summary['status'],summary['wall_s'],flush=True);return summary
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--name',required=True,choices=CASES);args=p.parse_args();result=run(args.name);raise SystemExit(0 if result['status']=='COMPLETED' else 2)
