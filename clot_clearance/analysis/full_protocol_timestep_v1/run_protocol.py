"""Exactly one native physical path; independent audits never advance positions."""
from native import *
import resource,csv,traceback

def main():
 freeze=json.loads((O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json').read_text());g=freeze['gates'];budget=freeze['budget'];dt=freeze['dt_s']
 assert sha(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json')==(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.sha256').read_text().strip()
 for n,h in freeze['driver_sha256'].items():assert sha(A/n)==h,(n,'modified after freeze')
 assert sha(LIB)==freeze['library_sha256']
 with (O/'FORMAL_CASE_STARTED.json').open('x') as f:json.dump(dict(pid=os.getpid(),dt_s=dt,restart=False,case_count=1),f)
 start=time.perf_counter();cpu0=time.process_time();ns=Native('full_dt025',dt=dt);ref=ns.ref;b=ref['bonds'];N=ns.N;NB=ns.NB;model=SofteningLPS(ref);contact=Contact(ref)
 ids=np.load(O/'provenance/341_boundary_ids.npz');left=ids['left_ids'];right=ids['right_ids'];free=ids['free_ids'];fi=free-1;bc=np.r_[left,right]-1
 ns.dynamics(left,right,free);s=ns.state();initial=s
 assert s['step']==0 and s['time']==0 and np.array_equal(s['x'],ref['X'])
 assert not np.any(s['d']) and not np.any(s['kappa']) and not np.any(s['velocity']) and s['scalars'][1]==0
 assert ns.l.extract_atom('mass')[1]==RHO*V
 assert np.array_equal(ns.field(4),ref['m']) and np.array_equal(ns.field(9),np.arange(1,NB+1))
 # One same-state run0 in this same instance, no trial deformation/restart matrix.
 ns.l.command('run 0');rr=ns.state();assert all(np.array_equal(s[k],rr[k]) for k in ['x','velocity','force','material_force','kappa','d']) and np.array_equal(s['scalars'][:6],rr['scalars'][:6])
 dump(O/'INITIAL_STATE_AND_LOADING_IDENTITY.json',dict(status='PASS',step=0,time_s=0,unloaded=True,zero_velocity=True,zero_kappa=True,zero_damage=True,initial_U_J=s['scalars'][0],initial_K_J=s['K'],initial_D_J=s['scalars'][1],initial_force_max_N=float(abs(s['force']).max()),original_X_exact=True,original_m_exact=True,bond_ids_exact=True,mass_kg=RHO*V,left_count=len(left),right_count=len(right),free_count=len(free),restart_used=False,read0_same_state_readonly=True,prior_trial_and_restart_evidence='Reused SHA-verified certificates; no new broad tests',reference_sha256=freeze['reference_sha256'],boundary=freeze['boundary']))
 rows=[];states={k:[] for k in ['x','velocity','force','material_force','scalars','counters']};hist_steps=[];hist_ids=[];hist_k=[];hist_d=[];events=[];keys=[];stop=[];contact_events=[];seen_failed=set();firstsoft=None;firstfail=None;sumdef=np.zeros(4);W=Jc=cwork=0.;scale=1e-15;max_dd_record=0.;maxerr={k:0. for k in ['force_N','U_J','D_increment_J','kappa','d','native_KE_J','contact_force_N','contact_U_J','force_sum_N','torque_Nm','fixed_m','loading_x_m','loading_v_m_s','decomposition_J','failed_tension_N','d_rollback','kappa_rollback']};minD=0.;maxpenetration=0.;peakrss=0;next_contact_first=True
 cols=['step','time_s','dt_s','W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','relative_R','energy_scale_J','max_speed_m_s','max_free_speed_m_s','max_free_disp_m','step_disp_dx','relative_endpoint_step_dx','max_dd','max_kappa','max_d','new_soft','new_fail','total_failed','live_bonds','contact_pairs','Px_Ns','Py_Ns','Pz_Ns','momentum_residual_Ns','x_equation_m','v_equation_m_s','Kfree_J','contact_work_J','signed_topology_J','commits','token','force_evaluations','VV_defect_J','spatial_defect_J','switch_work_J','contact_defect_J','remainder_J','max_contact_penetration_dx','relative_contact_step_dx']
 def keep(state,reason):
  step=state['step'];path=O/'keyframes'/f'step_{step:06d}.npz'
  if not path.exists():np.savez_compressed(path,**state,bond_ids=np.arange(1,NB+1),bonds=b,reference_X=ref['X'])
  keys.append(dict(step=step,time_s=state['time'],reason=reason,path=str(path.relative_to(O))))
 def record(state,row,writer,stream):
  rows.append(row)
  for k in states:states[k].append(state[k].copy())
  writer.writerow(row);stream.flush()
 def update_err(k,val):maxerr[k]=max(maxerr[k],float(val))
 def endpoint_ids_ok():
  nb,raw=ns.l.gather_bonds();bb=np.ctypeslib.as_array(raw,shape=(3*nb,)).reshape(-1,3)[:,1:].copy();bb.sort(axis=1)
  expected=b[s['d']<1]+1
  return nb==len(expected) and np.array_equal(np.sort(bb[:,0]*N+bb[:,1]),np.sort(expected[:,0]*N+expected[:,1]))
 keep(s,'initial');stream=(O/'step_scalars.csv').open('w');writer=csv.writer(stream);writer.writerow(cols)
 initialrow=[0,0,dt,0,s['scalars'][0],s['K'],0,s['Uc'],0,0,0,scale,0,0,0,0,0,0,0,0,0,0,0,NB,0,0,0,0,0,0,0,0,0,0,s['scalars'][4],s['scalars'][5],s['scalars'][6],0,0,0,0,0,0,0]
 assert len(initialrow)==len(cols);record(s,initialrow,writer,stream);prev=s;last_independent=model.evaluate(s['x'],s['d']);prevfc=np.zeros_like(s['x']);last_newfail=[]
 try:
  for step in range(1,freeze['steps']+1):
   ns.l.command('run 1 pre no post no');s=ns.state();dx=s['x']-prev['x'];dv=s['velocity']-prev['velocity'];sc=s['scalars'];dd=s['d']-prev['d'];newfail=np.flatnonzero((s['d']==1)&(prev['d']<1));newsoft=np.flatnonzero((s['d']>0)&(prev['d']==0));changed=np.flatnonzero(dd!=0)
   if len(changed):
    hist_steps.append(np.full(len(changed),step,np.int32));hist_ids.append((changed+1).astype(np.int32));hist_k.append(s['kappa'][changed]);hist_d.append(s['d'][changed])
   if firstsoft is None and len(newsoft):firstsoft=dict(step=step,time_s=s['time'],bond_ids=newsoft+1);keep(s,'first_softening')
   if firstfail is None and len(newfail):firstfail=dict(step=step,time_s=s['time'],bond_ids=newfail+1);keep(s,'first_complete_failure')
   stretch=np.linalg.norm(s['x'][b[:,1]]-s['x'][b[:,0]],axis=1)/ref['r0']-1
   if len(newsoft) or len(newfail):events.append(dict(step=step,time_s=s['time'],new_soft_ids=newsoft+1,new_failed_ids=newfail+1,failure_stretches=stretch[newfail],failure_kappa=s['kappa'][newfail]))
   if seen_failed.intersection((newfail+1).tolist()):stop.append('DUPLICATE_ORIGINAL_FAILURE')
   seen_failed.update((newfail+1).tolist())
   if len(newfail) and np.any(stretch[newfail]<SC-g['history_abs']):stop.append('FAILURE_BELOW_FROZEN_THRESHOLD')
   update_err('d_rollback',np.max(-dd));update_err('kappa_rollback',np.max(prev['kappa']-s['kappa']));minD=min(minD,float(sc[2]))
   # All independent calculations are read-only, using this step's native x,v.
   new=model.evaluate(s['x'],s['d']);old=model.evaluate(s['x'],prev['d']);kk,ddref=model.update(s['x'],prev['kappa'],s_complete=SC)
   update_err('force_N',abs(new['force']-s['material_force']).max());update_err('U_J',abs(new['U']-sc[0]));update_err('D_increment_J',abs(old['U']-new['U']-sc[2]));update_err('kappa',abs(kk-s['kappa']).max());update_err('d',abs(ddref-s['d']).max());update_err('native_KE_J',abs(.5*RHO*V*np.sum(s['velocity']**2)-s['K']))
   pairs=contact.tree(s['x']).query_pairs(DX,output_type='ndarray');pairs=pairs[~np.isin(pairs[:,0]*N+pairs[:,1],contact.codes[s['d']<1])]
   penetration=rel_contact=0.
   if len(pairs):
    dist=np.linalg.norm(s['x'][pairs[:,1]]-s['x'][pairs[:,0]],axis=1);penetration=float(np.max(1-dist/DX));rel_contact=float(np.max(np.linalg.norm(dx[pairs[:,1]]-dx[pairs[:,0]],axis=1)/DX));maxpenetration=max(maxpenetration,penetration)
    if np.min(dist/DX)<.10 or rel_contact>g['max_relative_endpoint_step_dx']:stop.append('CONTACT_GEOMETRIC_SAFETY')
    if next_contact_first:
     keep(s,'first_active_contact');contact_events.append(dict(step=step,time_s=s['time'],pairs=pairs+1,distances_m=dist,original_failed_ids=[int(contact.codes.tolist().index(int(i*N+j))+1) for i,j in pairs if i*N+j in contact.codes[s['d']==1]],penetration_dx=penetration));next_contact_first=False
   Uc,fc,nc=contact.evaluate(s['x'],s['d']);update_err('contact_force_N',abs(fc-(s['force']-s['material_force'])).max());update_err('contact_U_J',abs(Uc-s['Uc']))
   W+=float(-.5*np.sum((prev['force'][bc]+s['force'][bc])*dx[bc])+.5*RHO*V*np.sum(s['velocity'][bc]**2-prev['velocity'][bc]**2))
   dJ=0.
   if len(newfail):
    dist=np.linalg.norm(s['x'][b[newfail,1]]-s['x'][b[newfail,0]],axis=1);dJ=float(.5*contact.c*np.sum(np.minimum(dist-DX,0)**2))
   Jc+=dJ;cwork+=float(.5*np.sum((prevfc+fc)*dx));D=float(sc[1]-initial['scalars'][1]);U=float(sc[0]);Know=s['K'];Rnow=W-(U-initial['scalars'][0]+Know-initial['K']+D+s['Uc']-initial['Uc']-Jc);scale=max(scale,abs(W),abs(U-initial['scalars'][0])+abs(Know-initial['K'])+D+abs(s['Uc']-initial['Uc']));relR=abs(Rnow)/scale
   Kfree=.5*RHO*V*np.sum(s['velocity'][fi]**2);dKfree=.5*RHO*V*np.sum(s['velocity'][fi]**2-prev['velocity'][fi]**2)
   defects=np.array([dKfree-.5*np.sum((prev['force'][fi]+s['force'][fi])*dx[fi]),old['U']-last_independent['U']+.5*np.sum((prev['material_force']+old['force'])*dx),.5*np.sum((s['material_force']-old['force'])*dx),s['Uc']-prev['Uc']-dJ+.5*np.sum((prevfc+fc)*dx)])
   sumdef+=defects;remainder=Rnow+sumdef.sum();update_err('decomposition_J',abs(remainder))
   impulse_bc=RHO*V*np.sum(dv[bc],axis=0)-dt/2*np.sum(prev['force'][bc]+s['force'][bc],axis=0);mom=float(np.linalg.norm(RHO*V*np.sum(dv,axis=0)-impulse_bc));P=RHO*V*np.sum(s['velocity'],axis=0)
   xeq=float(abs(dx[fi]-dt*prev['velocity'][fi]-dt*dt/(2*RHO*V)*prev['force'][fi]).max());veq=float(abs(dv[fi]-dt/(2*RHO*V)*(prev['force'][fi]+s['force'][fi])).max())
   speed=float(np.linalg.norm(s['velocity'],axis=1).max());freespeed=float(np.linalg.norm(s['velocity'][fi],axis=1).max());freedisp=float(np.linalg.norm(s['x'][fi]-ref['X'][fi],axis=1).max());disp=float(np.linalg.norm(dx,axis=1).max()/DX);rdisp=float(np.linalg.norm(dx[b[:,1]]-dx[b[:,0]],axis=1).max()/DX);maxdd=float(np.max(dd));failed=s['d']==1
   tdir=ns.field(6);update_err('failed_tension_N',np.max(abs(tdir[np.r_[failed,failed]]),initial=0));update_err('force_sum_N',np.linalg.norm(s['force'].sum(axis=0)));update_err('torque_Nm',np.linalg.norm(np.cross(s['x'],s['force']).sum(axis=0)))
   t=s['time'];ql=min(t/.0003,1.);qu=0. if t<.00035 else min((t-.00035)/.0003,1.);u=30e-6*(np.cos(np.pi*qu)-np.cos(np.pi*ql));v=60e-6*np.pi/(2*.0003)*(np.sin(np.pi*ql)-np.sin(np.pi*qu));xb=ref['X'][bc].copy();vb=np.zeros_like(xb);xb[:len(left),0]-=u;xb[len(left):,0]+=u;vb[:len(left),0]=-v;vb[len(left):,0]=v
   update_err('loading_x_m',abs(s['x'][bc]-xb).max());update_err('loading_v_m_s',abs(s['velocity'][bc]-vb).max())
   row=[step,t,dt,W,U,Know,D,s['Uc'],Jc,Rnow,relR,scale,speed,freespeed,freedisp,disp,rdisp,maxdd,float(s['kappa'].max()),float(s['d'].max()),len(newsoft),len(newfail),int(failed.sum()),s['native_bonds'],nc,*P,mom,xeq,veq,Kfree,cwork,-D,sc[4],sc[5],sc[6],*sumdef,remainder,penetration,rel_contact]
   record(s,row,writer,stream)
   if not np.all(np.isfinite(np.r_[row,s['x'].ravel(),s['velocity'].ravel(),s['force'].ravel(),s['kappa'],s['d']])):stop.append('NONFINITE')
   if sc[4]!=step or sc[5]!=step or s['step']!=step or abs(t-step*dt)>1e-14:stop.append('STEP_COMMIT_IDENTITY')
   if s['native_bonds']!=NB-failed.sum() or not np.array_equal(ns.field(9),np.arange(1,NB+1)):stop.append('ORIGINAL_TOPOLOGY_IDENTITY')
   if maxerr['d_rollback']>g['history_abs'] or maxerr['kappa_rollback']>g['history_abs'] or sc[2]<g['no_negative_D_below_J']:stop.append('DAMAGE_IRREVERSIBILITY')
   if max(maxerr[k] for k in ['kappa','d'])>g['history_abs']:stop.append('HISTORY_REFERENCE_MISMATCH')
   if max(maxerr[k] for k in ['force_N','contact_force_N','failed_tension_N'])>g['force_abs_N']:stop.append('FORCE_REFERENCE_OR_FAILED_TENSION')
   if max(maxerr[k] for k in ['U_J','D_increment_J','native_KE_J','contact_U_J','decomposition_J'])>g['energy_abs_J']:stop.append('ENERGY_REFERENCE_OR_UNEXPLAINED_REMAINDER')
   if xeq>g['integrator_equation_abs_m'] or veq>g['velocity_abs_m_s'] or mom>g['momentum_absolute_Ns']:stop.append('INTEGRATOR_MOMENTUM')
   if maxerr['loading_x_m']>g['integrator_equation_abs_m'] or maxerr['loading_v_m_s']>g['velocity_abs_m_s']:stop.append('PHYSICAL_LOADING_IDENTITY')
   if relR>g['energy_relative_stop']:stop.append('ENERGY_RELATIVE_STOP')
   if speed>g['max_speed_m_s']:stop.append('MAX_SPEED')
   if disp>g['max_step_displacement_dx']:stop.append('NODE_STEP_DISPLACEMENT')
   if rdisp>g['max_relative_endpoint_step_dx']:stop.append('RELATIVE_ENDPOINT_STEP_DISPLACEMENT')
   if maxdd>g['max_damage_increment']:stop.append('MAX_DAMAGE_INCREMENT')
   rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024;peakrss=max(peakrss,rss)
   livebytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file())
   memoryarrays=sum(z.nbytes for arr in states.values() for z in arr)+sum(z.nbytes for arrays in [hist_steps,hist_ids,hist_k,hist_d] for z in arrays)
   if time.perf_counter()-start>budget['main_wall_s'] or rss>budget['RSS_stop_GiB']*1024**3 or livebytes+memoryarrays>budget['payload_limit_MiB']*1024**2:stop.append('RESOURCE_BUDGET')
   if maxdd>max_dd_record and maxdd>=.08:keep(s,'record_damage_increment_at_least_0.08')
   max_dd_record=max(max_dd_record,maxdd)
   if step in [1200,1400,2600,2800]:keep(s,'physical_phase_transition')
   if stop:
    keep(prev,'immediately_before_safety_stop');keep(s,'FIRST_SAFETY_STOP_NOT_ADMITTED');break
   prev=s;prevfc=fc;last_independent=new
   if step%200==0:print('STEP',step,'time_us',t*1e6,'failed',int(failed.sum()),'max_dd',maxdd,'energy_rel',relR,'wall_s',time.perf_counter()-start,flush=True)
 except Exception as exc:
  stop.append('EXCEPTION: '+repr(exc));(O/'EXCEPTION.txt').write_text(traceback.format_exc());keep(prev,'before_exception');keep(ns.state(),'exception_state')
 finally:
  elapsed=time.perf_counter()-start;cpu=time.process_time()-cpu0;stream.close();keep(s,'final_executed_state');ns.l.command(f'write_restart {O/"final_executed.restart"}')
  topology_match=endpoint_ids_ok();update_err('fixed_m',abs(ns.field(4)-ref['m']).max());maps=sorted(set(z.split()[-1] for z in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in z));ns.close()
 arr=np.array(rows);np.savez_compressed(O/'trajectory.npz',series=arr,columns=cols,**{k:np.asarray(z) for k,z in states.items()},reference_X=ref['X'],bonds=b,bond_ids=np.arange(1,NB+1),left_ids=left,right_ids=right,free_ids=free)
 np.savez_compressed(O/'damage_history_sparse.npz',step=np.concatenate(hist_steps) if hist_steps else np.array([],np.int32),bond_id=np.concatenate(hist_ids) if hist_ids else np.array([],np.int32),kappa=np.concatenate(hist_k) if hist_k else np.array([],np.float64),d=np.concatenate(hist_d) if hist_d else np.array([],np.float64))
 violations=np.flatnonzero(s['d']-prev['d']>.1) if stop else np.array([],int)
 detailed=[]
 for k in violations:
  nodes=b[k];detailed.append(dict(bond_id=int(k+1),node_ids=nodes+1,kappa_before=prev['kappa'][k],kappa_after=s['kappa'][k],d_before=prev['d'][k],d_after=s['d'][k],delta_d=s['d'][k]-prev['d'][k],stretch_before=np.linalg.norm(prev['x'][nodes[1]]-prev['x'][nodes[0]])/ref['r0'][k]-1,stretch_after=np.linalg.norm(s['x'][nodes[1]]-s['x'][nodes[0]])/ref['r0'][k]-1,velocity_before_m_s=prev['velocity'][nodes],velocity_after_m_s=s['velocity'][nodes],force_before_N=prev['force'][nodes],force_after_N=s['force'][nodes],material_force_before_N=prev['material_force'][nodes],material_force_after_N=s['material_force'][nodes]))
 dump(O/'NUMERICAL_SAFETY_STOP_RECORD.json',dict(status='FAILED_ADMISSION' if stop else 'NO_SAFETY_STOP',reasons=stop,first_stop_step=s['step'] if stop else None,first_stop_time_s=s['time'] if stop else None,last_admitted_step=prev['step'] if stop else s['step'],violating_state_committed_but_not_admitted=bool(stop),no_step_after_stop=True,violating_bonds=detailed,ledger=dict(zip(cols,arr[-1])),restart='final_executed.restart',freeze_sha256=sha(O/'FULL_PROTOCOL_TIMESTEP_DESIGN_FREEZE.json')))
 dump(O/'FIRST_SOFTENING_AND_FAILURE_EVENTS.json',dict(first_softening=firstsoft,first_complete_failure=firstfail,events=events,unique_failed_ids=sorted(seen_failed)))
 dump(O/'KEYFRAME_INDEX.json',keys)
 if contact_events:dump(O/'DYNAMIC_CONTACT_CLOSURE_AUDIT.json',dict(status='OBSERVED_CHECK_AUDITS',events=contact_events,max_penetration_dx=maxpenetration,max_force_error_N=maxerr['contact_force_N'],max_U_error_J=maxerr['contact_U_J'],history='All dynamic contact pair counts and penetration recorded per actual step. Active pairs exclude every alive original bond.',near_coincidence_guard_dx=.1))
 audit=dict(status='FAILED_ADMISSION' if stop else 'COMPLETED',executed_steps=s['step'],last_admitted_step=prev['step'] if stop else s['step'],final_time_s=s['time'],planned_steps=2800,dt_s=dt,stop_reasons=stop,max_errors=maxerr,min_D_increment_J=minD,final_ledger=dict(zip(cols,arr[-1])),maxima={c:float(arr[:,i].max()) for i,c in enumerate(cols)},native_endpoint_topology_matches=topology_match,actual_loaded_libraries={p:sha(p) for p in maps},independent_model_checked_every_step=True,full_states_saved=len(keys),actual_dynamic_cases=1,restart_continuations=0,contact_status='OBSERVED' if contact_events else 'UNOBSERVED')
 dump(O/'FULL_PROTOCOL_NATIVE_TIMESTEP_AUDIT.json',audit)
 dump(O/'RUN_RESOURCES.json',dict(main_wall_s=elapsed,main_CPU_s=cpu,peak_RSS_bytes=peakrss,wall_with_final_compression_s=time.perf_counter()-start,MPI_ranks=1,threads=1,actual_dynamic_cases=1,actual_physical_steps=s['step'],recompile_count=0,GPU_used=False,output_bytes_now=sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),budget=budget))
 print('FORMAL_DONE',json.dumps({k:audit[k] for k in ['status','executed_steps','final_time_s','stop_reasons','contact_status']}),flush=True);print('MAX_ERRORS',json.dumps(maxerr),flush=True)

if __name__=='__main__':main()
