from native import *
import resource,shutil,traceback

def run_case(name,root,left,right,free,steps,axis=0,amplitude=60e-6,ramp=.0003,hold=.00005,dt=.5e-6,wall_limit=300,restart=None,diagnostic=False):
 start=time.perf_counter();cpu0=time.process_time();ns=Native(name,root=root,restart=restart,dt=dt);ref=ns.ref;model=SofteningLPS(ref);contact=Contact(ref);N=ns.N;NB=ns.NB;bc=np.r_[left,right]-1;fi=np.asarray(free)-1;b=ref['bonds'];g=json.loads((O/'DYNAMIC_DAMAGE_COMMIT_DESIGN_FREEZE.json').read_text())['gates']
 if not restart and (Path(root)/'initial_damage.npz').exists():
  initial=np.load(Path(root)/'initial_damage.npz');ns.static(kappa=initial['kappa'],d=initial['d'])
 ns.dynamics(left,right,free,axis,amplitude,ramp,hold);prev=ns.state();start_state=prev;initial_step=prev['step'];U0=prev['scalars'][0];K0=prev['K'];Uc0=prev['Uc'];D0=prev['scalars'][1];W=0.;Jc=0.;contact_work=0.;max_energy_scale=1e-15;rows=[];stored={key:[] for key in ['x','velocity','force','material_force','kappa','d','scalars','counters']};events=[];checks=[];safety=[];stop=[];peak_rss=0;checkpoint=[];first_soft=None;first_failure=None;max_rollback=0.;max_force=0.;max_U=0.;max_momentum=0.;max_xeq=0.;max_veq=0.;max_ct=0.;max_ctU=0.;max_failed_t=0.;n_contact_steps=0;retry_trigger=None
 cols=['step','time_s','W_J','U_J','K_J','D_J','Uc_J','contact_topology_J','residual_J','relative_residual','vmax_m_s','step_disp_dx','relative_endpoint_step_dx','max_dd','new_failures','total_failed','commits','force_evaluations','native_bonds','momentum_residual_Ns','x_equation_m','v_equation_m_s','contact_work_J','signed_topology_J','Kfree_J','contact_pairs']
 def record_state(s,row):
  rows.append(row)
  # Full FP64 accepted histories for P3; sparse plus event states for larger cases.
  keep=N==341 or s['step']%10==0 or row[14]>0 or not rows
  if keep:
   for key in stored:stored[key].append(s[key].copy())
   saved_steps.append(s['step'])
 saved_steps=[]
 def readonly(s,trial=False):
  before=s['scalars'][:6].copy();ns.l.command('run 0');z=ns.state();ok=np.array_equal(before,z['scalars'][:6]) and np.array_equal(s['d'],z['d']) and np.array_equal(s['kappa'],z['kappa']);ferr=float(abs(s['force']-z['force']).max());rc=ns.lib.lpsvar_action(ns.l.lmp,2,s['step']+1,None,None);assert rc==-7
  if trial:
   ids=ns.gather('id');xx=s['x'].copy();xx[fi[0],axis]+=DX*.25;ns.l.numpy.extract_atom('x')[:N]=xx;ns.l.command('run 0');zz=ns.state();ok=ok and np.array_equal(s['d'],zz['d']) and np.array_equal(s['kappa'],zz['kappa']) and np.array_equal(before[1:],zz['scalars'][1:6]);ns.l.numpy.extract_atom('x')[:N]=s['x'];ns.l.command('run 0');z=ns.state();ferr=max(ferr,float(abs(z['force']-s['force']).max()))
  safety.append(dict(step=s['step'],run0_readonly=bool(ok),trial_coordinates_only=trial,explicit_commit_rejection=rc,restored_force_abs_N=ferr));assert ok and ferr<2e-17
  return z
 prev=readonly(prev,trial=bool(np.any(prev["d"]>0)))
 record_state(prev,[prev['step'],prev['time'],0,U0,K0,0,Uc0,0,0,0,float(np.linalg.norm(prev['velocity'],axis=1).max()),0,0,0,0,int(np.sum(prev['d']==1)),int(prev['scalars'][4]),int(prev['scalars'][6]),prev['native_bonds'],0,0,0,0,0,float(.5*RHO*V*np.sum(prev['velocity'][fi]**2)),0])
 try:
  for it in range(1,steps+1):
   ns.l.command('run 1 pre no post no');s=ns.state();dx=s['x']-prev['x'];dv=s['velocity']-prev['velocity'];sc=s['scalars'];newfail=np.flatnonzero((s['d']==1)&(prev['d']<1));newsoft=np.flatnonzero((s['d']>0)&(prev['d']==0));dd=float(np.max(s['d']-prev['d']));rollback=max(float(np.max(prev['d']-s['d'])),float(np.max(prev['kappa']-s['kappa'])));max_rollback=max(max_rollback,rollback)
   if len(newsoft) and first_soft is None:first_soft=dict(step=s['step'],time_s=s['time'],bond_ids=(newsoft+1).tolist())
   if len(newfail) and first_failure is None:first_failure=dict(step=s['step'],time_s=s['time'],bond_ids=(newfail+1).tolist())
   if len(newfail):events.append(dict(step=s['step'],time_s=s['time'],bond_ids=(newfail+1).tolist(),current_stretches=(np.linalg.norm(s['x'][b[newfail,1]]-s['x'][b[newfail,0]],axis=1)/ref['r0'][newfail]-1).tolist(),kappa=s['kappa'][newfail].tolist(),signed_U_jump_J=-float(sc[2])))
   W+=float(-.5*np.sum((prev['force'][bc]+s['force'][bc])*dx[bc])+.5*RHO*V*np.sum(s['velocity'][bc]**2-prev['velocity'][bc]**2))
   Uc,fc,nc=contact.evaluate(s['x'],s['d']);n_contact_steps+=int(nc>0);max_ct=max(max_ct,float(abs(fc-(s['force']-s['material_force'])).max()));max_ctU=max(max_ctU,abs(Uc-s['Uc']))
   if len(newfail):
    dist=np.linalg.norm(s['x'][b[newfail,1]]-s['x'][b[newfail,0]],axis=1);Jc+=float(.5*contact.c*np.sum(np.minimum(dist-DX,0)**2))
   contact_work+=float(.5*np.sum(((prev['force']-prev['material_force'])+fc)*dx))
   Dnow=float(sc[1]-D0);res=W-(sc[0]-U0+s['K']-K0+Dnow+s['Uc']-Uc0-Jc);max_energy_scale=max(max_energy_scale,abs(W),abs(sc[0]-U0)+abs(s['K']-K0)+Dnow+abs(s['Uc']-Uc0));rel=abs(res)/max_energy_scale
   impulse_bc=RHO*V*np.sum(dv[bc],axis=0)-dt/2*np.sum(prev['force'][bc]+s['force'][bc],axis=0);mom=float(np.linalg.norm(RHO*V*np.sum(dv,axis=0)-impulse_bc));max_momentum=max(max_momentum,mom)
   xeq=float(abs(dx[fi]-dt*prev['velocity'][fi]-dt*dt/(2*RHO*V)*prev['force'][fi]).max());veq=float(abs(dv[fi]-dt/(2*RHO*V)*(prev['force'][fi]+s['force'][fi])).max());max_xeq=max(max_xeq,xeq);max_veq=max(max_veq,veq)
   vmax=float(np.linalg.norm(s['velocity'],axis=1).max());disp=float(np.linalg.norm(dx,axis=1).max()/DX);rdisp=float(np.linalg.norm(dx[b[:,1]]-dx[b[:,0]],axis=1).max()/DX)
   tdir=ns.field(6);failed=s['d']==1;max_failed_t=max(max_failed_t,float(np.max(abs(tdir[np.r_[failed,failed]]),initial=0)))
   if it==1 or it%100==0 or len(newfail) and first_failure['step']==s['step'] or len(newsoft) and first_soft['step']==s['step']:
    q=model.evaluate(s['x'],s['d']);ferr=float(abs(q['force']-s['material_force']).max());uerr=abs(q['U']-sc[0]);max_force=max(max_force,ferr);max_U=max(max_U,uerr);checks.append(dict(step=s['step'],time_s=s['time'],force_abs_N=ferr,energy_abs_J=uerr));
   row=[s['step'],s['time'],W,sc[0],s['K'],Dnow,s['Uc'],Jc,res,rel,vmax,disp,rdisp,dd,len(newfail),int(failed.sum()),int(sc[4]),int(sc[6]),s['native_bonds'],mom,xeq,veq,contact_work,-Dnow,float(.5*RHO*V*np.sum(s['velocity'][fi]**2)),nc];record_state(s,row)
   if not np.all(np.isfinite(np.r_[row,s['x'].ravel(),s['velocity'].ravel(),s['kappa'],s['d']])):stop.append('NONFINITE')
   if sc[4]!=start_state['scalars'][4]+it or sc[5]!=s['step'] or s['step']!=initial_step+it:stop.append('STEP_COMMIT_IDENTITY')
   if s['native_bonds']!=NB-int(failed.sum()) or not np.array_equal(ns.field(9),np.arange(1,NB+1)):stop.append('ORIGINAL_TOPOLOGY_IDENTITY')
   if rollback>g['history_abs'] or sc[2]<g['no_negative_D_below_J']:stop.append('DAMAGE_IRREVERSIBILITY')
   if max_force>g['force_abs_N'] or max_U>g['energy_abs_J'] or max_ct>g['force_abs_N'] or max_ctU>g['energy_abs_J']:stop.append('FORCE_ENERGY_REFERENCE')
   if xeq>g['integrator_equation_abs_m'] or veq>g['velocity_abs_m_s'] or mom>g['momentum_absolute_Ns']:stop.append('INTEGRATOR_MOMENTUM')
   if rel>g['energy_relative_stop']:stop.append('ENERGY_RELATIVE_STOP')
   if vmax>g['max_speed_m_s'] or disp>g['max_step_displacement_dx'] or rdisp>g['max_relative_endpoint_step_dx'] or dd>g['max_damage_increment']:stop.append('TIMESTEP_SAFETY')
   if retry_trigger is None and (rel>g['energy_relative_risk'] or dd>g['max_damage_increment']):retry_trigger=dict(step=s['step'],time_s=s['time'],relative_energy=rel,damage_increment=dd,checkpoint=checkpoint[-1] if checkpoint else None)
   rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024;peak_rss=max(peak_rss,rss)
   if time.perf_counter()-start>wall_limit or rss>4*1024**3:stop.append('RESOURCE_BUDGET')
   # P4 target-area guard: stop on dominant full failures at the loaded end surfaces.
   if N>341 and len(newfail):
    mid=(ref['X'][b[failed,0]]+ref['X'][b[failed,1]])/2;end=np.abs(mid[:,axis])>np.max(abs(ref['X'][:,axis]))-HORIZON
    if end.sum()>max(10,.5*len(mid)):stop.append('LOADING_END_DOMINANT_FAILURE')
   if stop:
    np.savez_compressed(O/(name+'_STOP_STATE.npz'),**s);ns.l.command(f'write_restart {O/(name+"_STOP.restart")}');break
   if not diagnostic and (it%40==0 or first_soft and first_soft['step']==s['step']):
    cp=O/(name+f'_checkpoint_{it}.restart');ns.l.command('write_restart '+str(cp));checkpoint.append(dict(step=s['step'],time_s=s['time'],path=str(cp))); 
    # Keep only two recent checkpoints, unless selected for the registered dt/2 window.
    if len(checkpoint)>2:
     old=checkpoint[-3]
     if not retry_trigger or old!=retry_trigger['checkpoint']:Path(old['path']).unlink(missing_ok=True)
   if it==1 or first_soft and first_soft['step']==s['step'] or first_failure and first_failure['step']==s['step']:s=readonly(s,trial=bool(first_soft and first_soft['step']==s['step']))
   prev=s
   if it%200==0:print(name,'step',s['step'],'t',s['time'],'failed',int(failed.sum()),'Erel',rel,'wall',time.perf_counter()-start,flush=True)
 except Exception as exc:
  stop.append('EXCEPTION: '+repr(exc));(O/(name+'_exception.txt')).write_text(traceback.format_exc());np.savez_compressed(O/(name+'_EXCEPTION_STATE.npz'),**ns.state())
 finally:
  ns.l.command(f'write_restart {O/(name+"_final.restart")}');final=ns.state();ns.close()
 arr=np.asarray(rows);np.savez_compressed(O/(name+'_trajectory.npz'),series=arr,columns=cols,saved_steps=saved_steps,**{k:np.asarray(v) for k,v in stored.items()},reference_X=ref['X'],bonds=b,bond_ids=np.arange(1,NB+1),left_ids=left,right_ids=right,free_ids=free)
 summary=dict(name=name,status='STOPPED' if stop else 'COMPLETED',stop_reasons=stop,accepted_steps=len(rows)-1,dt_s=dt,initial_step=initial_step,final_time_s=float(arr[-1,1]),simulated_time_s=float(arr[-1,1]-arr[0,1]),N=N,NB=NB,free_nodes=len(free),free_motion_max_m=float(np.max(np.linalg.norm(final['x'][fi]-start_state['x'][fi],axis=1))),first_softening=first_soft,first_complete_failure=first_failure,final_failed_bonds=int(arr[-1,15]),force_abs_N=max_force,energy_abs_J=max_U,contact_force_abs_N=max_ct,contact_energy_abs_J=max_ctU,contact_steps=n_contact_steps,history_max_decrease=max_rollback,failed_endpoint_tension_N=max_failed_t,momentum_residual_max_Ns=max_momentum,position_equation_max_m=max_xeq,velocity_equation_max_m_s=max_veq,energy_relative_max=float(arr[:,9].max()),final_ledger=dict(zip(cols,arr[-1])),timesteps=dict(vmax_m_s=float(arr[:,10].max()),max_step_disp_dx=float(arr[:,11].max()),max_relative_endpoint_dx=float(arr[:,12].max()),max_damage_increment=float(arr[:,13].max())),sensitivity_risk=retry_trigger,wall_s=time.perf_counter()-start,CPU_s=time.process_time()-cpu0,peak_RSS_bytes=peak_rss,trajectory_bytes=(O/(name+'_trajectory.npz')).stat().st_size)
 dump(O/(name+'_AUDIT.json'),summary);dump(O/(name+'_EVENTS.json'),dict(first_softening=first_soft,first_complete_failure=first_failure,complete_failure_events=events));dump(O/(name+'_REFERENCE_CHECKS.json'),checks);dump(O/(name+'_RUN0_TRIAL.json'),safety)
 print('CASE_DONE',name,json.dumps(summary),flush=True);return summary

if __name__=='__main__':
 assert json.loads((O/'P2_FIXED_STATE_REGRESSION.json').read_text())['status']=='PASS'
 assert not (O/'DYNAMIC_341_NODE_AUDIT.json').exists(),'Formal main case already executed; do not replay'
 ids=np.load(O/'provenance/341_boundary_ids.npz');result=run_case('dynamic341',O/'provenance',ids['left_ids'],ids['right_ids'],ids['free_ids'],1400);dump(O/'DYNAMIC_341_NODE_AUDIT.json',result)
