"""One declared M5C run; passive diagnostics and preregistered accepted-state stops."""
from pathlib import Path
import sys,csv,json,math,numpy as np
from pd_clot.lammps_fragmentation.macro_lab_v1 import run
from metrics import Topology
from design import save,js
class Observer:
 def __init__(self,out,g,mask,sign,plan):
  self.out=out;self.g=g;self.X=g['X'];self.b=g['bonds'];self.mask=mask;self.sign=sign;self.plan=plan;self.rules=js(out.parent/'M5C_DESIGN_FREEZE.json');self.top=Topology(g,mask);self.previous=None;self.first=False;self.first_coherent=False;self.central_saved=False;self.separation_saved=False;self.peak=-np.inf;self.peak_state=None;self.peak_step=0;self.low=(mask[:,2]&(sign[:,2]<0));self.high=(mask[:,2]&(sign[:,2]>0));self.max_dx=0.;self.window_max_dx=0.;self.max_position_recurrence_error=0.;self.max_velocity_recurrence_error=0.;self.max_native_v_error=0.;self.max_boundary_position_error=0.;self.max_boundary_velocity_error=0.;self.max_native_x_error=0.;self.hold=0;self.actual_steps=0;self.stop_reason='STATIC_PREFLIGHT' if plan['steps']==0 else None;self.checkpoints=[];self.trace=[];self.event_markers={};self.commands=[]
  self.fs=(out/'boundary_dynamics_audit.csv').open('w');self.writer=csv.writer(self.fs);self.writer.writerow(['step','time_s','top_Rz_N','bottom_Rz_N','generalized_Q_N','actual_opening_m','max_speed_m_s','max_step_dx_over_dx','boundary_position_error_m','boundary_velocity_error_m_s','free_position_recurrence_error_m','free_velocity_recurrence_error_m_s','notch_live_util','grip_live_util','surface_live_util','interior_live_util','notch_pre_event_util','grip_pre_event_util','surface_pre_event_util','interior_pre_event_util','U_bulk_J','U_contact_J','K_free_J','W_boundary_J','D_numeric_J','signed_topology_J','residual_J','max_contact_force_N','contact_force_decomposition_valid'])
 def state(self,name,step,x,v,alive):np.savez_compressed(self.out/(name+'.npz'),step=step,time_s=step*self.plan['dt_s'],x_m=x,v_full_m_s=v,alive=alive)
 def record(self,step,x,v,fm,fd,alive,changed,F):
  assert self.previous is None and step==0 or self.previous is not None and step==self.previous[0]+1,'Duplicated/skipped accepted callback; cannot treat as valid trajectory'
  dt=self.plan['dt_s'];t=step*dt;d=.00068*.5*(1-math.cos(math.pi*t/.024));vel=.00068*.5*math.pi/.024*math.sin(math.pi*t/.024);pe=float(np.max(abs((x-self.X-d*self.sign)[self.mask])));ve=float(np.max(abs((v-vel*self.sign)[self.mask])));assert pe<1e-12 and ve<1e-10
  dxmax=xe=v_err=0.;force=fm+fd
  if self.previous is not None:
   st,xp,vp,ap,fp=self.previous;dxmax=float(np.linalg.norm(x-xp,axis=1).max()/50e-6);expected_x=xp+dt*vp+.5*dt*dt*fp/(1000*(50e-6)**3);expected_v=vp+.5*dt*(fp+force)/(1000*(50e-6)**3);xe=float(np.max(abs((x-expected_x)[~self.mask])));v_err=float(np.max(abs((v-expected_v)[~self.mask])));assert xe<1e-12 and v_err<1e-9,'Free-DOF velocity Verlet mismatch'
  self.max_dx=max(self.max_dx,dxmax);self.window_max_dx=max(self.window_max_dx,dxmax);self.max_position_recurrence_error=max(self.max_position_recurrence_error,xe);self.max_velocity_recurrence_error=max(self.max_velocity_recurrence_error,v_err);self.max_boundary_position_error=max(self.max_boundary_position_error,pe);self.max_boundary_velocity_error=max(self.max_boundary_velocity_error,ve)
  R=-float(fm[self.high,2].sum());bottom=-float(fm[self.low,2].sum());Q=-float(np.sum(fm*self.mask*self.sign))
  if R>self.peak:self.peak=R;self.peak_step=step;self.peak_state=(step,x.copy(),v.copy(),alive.copy())
  self.write_this=step%20==0 or changed
  if self.write_this:
   util=(np.linalg.norm(x[self.b[:,1]]-x[self.b[:,0]],axis=1)/self.top.r-1)/.15;now=alive.astype(bool);prior=self.previous[3].astype(bool) if self.previous is not None else now;maximum=lambda eligible:[float(util[(self.top.reg==j)&eligible].max()) if np.any((self.top.reg==j)&eligible) else 0 for j in range(4)];self.pending=[step,t,R,bottom,Q,float((x[self.high,2]-self.X[self.high,2]).mean()-(x[self.low,2]-self.X[self.low,2]).mean()),float(np.linalg.norm(v,axis=1).max()),dxmax,pe,ve,xe,v_err,*maximum(now),*maximum(prior)]
  if changed and not self.first:
   self.state('pre_first_break',*self.previous[:4]);self.first=True;self.event_markers['first_break']=dict(step=step,time_s=t)
  if changed and not self.first_coherent:
   c=self.top.cluster(alive)
   if c['largest_voxels']>=10:self.first_coherent=True;self.event_markers['first_coherent_cluster']=dict(step=step,time_s=t,**c);self.state('first_coherent_cluster',step,x,v,alive)
  if step%1000==0:self.state(f'frame_{step:06d}',step,x,v,alive)
  self.previous=(step,x.copy(),v.copy(),alive.copy(),force.copy());self.actual_steps=step
  if step%1000==0:self.fs.flush()
 def record_energy(self,step,U,Uc,K,W,Dn,Db,fm,fp,changed):
  if self.write_this:
   valid=not changed or Uc<1e-25;cf=0. if Uc<1e-25 else float(np.linalg.norm(fm-fp,axis=1).max()) if not changed else float('nan');self.writer.writerow([*self.pending,U,Uc,K,W,Dn,Db,W-U-Uc-K-Dn-Db,cf,int(valid)])
 def advance(self,L,steps):
  while self.actual_steps<steps and self.stop_reason is None:
   remaining=min(100,steps-self.actual_steps);cmd=f'run {remaining} pre '+('yes' if self.actual_steps==0 else 'no')+' post no';self.commands.append(cmd);L.command(cmd)
   if (self.out/'CALLBACK_ERROR.txt').exists():raise RuntimeError('Callback diagnostic failed; stop without successful execution')
   step,x,v,alive,force=self.previous;assert step==int(L.extract_global('ntimestep'));nx=L.numpy.extract_atom('x')[:len(x)];nv=L.numpy.extract_atom('v')[:len(v)];self.max_native_x_error=max(self.max_native_x_error,float(np.max(abs(nx-x))));self.max_native_v_error=max(self.max_native_v_error,float(np.max(abs((nv-v)[~self.mask]))));assert self.max_native_x_error<1e-12 and self.max_native_v_error<1e-9
   q=self.top.summary(alive);q.update(step=step,time_s=step*self.plan['dt_s'],top_Rz_N=self.pending[2] if int(self.pending[0])==step else -float(L.numpy.extract_atom('f')[self.high,2].sum()),running_positive_peak_N=self.peak,window_max_step_dx_over_dx=self.window_max_dx)
   sizes=q['component_sizes'];separated=len(sizes)>1 and sizes[1]>=.2*len(x) and sum(sizes[:2])>=.8*len(x)
   if q['central_advancing_strips_2dx']>=2 and not self.central_saved:self.central_saved=True;self.state('first_central_advance',step,x,v,alive);self.event_markers['first_central_advance']=q.copy()
   if separated and not self.separation_saved:self.separation_saved=True;self.state('first_large_separation',step,x,v,alive);self.event_markers['first_large_separation']=q.copy()
   lowload=separated and self.peak>50*.0006**2 and q['top_Rz_N']<.1*self.peak;self.hold=self.hold+100 if lowload else 0
   if self.window_max_dx>.2:self.stop_reason='NUMERICAL_STEP_DISPLACEMENT_LIMIT'
   elif q['broken']>=1000 and q['grip_fraction']>=.5:self.stop_reason='SEVERE_GRIP_DAMAGE_CONTAMINATION'
   elif q['broken']>=1000 and q['target_band_fraction']<.5 and q['off_plane_fraction']>.6:self.stop_reason='SINGLE_TARGET_BAND_INTERPRETATION_LOST'
   elif self.hold>=2000:self.stop_reason='POST_SEPARATION_LOW_LOAD_OBSERVATION_COMPLETE'
   elif step==steps:self.stop_reason='MAXIMUM_24MS'
   q['stop_reason']=self.stop_reason;self.checkpoints.append(q);self.window_max_dx=0
   if step%1000==0 or self.stop_reason:save(self.out/'PROGRESS.json',dict(step=step,time_s=q['time_s'],broken=q['broken'],stop_reason=self.stop_reason))
  save(self.out/'STOP_DECISION.json',dict(reason=self.stop_reason,actual_steps=self.actual_steps,time_s=self.actual_steps*self.plan['dt_s'],checking_interval_steps=100,design_freeze=str(self.out.parent/'M5C_DESIGN_FREEZE.json'),no_future_data=True,last_checkpoint=self.checkpoints[-1]));save(self.out/'CHECKPOINT_TOPOLOGY.json',self.checkpoints)
 def close(self):
  self.fs.close()
  if self.peak_state is not None:self.state('peak_tensile',*self.peak_state)
  save(self.out/'OBSERVER_CHECKS.json',dict(actual_steps=self.actual_steps,max_boundary_position_error_m=self.max_boundary_position_error,max_boundary_velocity_error_m_s=self.max_boundary_velocity_error,max_free_position_recurrence_error_m=self.max_position_recurrence_error,max_free_velocity_recurrence_error_m_s=self.max_velocity_recurrence_error,max_native_end_block_x_error_m=self.max_native_x_error,max_native_end_block_free_v_error_m_s=self.max_native_v_error,max_step_dx_over_dx=self.max_dx,events=self.event_markers,stop_reason=self.stop_reason,contact_force_note='0 when native conservative contact energy zero; otherwise native total minus independent bulk only when no topology change. During nonzero-contact break steps decomposition unavailable, not silently called contact.'))
  (self.out/'run_commands.txt').write_text('\n'.join(self.commands)+'\n')
if __name__=='__main__':run(Path(sys.argv[1]).resolve(),'M5C',static='--static' in sys.argv,observer_factory=Observer)
