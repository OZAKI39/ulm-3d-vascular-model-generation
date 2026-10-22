"""Analyze only saved accepted states; preserve failures and conditional AI gates."""
import time,json
import numpy as np,h5py
from scipy.spatial.distance import pdist
from mevo_mb.fast_sim.f10g_contact import resistance
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.control.scheduled_pzt import Schedule
from .common import *
from .state import export,series

def run_key(initial,action,duration=.002,variant='DT1'):return f'{initial}__{action}__{round(duration*1e6)}us__{variant}'
def main():
 start=time.perf_counter();out=BASE/'analysis';out.mkdir(exist_ok=True);dataset=BASE/'dataset';dataset.mkdir(exist_ok=True);runs={};cfg=read(CONFIG);gates={}
 for p in sorted(BASE.glob('*/runs/*/RESULT.json')):
  r=read(p);name=p.parent.name;raw=p.parent/'trajectory.h5';assert sha(raw)==r['trajectory_sha256'];m=export(p.parent,dataset/(name+'.h5'));t=m['time_s'];x=m['centered_positions_m']+m['centroid_m'][:,None]
  with h5py.File(raw) as h:
   x=h['positions_m'][:];fields={k:h['force_samples'][k][:] for k in h['force_samples']};steps=h['step_metrics'][:];extra=h['extra_metrics'][:];step_t=h['step_time_duration_segment'][:]
  pressure_minmax=[float(extra[:,3].min()),float(extra[:,4].max())] if len(extra) else None
  balance=[];drag=[];drag_times=[]
  if fields:
   for i,tau in enumerate(fields['time_s']):
    j=int(np.searchsorted(t,tau-1e-13))
    if j>=len(t)-1:continue
    zeta=6*np.pi*F2Physics(radius_m=R,frequency_Hz=cfg["frequency_Hz"]).viscosity_Pa_s*R;A=resistance(x[j],R)['A'];force_drag=-zeta*(A@m['velocity_m_s'][j].ravel()).reshape(32,3);drag.append(force_drag);drag_times.append(tau);balance.append(float(np.linalg.norm(force_drag+fields['total_force_N'][i]+fields['contact_force_N'][i],axis=1).max()))
  # A fixed-ID network ledger; new proximity edges are not physical fusion.
  changed=np.diff(m['edges'].astype(int),axis=0);new=(~m['edges'][0])&m['edges'][-1];edge_persistence=[]
  for idx in np.flatnonzero(new):
   last_off=np.flatnonzero(~m['edges'][:,idx]);since=float(t[last_off[-1]+1]) if len(last_off) else 0.;edge_persistence.append(float(t[-1]-since))
  with h5py.File(dataset/(name+'.h5'),'a') as h:
   action=Schedule.from_record(read(p.parent/'INPUT.json')['schedule']);segments=action.segments;sel=np.minimum(np.searchsorted([s.start_s for s in segments],t,side='right')-1,len(segments)-1)
   for key,value in dict(segment_id=np.array([segments[i].segment_id for i in sel]),amplitude_Pa_equivalent=np.array([segments[i].amplitude_Pa_equivalent for i in sel]),phase_rad=np.array([segments[i].phase_rad for i in sel]),enabled=np.array([segments[i].enabled for i in sel])).items():h.create_dataset(key,data=value,compression='gzip')
   fs=h.create_group('force_samples')
   for key,value in fields.items():fs.create_dataset(key,data=value,compression='gzip')
   if drag:
    fs.create_dataset('derived_drag_N',data=drag,compression='gzip');fs.create_dataset('derived_drag_time_s',data=drag_times)
   h.create_dataset('step_metrics',data=steps,compression='gzip');h.create_dataset('step_time_duration_segment',data=step_t,compression='gzip')
   h.attrs['execution_commit']=r.get('execution_commit',read(p.parent/'INPUT.json')['identity']['execution_commit'])
  endpoint_Q=m['xy_tensor_m2'][-1]-m['xy_tensor_m2'][0];y=np.r_[m['centroid_m'][-1,:2]-m['centroid_m'][0,:2],m['xy_rg_m'][-1]-m['xy_rg_m'][0],endpoint_Q[0,0],endpoint_Q[0,1],endpoint_Q[1,1]]
  r.update(initial=p.parent.name.split('__')[0],action=p.parent.name.split('__')[1],path=str(p.parent.relative_to(ROOT)),samples=len(t),endpoint_state_delta=y,COM_displacement_m=m['centroid_m'][-1]-m['centroid_m'][0],Rg_initial_final_m=m['rg_m'][[0,-1]],xy_Rg_initial_final_m=m['xy_rg_m'][[0,-1]],nn_initial_final_m=m['nn_m'][[0,-1]],pair_quantiles_initial_final_m=m['pair_quantiles_m'][[0,-1]],components_initial_final=m['components'][[0,-1]],largest_component_initial_final=m['largest_component'][[0,-1]],edge_initial_final=m['edge_count'][[0,-1]],new_edges_final=int(new.sum()),new_edges_persistent_0p2ms=int(np.sum(np.array(edge_persistence)>=.0002-1e-13)),new_edges_persistent_2ms=int(np.sum(np.array(edge_persistence)>=.002-1e-13)),edge_events=int(np.sum(changed!=0)),contact_initial_final=m['contact_count'][[0,-1]],relative_Rg_change=float(m['rg_m'][-1]/m['rg_m'][0]-1),xy_aspect_initial_final=m['xy_aspect'][[0,-1]],shape_tensor_initial_final_m2=m['tensor_m2'][[0,-1]],axis_identifiable_initial_final=m['axis_identifiable'][[0,-1]],centered_max_displacement_m=float(np.linalg.norm(m['centered_positions_m'][-1]-m['centered_positions_m'][0],axis=1).max()),rigid_residual_final_m=float(m['kabsch_from_initial_m'][-1]),relative_radial_mean_m_s=float(np.nanmean(m['radial_mean_m_s'])) if len(t)>1 else None,inward_interval_fraction=float(np.mean(m['radial_mean_m_s']<0)) if len(t)>1 else None,actual_pressure_Pa=pressure_minmax,force_balance_sample_max_N=max(balance,default=None),dataset_sha256=sha(dataset/(name+'.h5')))
  r['force_summary']={key:dict(mean_norm_N=float(np.linalg.norm(fields[key],axis=-1).mean()),max_norm_N=float(np.linalg.norm(fields[key],axis=-1).max()),mean_vector_N=fields[key].mean((0,1))) for key in ['primary_force_N','secondary_force_N','total_force_N','contact_force_N'] if key in fields}
  if len(steps):r['measured_compute_s']={name:float(steps[:,col].sum()) for name,col in [('full_step',6),('QP_solver',7),('resistance',8)]}
  runs[name]=r
 for p in BASE.glob('checks/*.json'):gates[p.stem]=read(p)
 # Every response gets a matching measured time-step/cache error, not a global noise guess.
 for name,r in runs.items():
  if not name.endswith('DT1'):continue
  gate=gates.get(name[:-5]);comps=gate['comparisons'] if gate and gate.get('admitted') else []
  r['refinement_admitted']=bool(comps);noise=max([s['Rg_m'] for s in comps],default=float('inf'))
  r['Rg_numerical_difference_m']=noise if np.isfinite(noise) else None
  dr=float(r['Rg_initial_final_m'][1]-r['Rg_initial_final_m'][0]);r['relative_convergence_resolved']=bool(comps and dr< -cfg['signal_noise_multiplier']*noise and r['nn_initial_final_m'][1]<r['nn_initial_final_m'][0])
  r['shape_response_resolved']=bool(comps and abs(dr)>cfg['signal_noise_multiplier']*noise)
  # Compare endpoint networks and event-time sensitivity separately from continuous shape.
  graph_agreement=False;edge_time_error=None;feature_noise=[]
  if comps:
   with h5py.File(dataset/(name+'.h5')) as h:
    e0=h['edges'][:];t0=h['time_s'][:];Q=h['xy_tensor_m2'][:];za=np.c_[h['centroid_m'][:,:2],h['xy_rg_m'][:],Q[:,0,0],Q[:,0,1],Q[:,1,1]]
   graph_agreement=True;errors=[]
   for variant in ['DT05','FINE']:
    other=name[:-3]+variant
    with h5py.File(dataset/(other+'.h5')) as h:
     ee=h['edges'][:];tt=h['time_s'][:];Q=h['xy_tensor_m2'][:];zb=np.c_[h['centroid_m'][:,:2],h['xy_rg_m'][:],Q[:,0,0],Q[:,0,1],Q[:,1,1]]
    ix=np.searchsorted(tt,t0-1e-13);assert np.max(abs(tt[ix]-t0))<1e-13;feature_noise.append(np.max(abs(za-zb[ix]),axis=0))
    graph_agreement &= bool(np.array_equal(e0[0],ee[0]) and np.array_equal(e0[-1],ee[-1]))
    for edge in range(e0.shape[1]):
     ii=np.flatnonzero(np.diff(e0[:,edge].astype(int)))+1;jj=np.flatnonzero(np.diff(ee[:,edge].astype(int)))+1
     if len(ii)!=len(jj):graph_agreement=False
     elif len(ii):errors.append(float(np.max(abs(t0[ii]-tt[jj]))))
   edge_time_error=max(errors,default=0.)
  r['state_feature_numerical_difference']=np.max(feature_noise,axis=0).tolist() if feature_noise else None
  r['graph_endpoint_event_count_agrees']=graph_agreement;r['edge_event_max_time_difference_s']=edge_time_error
  r['formation_geometric_exploratory']=bool(r['relative_convergence_resolved'] and dr/r['Rg_initial_final_m'][0]<=-.01 and r['components_initial_final'][1]<r['components_initial_final'][0] and r['new_edges_persistent_0p2ms']>0 and graph_agreement and edge_time_error<=policy()['instant_edge_time_tolerance_s'])
 write(out/'RUN_ANALYSIS.json',runs);write(out/'NUMERICAL_CHECKS.json',gates)
 # Unit-normalized finite differences, SI derivatives retained with their units.
 center=run_key('medium_train','PAIR');derivatives=[]
 for low,high,width,unit in [('AMP-5','AMP+5',.10*1036.2016581041792,'Pa-equivalent'),('PHASE-3','PHASE+3',6.,'degree'),('PHASE-6','PHASE+6',12.,'degree')]:
  aa=runs.get(run_key('medium_train',low));bb=runs.get(run_key('medium_train',high))
  if aa and bb and aa.get('refinement_admitted') and bb.get('refinement_admitted'):
   delta=np.asarray(bb['endpoint_state_delta'])-np.asarray(aa['endpoint_state_delta']);derivatives.append(dict(parameter=unit,low=low,high=high,width=width,derivative_SI_per_unit=delta/width,normalized_response=delta/np.array([1e-7,1e-7,1e-7,1e-12,1e-12,1e-12]),noise_Rg_m=sum(v['state_feature_numerical_difference'][2] for v in [aa,bb]),shape_difference_m=float(delta[2])))
 mat=np.array([v['normalized_response'] for v in derivatives[:2]]).T if len(derivatives)>=2 else np.empty((6,0));sv=np.linalg.svd(mat,compute_uv=False);rank=int(np.sum(sv>max(sv[0]*1e-6,1e-8))) if len(sv) else 0
 write(out/'ACTION_RESPONSE.json',dict(columns=derivatives,units=['m','m','m','m2','m2','m2'],state_order=['cx','cy','Rg_xy','Qxx','Qxy','Qyy'],singular_values=sv,local_rank=rank,condition_number=float(sv[0]/sv[-1]) if len(sv) and sv[-1]>0 else None,parameter_scope='Only X-pair amplitude and relative phase are refined locally; Y one-hot is separately proxy-specific. Rank is local, not nonlinear/global control.',phase_scale_stability= 'Compare ±3 and ±6 degree derivatives; no extrapolation'))
 # Entire-trajectory groups only. These deterministic baselines are not a world model.
 train=[r for r in runs.values() if r['initial'] in cfg['primary_states'] and r['action'] in ['OFF','XM','XP','PAIR'] and r['requested_s']==.002 and r['name'].endswith('DT1') and r['status']=='PASS']
 hold=[r for r in runs.values() if r['initial'] in cfg['heldout_states'] and r['requested_s']==.002 and r['name'].endswith('DT1') and r['status']=='PASS']
 baselines=[]
 for r in hold:
  with h5py.File(dataset/(r['name']+'.h5')) as h:
   ts=h['time_s'][:];z=np.c_[h['centroid_m'][:,:2],h['xy_rg_m'][:]];k=int(np.searchsorted(ts,.0005))
  actual=z[-1]-z[k];same=[v for v in train if v['action']==r['action']];mu=np.mean([np.array(v['endpoint_state_delta'])[:3] for v in same],axis=0)*((ts[-1]-ts[k])/.002) if same else np.zeros(3)
  predictions={'state_hold':np.zeros(3),'constant_velocity_first_0p5ms':(z[k]-z[0])/ts[k]*(ts[-1]-ts[k]),'training_action_mean':mu}
  baselines.append(dict(name=r['name'],split=read(BASE/'INITIAL_STATE_AUDIT.json')[r['initial']]['split'],observed_prefix_s=float(ts[k]),prediction_horizon_s=float(ts[-1]-ts[k]),errors_m={key:(value-actual).tolist() for key,value in predictions.items()},scope='Same simulated horizon only; heldout initial configuration, not time-sample random split. Actions OFF/PAIR were seen during training.'))
 action_rows=[runs[run_key('medium_train',a)] for a in ['XM','XP','PAIR','YM','YP'] if run_key('medium_train',a) in runs and runs[run_key('medium_train',a)].get('refinement_admitted')]
 action_noise=np.max([r['state_feature_numerical_difference'] for r in action_rows],axis=0) if action_rows else np.full(6,np.inf)
 action_deltas=np.array([r['endpoint_state_delta'] for r in action_rows]).T if action_rows else np.empty((6,0))
 resolved=np.where(abs(action_deltas)>5*action_noise[:,None],action_deltas,0)/np.array([1e-7,1e-7,1e-7,1e-12,1e-12,1e-12])[:,None]
 discrete_sv=np.linalg.svd(resolved,compute_uv=False);discrete_rank=int(np.sum(discrete_sv>max(discrete_sv[0]*1e-6,1e-8))) if len(discrete_sv) else 0
 write(out/'DISCRETE_ACTION_SPAN.json',dict(actions=[r['action'] for r in action_rows],feature_noise=action_noise,normalized_resolved_deltas=resolved,singular_values=discrete_sv,rank=discrete_rank,interpretation='Finite action contrasts above measured component-wise dt/cache differences. Not a local Jacobian and not proof of nonlinear control.'))
 heldout_action_and_initial_available=False
 ready=bool(discrete_rank>=3 and heldout_action_and_initial_available and len(hold)>=4)
 write(out/'WORLD_MODEL_READINESS.json',dict(status='WORLD_MODEL_DATA_READY' if ready else 'WORLD_MODEL_DATA_NOT_READY',joint_position_shape_ready=ready,discrete_action_span_rank=discrete_rank,local_action_rank=rank,independent_training_trajectories=len(train),heldout_initial_trajectories=len(hold),heldout_action_validation='NOT_RUN',world_model_trained=False,MPPI='MPPI_CONTROL_NOT_RUN',reason='No independent held-out action model validation and no demonstrated independent XY+shape authority; simple baselines only.',baselines=baselines,training_ready=False))
 campaigns=[read(p) for p in BASE.glob('*/CAMPAIGN.json')];write(out/'COSTS.json',dict(remote_campaign_wall_s=sum(v['wall_s'] for v in campaigns),saved_integration_s=sum(v['integration_s'] for v in runs.values()),analysis_wall_s=time.perf_counter()-start,completed_trajectories=sum(v['status']=='PASS' for v in runs.values()),failed_trajectories=sum(v['status']!='PASS' for v in runs.values()),longest_actual_s=max([v['completed_s'] for v in runs.values()],default=0),campaigns=[dict(stage=v['stage'],wall_s=v['wall_s'],status=v['status']) for v in campaigns]))
 print('Analyzed',len(runs),'trajectories; local response rank',rank,'; joint world model ready',ready)
if __name__=='__main__':main()
