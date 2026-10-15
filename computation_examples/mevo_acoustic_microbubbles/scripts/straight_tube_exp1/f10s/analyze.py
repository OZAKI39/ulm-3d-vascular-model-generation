"""Actual full-window metrics with original F10P–R quality thresholds."""
import time
import numpy as np
import h5py
from dataclasses import replace
from mevo_mb.control.scheduled_pzt import Schedule
from mevo_mb.fast_sim.nearfield import pair_reference
from scripts.straight_tube_exp1.f10o.analyze import load,metrics,topology_record
from scripts.straight_tube_exp1.f10p.analyze import comparison
from .admission import layered
from .common import *

def main():
 start=time.perf_counter();ident=identity();out=BASE/'analysis';out.mkdir(exist_ok=True);p=policy();cfg=read(CONFIG);c=context()
 runs={};paths={};comps={};gates={};costs=[]
 for cp in sorted(BASE.glob('*/CAMPAIGN.json')):
  costs.append(read(cp))
  for rp in sorted((cp.parent/'runs').glob('*/RESULT.json')):
   name=rp.parent.name;r=read(rp);paths[name]=rp.parent;inp=read(rp.parent/'INPUT.json');t,x=load(rp.parent)
   assert sha(rp.parent/'trajectory.h5')==r['trajectory_sha256']
   np.testing.assert_array_equal(x[0],c['x0'])
   m=metrics(t,x);top,labs,ets,eids=topology_record(t,x)
   top['initial_network_persists_entire_window']=bool(abs(t[-1]-r['requested_s'])<1e-13 and top['max_components']==1)
   top['requested_window_s']=r['requested_s']
   err=np.linalg.norm(m['centroid_m']-cfg['target_m'],axis=1)
   np.savez_compressed(out/(name+'.npz'),**m,labels=labs,edge_times_s=ets,edge_ids=eids,target_error_m=err)
   with h5py.File(rp.parent/'trajectory.h5') as h:
    rows=h['step_metrics'][:];extra=h['extra_metrics'][:];fs={k:h['force_samples'][k][:] for k in h['force_samples']}
    assert h.attrs['flow']=='STRICT_ZERO'
   assert np.max(np.linalg.norm(x[:,:,1:],axis=-1))+p['radius_m']<.0008
   assert np.max(abs(x[:,:,0]))+p['radius_m']<.01625
   i,j=np.triu_indices(32,1);actual_overlap=float(np.maximum(0,2*p['radius_m']-np.linalg.norm(x[:,i]-x[:,j],axis=-1)).max())
   r.update(topology=top,rg_initial_final_m=m['rg_m'][[0,-1]],aspect_initial_final=m['aspect'][[0,-1]],
    Q_initial_final_m2=m['tensor_m2'][[0,-1]],xy_initial_final=m['directional_xy'][[0,-1]],
    target_error_initial_final_m=err[[0,-1]],closest_time_s=t[int(np.argmin(err))],closest_error_m=float(err.min()),
    endpoint_overshoot_x=bool(m['centroid_m'][-1,0]>cfg['target_m'][0]),
    target_error_increased_any_step=bool(np.any(np.diff(err)>1e-16)),
    actual_pressure_min_max_Pa=[float(extra[:,3].min()),float(extra[:,4].max())],actual_overlap_from_all_positions_m=actual_overlap,
    solver_s=float(rows[:,7].sum()),resistance_s=float(rows[:,8].sum()),
    force_sample_time_range_s=fs['time_s'][[0,-1]],
    force_samples_summary={key:dict(mean_norm_N=float(np.linalg.norm(fs[key],axis=-1).mean()),max_norm_N=float(np.linalg.norm(fs[key],axis=-1).max()),mean_vector_N=fs[key].mean((0,1))) for key in ['primary_force_N','secondary_force_N','contact_force_N']},
    mean_COM_velocity_m_s=(m['centroid_m'][-1]-m['centroid_m'][0])/t[-1],
    centered_fixed_ID_shape_change_m=float(np.linalg.norm(m['centered_positions_m'][-1]-m['centered_positions_m'][0],axis=1).max()),
    Rg_relative_change=float(m['rg_m'][-1]/m['rg_m'][0]-1),
    Q_relative_change=float(np.linalg.norm(m['tensor_m2'][-1]-m['tensor_m2'][0])/np.linalg.norm(m['tensor_m2'][0])))
   if name.endswith('DT1'):
    e=Schedule.from_record(inp['schedule']).segments[0].event();u=c['basis'].vector(e['amplitude_Pa'],e['phase_rad'],e['enabled'])
    snapshots=[];phys=c['physics'];chi=c['response'].query(np.array([complex(c['response'].knots[0])]))[0,0]/c['response'].knots[0]
    for k in [0,-1]:
     pp,gg=c['basis'].query_basis(x[k]);pressure=u@pp;gradient=np.einsum('i,inj->nj',u,gg);v=c['response'].query(pressure)
     f1=-.5*np.real(v[:,0,None]*gradient.conj());f2=pair_reference(x[k],v,pressure,phys,chi);far=pair_reference(x[k],v,pressure,replace(phys,force_model='f1_smooth_near_field'),chi)
     snapshots.append(dict(time_s=t[k],positions_m=x[k],pressure_Pa=pressure,primary_N=f1,secondary_N=f2,F2_increment_N=f2-far))
    write(out/(name+'_FROZEN_FORCE.json'),snapshots)
    r['endpoint_F2_increment_mean_max_N']=[float(np.linalg.norm(snapshots[-1]['F2_increment_N'],axis=-1).mean()),float(np.linalg.norm(snapshots[-1]['F2_increment_N'],axis=-1).max())]
   runs[name]=r
 for name,r in runs.items():
  if not name.endswith('_DT1'):continue
  stem=name[:-4];tops=[r['topology']];cs=[]
  for suffix in ['DT05','FINE']:
   other=stem+'_'+suffix
   if other not in paths:cs.append(None);continue
   s,a=comparison(paths[name],paths[other],'Unchanged full physics; '+suffix)
   ta,_=load(paths[name]);tb,_=load(paths[other])
   s['complete_window']=bool(abs(ta[-1]-r['requested_s'])<1e-13 and abs(tb[-1]-r['requested_s'])<1e-13 and s['common_frames']==len(ta))
   s['requested_window_s']=r['requested_s'];s['fixed_id_pass']=bool(s['complete_window'] and s['fixed_id_m']<=p['fixed_id_tolerance_m'])
   comps[name+'__'+suffix]=s;cs.append(s);tops.append(runs[other]['topology'])
   np.savez_compressed(out/('compare_'+name+'__'+suffix+'.npz'),**a)
  g=layered(r,*cs,tops,p,r['requested_s'])
  g.update(target_definition='X target applies to A only; B is direction verification, no target registered',TARGET_REACHED='NOT_RUN',STRICT_ENDPOINT_REACHED='NOT_RUN')
  if stem=='X_PAIR_3MS':
   g['TARGET_REACHED']='PASS' if g['simulation_group_position_eligible'] and r['target_error_initial_final_m'][-1]<=cfg['target_tolerance_m'] else 'FAIL'
   g['STRICT_ENDPOINT_REACHED']='PASS' if g['simulation_group_position_eligible'] and r['target_error_initial_final_m'][-1]<=cfg['strict_endpoint_tolerance_m'] else 'FAIL'
  g['shape_change_exceeds_Rg_dt_cache_difference']=bool(abs(np.diff(r['rg_initial_final_m'])[0])>max((v['Rg_m'] for v in cs if v),default=np.inf))
  g['shape_cost_Rg_relative']=r['Rg_relative_change'];g['max_components']=r['topology']['max_components'];gates[name]=g
 write(out/'RUN_ANALYSIS.json',runs);write(out/'COMPARISONS.json',comps)
 write(out/'FINAL_GATE.json',dict(groups=gates,TRAP_FORCE_STABILITY_CHECKED='PASS',PDMS_FIELD_ADMITTED='FAIL',PDMS_FIELD='NOT_RUN',EXPERIMENTALLY_CALIBRATED='NOT_VERIFIED',TRAINING_READY=False,labels=cfg['labels'],training_ready=False))
 write(out/'COSTS.json',dict(remote_process_s=sum(v['wall_s'] for v in costs),local_analysis_s=time.perf_counter()-start,
   integration_s=sum(v['integration_s'] for v in runs.values()),CPU_solver_s=sum(v['solver_s'] for v in runs.values()),CPU_resistance_s=sum(v['resistance_s'] for v in runs.values()),GPU_kernel_only_s=None,
   identity=ident,max_real_particle_time_s=max(v['completed_s'] for v in runs.values()),note='Remote campaign wall time includes checks, startup and writes; analysis and diagnostics counted separately. GPU kernel-only time not measured.'))
 print({k:{a:b for a,b in v.items() if a.startswith(('BASIC','GROUP','PARTICLE','TARGET','STRICT'))} for k,v in gates.items()})
if __name__=='__main__':main()
