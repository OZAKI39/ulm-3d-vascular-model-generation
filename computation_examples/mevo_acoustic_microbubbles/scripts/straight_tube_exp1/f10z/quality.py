"""Evidence summaries; no change to any simulated state or historical policy."""
import json,resource,time
import numpy as np,h5py
from .common import *
def main():
 runs=read(BASE/'analysis/RUN_ANALYSIS.json');checks=read(BASE/'analysis/NUMERICAL_CHECKS.json');audit=read(BASE/'INITIAL_STATE_AUDIT.json')
 nominal=[r for r in runs.values() if r['name'].endswith('DT1') and r['requested_s']==.002]
 offs=[r for r in runs.values() if r['action']=='OFF']
 rows=[]
 for r in nominal:
  with h5py.File(BASE/'dataset'/(r['name']+'.h5')) as h:
   axes=h['principal_axis'][:] if 'principal_axis' in h else None
   # Use stored covariance, with sign-insensitive directions and a degeneracy check.
   Q=h['tensor_m2'][[0,-1]];eig,vec=np.linalg.eigh(Q);ok=np.all((eig[:,-1]-eig[:,-2])>1e-3*eig[:,-1])
   angle=float(np.degrees(np.arccos(np.clip(abs(vec[0,:,-1]@vec[1,:,-1]),0,1)))) if ok else None
   aspect=np.sqrt(eig[:,-1]/np.maximum(eig[:,0],1e-30))
  rows.append(dict(name=r['name'],action=r['action'],initial=r['initial'],COM_um=np.asarray(r['COM_displacement_m'])*1e6,delta_Rg_nm=np.diff(r['Rg_initial_final_m'])[0]*1e9,delta_NN_um=np.diff(r['nn_initial_final_m'])[0]*1e6,aspect3d_initial_final=aspect,axis_change_deg=angle,Q_relative_change=float(np.linalg.norm(Q[1]-Q[0])/np.linalg.norm(Q[0])),Rg_error_nm=(r.get('Rg_numerical_difference_m') or 0)*1e9 if r.get('refinement_admitted') else None,refined=r.get('refinement_admitted',False),components=r['components_initial_final'],contacts=r['contact_initial_final'],largest_component=r['largest_component_initial_final']))
 force=read(BASE/'a_formal_v1/FORCE_DIAGNOSTICS.json');mechanisms=[]
 for r in force['rows']:
  if r['status']!='PASS':continue
  x=np.array(r['position_m']);rad=x-x.mean(0);unit=rad/np.linalg.norm(rad,axis=1)[:,None];vel={}
  for k in ['primary_only_frozen_velocity_m_s','secondary_only_frozen_velocity_m_s','unconstrained_velocity_m_s']:
   v=np.array(r[k]);rel=v-v.mean(0);vel[k]=dict(COM_um_s=v.mean(0)*1e6,mean_inward_radial_um_s=float(np.mean(np.sum(rel*unit,axis=1)))*1e6)
  mechanisms.append(dict(initial=r['initial'],action=r['action'],velocity_decomposition=vel,net_F1_N=r['net_primary_N'],net_F2_N=r['net_secondary_N'],gradient_error=r['gradient_finite_difference_relative']))
 monitor=[]
 for p in BASE.glob('*/resources.json'):
  rows0=read(p)
  if rows0:monitor+=rows0
 campaigns=[read(p) for p in BASE.glob('*/CAMPAIGN.json')]
 resources=dict(monitor_samples=len(monitor),missing_stage_A_time_series=True,reason='A monitor thread lacked psutil; physics completed; B/C use stdlib /proc monitor. Static audit and per-step timing retained.',
   CPU_user_s=sum(c['CPU_user_s'] for c in campaigns),CPU_system_s=sum(c['CPU_system_s'] for c in campaigns),
   GPU_utilization_percent_max=max((r.get('gpu_percent',0) for r in monitor),default=None),
   GPU_utilization_percent_mean=float(np.mean([r['gpu_percent'] for r in monitor if 'gpu_percent' in r])) if monitor else None,
   GPU_memory_MiB_max=max((r.get('gpu_MiB',0) for r in monitor),default=None),
   process_RSS_peak_MiB=max(r['rss_peak_MiB'] for r in runs.values()),
   measured_step_s=sum(r.get('measured_compute_s',{}).get('full_step',0) for r in runs.values()),
   measured_QP_s=sum(r.get('measured_compute_s',{}).get('QP_solver',0) for r in runs.values()),
   measured_resistance_s=sum(r.get('measured_compute_s',{}).get('resistance',0) for r in runs.values()))
 maximum={k:max(r[k] for r in runs.values()) for k in ['max_overlap_m','max_contact_total_m','max_wall_correction_m','max_force_residual_N','max_actual_complementarity_m']}
 comparative={metric:max(c[metric] for g in checks.values() for c in g.get('comparisons',[])) for metric in ['COM_m','Rg_m','fixed_id_m','Q_relative','aspect_relative','matched_RMS_over_R']}
 independent={(r['initial'],r['action']) for r in nominal}
 out=dict(nominal_2ms=rows,independent_initial_action_groups=len(independent),numerical_replicates_not_independent=True,OFF_max_displacement_m=max(np.linalg.norm(r['COM_displacement_m']) for r in offs),OFF_max_Rg_change_m=max(abs(np.diff(r['Rg_initial_final_m'])[0]) for r in offs),numerical_maximum=maximum,comparison_maximum=comparative,refinement_groups=len(checks),refinement_admitted=sum(bool(c['admitted']) for c in checks.values()),force_diagnostics_pass=sum(r['status']=='PASS' for r in force['rows']),force_diagnostics_total=len(force['rows']),maximum_gradient_FD_relative=max(r['gradient_finite_difference_relative'] for r in force['rows'] if r['status']=='PASS'),mechanisms=mechanisms,resources=resources)
 write(BASE/'analysis/QUALITY_SUMMARY.json',out)
 print(json.dumps({k:v for k,v in out.items() if k not in ['nominal_2ms','mechanisms']},indent=2))
if __name__=='__main__':main()

