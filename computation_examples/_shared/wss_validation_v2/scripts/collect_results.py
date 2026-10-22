"""Collect only executed cases; comparisons never substitute for CFD solutions."""
from pathlib import Path
import csv,json
import numpy as np
from case_common import *
from analyze_vessel import csvout

def read(p):return list(csv.DictReader(Path(p).open()))
def aggregate(pattern,name):
 rows=[]
 for p in sorted(V.glob(pattern)):
  new=read(p)
  if name=='mesh_geometry.csv':
   case=p.parents[1]
   status='SKIPPED_BY_USER' if (case/'reports/user_cancellation.json').exists() else 'QUALIFIED_REUSED_CFD' if (case/'reports/baseline_reuse.json').exists() else json.loads((case/'reports/execution.json').read_text()).get('status','NOT_RUN') if (case/'reports/execution.json').exists() else 'NOT_YET_COMPLETED'
   for row in new:row['CFD_status']=status
  rows.extend(new)
 if rows:csvout(V/'data'/name,rows)
 return rows

def compare_fields(c1,c2):
 a=np.load(c1/'SV_MESH/mesh_arrays.npz');b=np.load(c2/'SV_MESH/mesh_arrays.npz')
 for k in ['points_m','tetra','boundary_triangles','facet_tags']:assert np.array_equal(a[k],b[k]),k
 x=a['points_m'];t=a['tetra'];wall=a['boundary_triangles'][a['facet_tags']==1];own=wss.boundary_owners(t,wall);_,area,n=wss.wall_geometry(x,t,wall,own);vol=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])/6
 fa=np.load(c1/'frozen_flow/flow_arrays_si.npz');fb=np.load(c2/'frozen_flow/flow_arrays_si.npz');ua=fa['velocity_m_s'];ub=fb['velocity_m_s'];pa=fa['pressure_pa'];pb=fb['pressure_pa']
 def l2(z):
  arr=z[t]
  if z.ndim==2:vv=np.sum(arr*arr,axis=(1,2))+np.sum(arr.sum(axis=1)**2,axis=-1)
  else:vv=np.sum(arr*arr,axis=1)+arr.sum(axis=1)**2
  return float(np.sqrt(np.dot(vol,vv/20)))
 wa=np.linalg.norm(wss.tangential_traction(wss.p1_gradients(x,t[own],ua),n,MU),axis=1);wb=np.linalg.norm(wss.tangential_traction(wss.p1_gradients(x,t[own],ub),n,MU),axis=1)
 return dict(reference_case=c1.name,comparison_case=c2.name,mesh_identity=True,velocity_volume_L2_relative_difference=l2(ub-ua)/l2(ua),pressure_volume_L2_relative_difference=l2(pb-pa)/l2(pa),wss_area_L2_relative_difference=float(np.sqrt(np.average((wb-wa)**2,weights=area)/np.average(wa**2,weights=area))),wss_mean_signed_difference_pct=float(100*(np.average(wb,weights=area)/np.average(wa,weights=area)-1)),wss_max_abs_difference_Pa=float(np.max(abs(wb-wa))))

aggregate('stage2/pipe_*/reports/pipe_validation.csv','pipe_validation.csv')
aggregate('stage[234]/*/reports/mesh_geometry.csv','mesh_geometry.csv')
regions=aggregate('stage[34]/*/reports/region_wss.csv','vessel_region_wss.csv')
for name in ['boundary_flows','internal_sections','adjacent_jump_summary','known_jump_pairs','outlet2_profile','fixed_baseline_mask_wss','threshold_overlap','cap_traction']:aggregate('stage[34]/*/reports/'+name+'.csv',name+'.csv')
coarse='vessel_baseline';medium='vessel_medium';cases={r['case'] for r in regions}
if {coarse,medium}.issubset(cases):
 rows=[];lookup={(r['case'],r['region']):r for r in regions}
 for region in sorted({r['region'] for r in regions if r['case']==coarse}):
  a,b=[lookup[(case,region)] for case in [coarse,medium]]
  for metric in ['mean_Pa','p05_Pa','p50_Pa','p95_Pa','min_Pa','max_Pa','below1_area_pct','below2_area_pct','above30_area_pct']:
   av,bv=[float(v[metric]) for v in [a,b]];delta=bv-av;rel=100*delta/av if av else ''
   rows.append(dict(region=region,metric=metric,unit='percentage_points' if 'area_pct' in metric else 'Pa',statistical_weight='wall_triangle_area',baseline_case=coarse,medium_case=medium,baseline_value=av,medium_value=bv,coarse_medium_absolute_change=delta,coarse_medium_relative_change_pct=rel,analysis='two_mesh_sensitivity_only',fine_CFD_status='SKIPPED_BY_USER',not_a_mesh_independence_or_error_estimate=True,reference='original_grid_for_relative_difference'))
 csvout(V/'data/vessel_mesh_comparison.csv',rows)
bcgroup=json.loads((V/'stage4/boundary_group.json').read_text()) if (V/'stage4/boundary_group.json').exists() else dict(baseline_case=medium,mesh=medium)
bcbase=bcgroup['baseline_case'];bcmesh=bcgroup['mesh']
if {bcbase,'O2_minus1pct','O2_plus1pct'}.issubset(cases):
 lookup={(r['case'],r['region']):r for r in regions};rows=[]
 for region in sorted({r['region'] for r in regions if r['case']==bcbase}):
  a,m,p=[lookup[(case,region)] for case in [bcbase,'O2_minus1pct','O2_plus1pct']]
  for metric in ['mean_Pa','p05_Pa','p50_Pa','p95_Pa','min_Pa','max_Pa','below1_area_pct','below2_area_pct','above30_area_pct']:
   av,mv,pv=[float(v[metric]) for v in [a,m,p]];dm=mv-av;dp=pv-av;den=max(abs(dm),abs(dp));rows.append(dict(mesh=bcmesh,baseline_case=bcbase,region=region,metric=metric,unit='percentage_points' if 'area_pct' in metric else 'Pa',baseline_value=av,minus1pct_value=mv,plus1pct_value=pv,minus_absolute_change=dm,plus_absolute_change=dp,minus_relative_change_pct=100*dm/av if av else '',plus_relative_change_pct=100*dp/av if av else '',signed_symmetry_defect=dm+dp,symmetry_defect_over_larger_response=abs(dm+dp)/den if den else 0,statistical_weight='wall_triangle_area',reference=bcbase+'_same_mesh_only_O2_changed'))
 csvout(V/'data/boundary_region_responses.csv',rows)
 flows=read(V/'data/boundary_flows.csv');lk={(r['case'],r['boundary']):r for r in flows};rr=[]
 for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
  a,m,p=[lk[(case,role)] for case in [bcbase,'O2_minus1pct','O2_plus1pct']]
  for metric in ['Q_outward_m3_s','flow_relative_to_inlet_pct']:
   av,mv,pv=[float(v[metric]) for v in [a,m,p]];rr.append(dict(mesh=bcmesh,baseline_case=bcbase,boundary=role,metric=metric,unit='m3/s' if metric.startswith('Q_') else 'percentage_points',baseline_value=av,minus1pct_value=mv,plus1pct_value=pv,minus_absolute_change=mv-av,plus_absolute_change=pv-av,minus_relative_change_pct=100*(mv/av-1),plus_relative_change_pct=100*(pv/av-1),signed_symmetry_defect=mv+pv-2*av,reference=bcbase+'_same_mesh_only_O2_changed',statistical_weight='exact_P1_boundary_triangle_flux'))
 csvout(V/'data/boundary_flow_responses.csv',rr)
for a,b,name in [(V/'stage2/pipe_nr16_cpu_mpi8',V/'stage2/pipe_nr16_cpu_mpi8_halfdt','pipe_timestep_comparison.json'),(V/'stage3/vessel_baseline',V/'stage3/vessel_baseline_cpu_mpi8','vessel_backend_comparison.json'),(V/'stage3/vessel_baseline',V/'stage3/vessel_baseline_gpu_mpi1_halfdt','vessel_timestep_comparison.json')]:
 if any((a/'reports'/n).exists() for n in ['pipe_validation.csv','region_wss.csv']) and any((b/'reports'/n).exists() for n in ['pipe_validation.csv','region_wss.csv']):dump(V/'data'/name,compare_fields(a,b))
print('Collected available executed-case evidence.')
