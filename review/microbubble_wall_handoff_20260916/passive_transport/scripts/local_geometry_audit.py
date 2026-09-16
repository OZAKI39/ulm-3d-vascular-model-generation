"""Reconstruct every accepted RK2 step independently and audit original STL geometry."""
from pathlib import Path
import sys,json,itertools
import numpy as np
R=Path(sys.argv[1]);sys.path.insert(0,str(R/'src'))
from finalize_passive_transport_v0 import read,trajectory,xyz,vel,integrity
from reference_passive_transport import ReferenceField,WallReference,velocity
assert integrity(R)
ct=read(R/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json');selection=read(R/'contracts/PASSIVE_TRANSPORT_TIMESTEP_CONTRACT.json');label='c025' if selection['selected_c_adv']==.25 else 'c0125';field=ReferenceField(R/'fields/FROZEN_FLOW_FIELD_V0.h5');wall=WallReference(R/'provenance/closed_geometry_m.stl')
names=['case2_'+s for s in ['c050','c025','c0125']]
if (R/'cases'/('case5_'+label)/'RUN_METRICS.json').exists():names.append('case5_'+label)
results=[]
def pair_min(a,b,radii):
 best=float('inf');pair=None
 for i,j in itertools.combinations(range(len(a)),2):
  rel=a[i]-a[j];change=(b[i]-a[i])-(b[j]-a[j]);vv=np.dot(change,change);t=np.clip(-np.dot(rel,change)/vv,0.,1.) if vv>0 else 0.;gap=np.linalg.norm(rel+t*change)-radii[i]-radii[j]
  if gap<best:best=float(gap);pair=[i+1,j+1]
 return best,pair
for name in names:
 sp=read(R/'cases'/name/'CASE_CONTRACT.json');t=trajectory(R/'cases'/name);end=int(t['step'].max());x=np.array(sp['initial_positions_m']);radii=np.array(sp['diameters_um'])*.5e-6;dt=sp['dt_s'];threshold=radii+3*ct['dx_m'];max_x=0.;max_v=0.;minimum=float('inf');minimum_pair=float('inf');stored={int(step):xyz(t[t['step']==step]) for step in np.unique(t['step'])};max_distance_error=0.;segments_safe=True;midpoints_safe=True;initial_inside=bool(wall.enclosed(x).all());initial_gap=np.array([wall.gap(q) for q in x])-radii
 print('REPLAY',name,'STEPS',end,flush=True)
 for step in range(end+1):
  distances=np.array([wall.gap(q) for q in x]);minimum=min(minimum,float(np.min(distances-radii)))
  if step in stored:
   rows=t[t['step']==step];max_x=max(max_x,float(np.max(np.linalg.norm(x-stored[step],axis=1))));max_v=max(max_v,float(np.max(np.linalg.norm(velocity(field,x)-vel(rows),axis=1))));cpp_d=rows['center_wall_distance_um']*1e-6;actual_d=np.array([wall.gap(q) for q in stored[step]]);max_distance_error=max(max_distance_error,float(np.max(np.abs(cpp_d-actual_d))))
  if step==end:break
  k1=velocity(field,x);mid=x+.5*dt*k1;k2=velocity(field,mid);new=x+dt*k2;velocity(field,new)
  for i in range(len(x)):
   midpoints_safe &= wall.gap(mid[i])>=threshold[i];segments_safe &= wall.segment_safe(x[i],new[i],threshold[i])
  if len(x)>1:
   gap,pair=pair_min(x,new,radii);midgap,_=pair_min(mid,mid,radii);minimum_pair=min(minimum_pair,gap,midgap)
  x=new
 # Independently explain the halted next step, without accepting it.
 reason=read(R/'cases'/name/'trajectory_state_rank0.json')['termination_reason'];next_evidence={}
 if reason!='MAX_PHYSICAL_TIME':
  mid=x+.5*dt*velocity(field,x);new=x+dt*velocity(field,mid);velocity(field,new)
  if reason=='SPHERE_OVERLAP_SAFETY_STOP':
   gap,pair=pair_min(x,new,radii);next_evidence={'predicted_next_segment_min_pair_gap_m':gap,'pair_LAMMPS_ids':pair,'pair_source_bubble_ids':[sp['source_bubble_ids'][i-1] for i in pair]};assert gap<0
  elif reason=='WALL_SAFETY_STOP':
   stops=[wall.gap(mid[i])<threshold[i] or not wall.segment_safe(x[i],new[i],threshold[i]) for i in range(len(x))];assert any(stops);next_evidence={'wall_safety_violation_reproduced':True}
  else:raise RuntimeError('Unexpected terminal reason '+reason)
 accepted_inside=bool(wall.enclosed(xyz(t)).all());checks={'initial_inside':initial_inside,'all_output_centers_inside':accepted_inside,'initial_margin_5dx':bool(np.all(initial_gap>=5*ct['dx_m'])),'all_internal_wall_gaps_ge_3dx':minimum>=3*ct['dx_m'],'all_midpoints_safe':bool(midpoints_safe),'all_continuous_segments_safe':bool(segments_safe),'Cplusplus_distance_vs_VTK':max_distance_error<=ct['gates']['geometry_distance_m'],'independent_RK2_replay':max_x<=1e-12 and max_v<=ct['gates']['cpp_python_velocity_m_s']}
 if len(x)>1:checks['no_overlap_any_accepted_segment_or_midpoint']=minimum_pair>=0
 r={'name':name,'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'accepted_steps_replayed':end,'max_independent_position_difference_m':max_x,'max_independent_velocity_difference_m_s':max_v,'max_wall_distance_difference_m':max_distance_error,'min_internal_surface_wall_gap_m':minimum,'min_accepted_pair_gap_m':minimum_pair if len(x)>1 else None,'termination_reason':reason,'next_step_rejection_evidence':next_evidence};results.append(r);print(json.dumps(r),flush=True)
result={'status':'PASS' if all(r['status']=='PASS' for r in results) else 'FAIL','method':'Independent Python trilinear + RK2 replay of every internal step, VTK original STL distances, enclosure, Lipschitz segment certificate and swept pair distances. No C++ calls.','cases':results}
(R/'LOCAL_GEOMETRY_AND_REPLAY_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');assert result['status']=='PASS'
