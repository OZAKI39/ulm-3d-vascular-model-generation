#!/usr/bin/env python3
"""Read saved diagnostic continuations with a curved-flow local passage criterion.

The original fixed final-FEM direction projection is retained as exploratory evidence.
It is not an appropriate downstream coordinate after a bend. No path is reintegrated.
"""
from pathlib import Path
from copy import deepcopy
import numpy as np
from particle_3d.routing_stationary_audit import *
from particle_3d.particle8_replay import REPO
from particle_3d.convex_triangle import triangle_closest_many
R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';D=R/'data'

def main():
 env=environment();source=D/'stationary_audit_fixed_direction_reference.json'
 if not source.exists():dump(source,read(D/'stationary_audit.json'))
 original=read(source);result=[]
 policy=dict(role='DIAGNOSTIC_ONLY',revision='LOCAL_CURVED_FLOW_PASSAGE_V2',rationale='A fixed upstream FEM direction is not a downstream coordinate after a vessel bend.',
     passage='center outside ball of 3 original radii around original jam, cumulative local-FEM downstream displacement >3 original radii, five consecutive moving nominal steps, positive clearance from every original contact triangle',
     outlet_alternative='original triangle classifier verifies terminal outlet AND original contact triangles are cleared AND cumulative local downstream displacement >3 radii',
     original_paths_retained=True,no_reintegration=True,no_production_change=True)
 for old in original:
  r=deepcopy(old);final=np.array(r['final_position_m']);radius=r['radius_m'];pid=r['particle_id'];trials=[];supplement=[];already_passed=False
  for t in r['trials']:
   file=R/'diagnostic_outputs/virtual_radius'/Path(t['trajectory_path']).name;a=np.load(file)['samples'];cumulative=0.;run=0;passage=None;previous=0
   # Nominal outputs preserve the exact same dt grid used by the original continuation.
   ids=[i for i in range(1,len(a)) if abs(a[i,0]/DT-round(a[i,0]/DT))<64*np.finfo(float).eps*max(1,abs(a[i,0]/DT))]
   if len(a)-1 not in ids:ids.append(len(a)-1)
   for i in ids:
    p=a[previous,1:4];x=a[i,1:4];delta=x-p;f=env.field.sample(p);v=f.velocity_m_s;speed=np.linalg.norm(v)
    progress=float(delta@v/speed) if speed else 0.;cumulative+=progress;run=run+1 if progress>0 and np.linalg.norm(delta)>16*env.wall.roundoff_m else 0;previous=i
    if cumulative<=3*radius:continue
    near,_=triangle_closest_many(x,env.wall.triangles[r['contact_triangle_ids']]);gaps=np.linalg.norm(x-near,axis=1)-t['virtual_radius_m'];clear=bool(np.all(gaps>2e-9+16*env.wall.roundoff_m))
    radial=float(np.linalg.norm(x-final));hit=env.classifier.first_event(a[i-1,1:4],x) if i>0 else None
    if clear and ((radial>3*radius and run>=5) or (hit is not None and hit.role.startswith('OUTLET_'))):
     passage=dict(sample_index=i,time_s=float(a[i,0]),radial_distance_m=radial,cumulative_local_downstream_m=cumulative,consecutive_moving_nominal_steps=run,original_contact_gaps_m=gaps.tolist(),outlet=hit.role if hit else None);break
   t['original_fixed_direction_status']=t['status'];t['local_passage_witness']=passage;t['passage_criterion_revision']=policy['revision'];t['verified_cumulative_local_downstream_m']=cumulative
   if passage:t['status']='PASSED_LOCAL_HOTSPOT'
   if already_passed:supplement.append(t)
   else:
    trials.append(t)
    if t['radius_ratio']<1 and passage:already_passed=True
  r['trials']=trials;r['supplementary_trials_already_computed_under_fixed_direction_rule']=supplement
  passes=[t for t in trials if t['radius_ratio']<1 and t['status']=='PASSED_LOCAL_HOTSPOT'];unresolved=any('UNRESOLVED' in t['status'] for t in trials)
  category='NUMERICAL_OR_GEOMETRIC_UNRESOLVED';sensitivity='UNRESOLVED';critical=None
  consistent=r['reference_jam_reproduced'] and r['retained_rows_independent'] and r['downstream']['feasible_downstream_direction'] is False
  if consistent and not unresolved and passes:
   p=passes[0];ratio=p['radius_ratio'];fail=trials[trials.index(p)-1]['radius_ratio'];sensitivity='HIGHLY_RADIUS_SENSITIVE' if ratio>=.95 else 'INTERMEDIATE';category='RIGID_MODEL_HIGH_SENSITIVITY' if ratio>=.95 else 'RIGID_GEOMETRIC_JAM';critical=dict(highest_tested_passing=ratio,lowest_tested_failing=fail,interval=[ratio,fail],role='TESTED_LOCAL_PASSAGE_BRACKET_NOT_EXACT_THRESHOLD; MONOTONICITY_NOT_PROVEN')
  elif consistent and not unresolved and all(t['status']=='STATIONARY' for t in trials):
   category='RIGID_GEOMETRIC_JAM';sensitivity='ROBUST_RIGID_SIZE_EXCLUSION';critical=dict(highest_tested_passing=None,lowest_tested_failing=.9,interval=None,role='NO_PASSAGE_AT_TESTED_RATIOS')
  r.update(classification=category,radius_sensitivity=sensitivity,critical_radius_ratio=critical,passage_criterion_revision=policy['revision']);result.append(r)
  print(pid,category,[(t['radius_ratio'],t['status']) for t in trials])
 dump(D/'radius_passage_policy.json',policy);dump(D/'stationary_audit.json',result)
if __name__=='__main__':main()
