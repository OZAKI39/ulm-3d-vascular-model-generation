#!/usr/bin/env python3
"""Read saved production paths and re-solve failed FINAL states, without advancing time."""
from collections import Counter
import argparse
import numpy as np
from analyze_particle9a1 import ROOT,R,D,O,read,csvwrite,summarize,fixed
from particle_3d.particle81_simulation import environment,dump
from particle_3d.particle_shapes import Sphere
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle9a_motion import augment_planar_system
from particle_3d.wall_gap import wall_gap
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.resistance_solver import solve_resistance

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--partial',action='store_true');args=parser.parse_args()
 if not args.partial and not (O/'provenance/production_completed.json').exists():raise ValueError('Full production audit requires completed batch receipt')
 prefix='production_partial' if args.partial else 'production'
 env=environment();ledger=read(O/'admission/birth_ledger.json');policy=NearFieldRegularizationV1();rows=[];failed=[]
 for event in ledger['events']:
  pid=event['particle_id'];meta=O/'P9A1/trajectories'/f'mb_{pid:06d}.json'
  if not meta.exists():continue
  m=read(meta);a=np.load(meta.with_suffix('.npz'))['samples'];row=summarize(m,a,[],'P9A1',env)
  # No diagnostic histories inferred when production logging is off.
  for k in ['contact_max_condition','duplicate_constraints_dropped','handoff_events','rejection_reasons']:row.pop(k,None)
  if not m['completed']:
   x=a[-1,1:4];s={pid:Sphere(x,m['radius_m'])};field=env.field.sample(x);record=dict(particle_id=pid,failure_detail=m['failure_detail'],time_s=float(a[-1,0]),position=x.tolist(),inside_lumen=bool(field.inside_lumen),state_role='FINAL_SAVED_STATE_REASSEMBLY_NO_TIME_ADVANCE',supported_stationary_contact=False)
   if field.inside_lumen:
    u={pid:np.r_[field.velocity_m_s,.5*field.vorticity_s_inv]};base=assemble_v1(s,u,env.mu,env.wall,policy=policy)
    sys=augment_planar_system(base,s,env.mu,env.wall,lambda _:field.velocity_gradient_s_inv,wall_gap);cs,cr=handoff_constraints(s,env.wall,policy)
    try:
     sol=solve_resistance(sys,constraints=cs);rank=sol.record.get('contact_redundancy',{}).get('rank_after',0);budget=sol.record.get('contact_kkt',{}).get('velocity_budget_m_s',0.)
     stationary='STATIONARY' in str(m['failure_detail']) or 'ROUNDOFF_SCALE_STAGNATION' in str(m['failure_detail'])
     record.update(solver=sol.record,contacts=cr,free=u[pid],solved=sol.velocity,rank=rank,
       supported_stationary_contact=bool(stationary and rank==3 and np.max(abs(sol.velocity[:3]))<=budget and np.all(np.array(sol.record['multipliers'])>=0)))
    except Exception as err:record['solve_error']=str(err)
   failed.append(record);row['stationary_contact_supported']=record['supported_stationary_contact']
  rows.append(row)
 csvwrite(D/(prefix+'_outcomes.csv'),rows);dump(D/(prefix+'_failed_final_states.json'),failed)
 result=dict(generated_count=len(rows),completed_count=sum(r['completed'] for r in rows),outlet_distribution=dict(Counter(r['outlet'] or 'NO_EXIT' for r in rows)),
  outcome_counts=dict(Counter(r['outcome'] for r in rows)),failure_counts=dict(Counter(str(r['failure_detail']) for r in rows if not r['completed'])),
  failure_ids=[r['particle_id'] for r in rows if not r['completed']],supported_stationary_failure_count=sum(f['supported_stationary_contact'] for f in failed),
  unresolved_failure_ids=[f['particle_id'] for f in failed if not f['supported_stationary_contact']],
  penetration_count=sum(r['minimum_gap_m'] < -env.wall.roundoff_m for r in rows),handoff_violation_count=sum(r['minimum_g_nf_m'] < -env.wall.roundoff_m for r in rows),
  all_500_present=len(rows)==500,classification_role='ALL_INCOMPLETE_PATHS_RETAINED; STATIONARITY_IS_NOT_PHYSIOLOGICAL_TRAPPING_PROOF')
 result['outlet_distribution']={role:result['outlet_distribution'].get(role,0) for role in ['OUTLET_01','OUTLET_02','OUTLET_03','NO_EXIT']}
 result['outcome_counts']={role:result['outcome_counts'].get(role,0) for role in ['COMPLETED','INLET_ESCAPE','HANDOFF_OR_SAFETY_STOP','TIME_LIMIT','OTHER','OUTSIDE_OTHER']}
 dump(D/(prefix+'_summary.json'),result);print(result)
if __name__=='__main__':main()
