#!/usr/bin/env python3
"""Read-only P9-A.1 gates and independently regenerable figure source tables."""
import ast,csv,gzip,importlib.util,json,types,xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
import numpy as np
from scipy import sparse
from particle_3d.particle81_simulation import environment,dump
from particle_3d.resistance_assembly import ResistanceSystem
from particle_3d.resistance_solver import solve_resistance
from particle_3d.kinematic_contact import ContactConstraint
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block,legacy_planar_coefficients_si
from particle_3d.validation_boundary import ValidationBoundaryClassifier
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'particle_3d/reports/particle9a1_2mmps';D=R/'data';O=ROOT/'particle_3d/outputs/particle9a1_2mmps';OLD=ROOT/'particle_3d/reports/particle9a_2mmps_diagnosis/data';IDS=[3,5,7,10,11,12,14,15,17,18,19,22]

def read(p):return json.loads(p.read_text())
def lines(p):
 with gzip.open(p,'rt') as f:return [json.loads(x) for x in f]
def csvwrite(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def original(model,pid):
 stem=OLD/'replays'/model/'trajectories'/f'mb_{pid:06d}'
 return read(stem.with_suffix('.json')),np.load(stem.with_suffix('.npz'))['samples']
def fixed(model,pid):
 stem=O/model/'trajectories'/f'mb_{pid:06d}'
 return read(stem.with_suffix('.json')),np.load(stem.with_suffix('.npz'))['samples'],lines(stem.with_suffix('.audit.jsonl.gz')) if stem.with_suffix('.audit.jsonl.gz').exists() else []
def category(m,a,env):
 if m['completed']:return 'COMPLETED'
 detail=str(m.get('failure_detail'))
 if 'CENTER_OUTSIDE' in detail:
  c=ValidationBoundaryClassifier({'INLET':env.boundaries['INLET']})
  for x,y in zip(a[:-1,1:4],a[1:,1:4]):
   hit=c.first_event(x,y)
   if hit is not None:return 'INLET_ESCAPE'
  return 'OUTSIDE_OTHER'
 if 'SUBDIVISION' in detail or 'REFINEMENT' in m['end_reason'] or 'GUARD' in detail:return 'HANDOFF_OR_SAFETY_STOP'
 if m['end_reason']=='PHYSICAL_RESIDENCE_HORIZON_REACHED':return 'TIME_LIMIT'
 return 'OTHER'
def canonical():
 spec=importlib.util.spec_from_file_location('old',R/'reference/old_planar_wall_hydrodynamics.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 tree=ast.parse((ROOT/'particle_3d/reports/particle9a_2mmps/reference/legacy_particle_mobility.py').read_text());tree.body=[n for n in tree.body if not isinstance(n,ast.Try)]
 for node in ast.walk(tree):
  if isinstance(node,ast.FunctionDef):node.decorator_list=[]
 legacy=types.ModuleType('legacy');exec(compile(tree,'legacy','exec'),legacy.__dict__)
 n=np.array([.6,0,.8]);t=np.array([-.8,0,.6]);gamma=500.;mu=.00345312;a=1e-6;g=gamma*np.outer(t,n);rows=[]
 for xi in np.geomspace(.001,1,61):
  u=np.r_[gamma*a*(1+xi)*t, .5*gamma*np.cross(n,t)]
  before=old.planar_wall_affine_block(a,mu,xi*a,n,g,u)[2];after=planar_wall_affine_block(a,mu,xi*a,n,g,u)[2]
  vx,vz,wy,w,_=legacy.background_hydrodynamic_velocity_scalar(u[0]*1e6,u[2]*1e6,g[0,0],g[0,2],g[2,0],g[2,2],n[0],n[2],a*1e6,xi*a*1e6,mu)
  fs,ty,_=legacy.wall_shear_load_scalar(mu,a*1e6,xi*a*1e6,gamma);err=legacy_planar_coefficients_si(xi)['projection_entry_error'];base=1/(6*np.pi*mu*a)
  rows.append(dict(gap_ratio=xi,old_Vt_m_s=np.array(before['target_tangential_velocity_xyz'])@t,new_Vt_m_s=np.array(after['target_tangential_velocity_xyz'])@t,reference_2d_Vt_m_s=(vx*t[0]+vz*t[2])*1e-6,
   old_omega_s_inv=before['target_tangential_omega_xyz'][1],new_omega_s_inv=after['target_tangential_omega_xyz'][1],reference_2d_omega_s_inv=wy,
   reciprocity_projection_Vt_bound_m_s=w*base*err*abs(ty*1e-18/a)+1e-15,reciprocity_projection_omega_bound_s_inv=w*base*err*abs(fs*1e-12)/a+1e-9))
 csvwrite(D/'canonical.csv',rows)
 result=dict(max_new_old_Vt_m_s=max(abs(r['new_Vt_m_s']-r['old_Vt_m_s']) for r in rows),max_new_old_omega_s_inv=max(abs(r['new_omega_s_inv']-r['old_omega_s_inv']) for r in rows),
  max_new_raw_2d_Vt_m_s=max(abs(r['new_Vt_m_s']-r['reference_2d_Vt_m_s']) for r in rows),max_new_raw_2d_omega_s_inv=max(abs(r['new_omega_s_inv']-r['reference_2d_omega_s_inv']) for r in rows),
  old_parity_pass=all(abs(r['new_Vt_m_s']-r['old_Vt_m_s'])<1e-15 and abs(r['new_omega_s_inv']-r['old_omega_s_inv'])<1e-9 for r in rows),
  raw_2d_with_published_reciprocity_projection_pass=all(abs(r['new_Vt_m_s']-r['reference_2d_Vt_m_s'])<=r['reciprocity_projection_Vt_bound_m_s'] and abs(r['new_omega_s_inv']-r['reference_2d_omega_s_inv'])<=r['reciprocity_projection_omega_bound_s_inv'] for r in rows))
 dump(D/'canonical_summary.json',result);return result

def contact_fixture():
 f=read(R/'reference/id7_duplicate_contact.json');sys=ResistanceSystem((7,),sparse.csr_matrix(f['R']),np.array(f['b']),np.array(f['self_diagonal']),np.array(f['free']),[],[])
 cs=[ContactConstraint(tuple(r['canonical_id']),7,None,r['normal'],r['particle_point_m']) for r in f['contacts']];s=solve_resistance(sys,constraints=cs);record=s.record
 dump(D/'id7_fixed_fixture.json',dict(source='reference/id7_duplicate_contact.json',velocity=s.velocity,solver=record))
 rows=[dict(model='Original',constraint_count=2,retained_count=2,numerical_rank=2,rank_definition='OLD_DUAL_LINEAR_EPS_RANK',condition=f['original']['contact_kkt']['condition'],minimum_normal_speed_m_s=min(f['original']['normal_speeds_m_s'])),
 dict(model='Fixed',constraint_count=2,retained_count=record['contact_redundancy']['retained_constraint_count'],numerical_rank=record['contact_redundancy']['rank_after'],rank_definition='GRAM_RESOLUTION_SVD',condition=record['contact_kkt']['condition'],minimum_normal_speed_m_s=min(record['normal_speeds_m_s']))]
 csvwrite(D/'contact_id7.csv',rows);return record

def summarize(m,a,rows,model,env):
 rejected=Counter(r.get('error') for r in rows if not r['accepted']);contact=[r['solver'] for r in rows if r.get('solver',{}).get('contact_count')]
 return dict(particle_id=m['particle_id'],model=model,outcome=category(m,a,env),completed=m['completed'],outlet=m['exit_outlet'],failure_detail=m['failure_detail'],
  elapsed_s=float(a[-1,0]),path_length_m=m['path_length_m'],mean_speed_m_s=m['path_length_m']/a[-1,0] if a[-1,0] else 0.,minimum_gap_m=float(a[:,14].min()),minimum_g_nf_m=float(a[:,15].min()),
  accepted_steps=m['accepted_steps'],rejected_trials=m['rejected_trials'],rejection_reasons=dict(rejected),handoff_events=sum(r.get('handoff_event') is not None for r in rows),
  contact_max_condition=max([c['contact_kkt']['condition'] for c in contact],default=1.),duplicate_constraints_dropped=sum(len(c.get('contact_redundancy',{}).get('dropped',[])) for c in contact),
  wall_seconds=m['wall_seconds'],position_identity_max_m=m.get('maximum_position_identity_error_m'))

def main():
 env=environment();c=canonical();contact=contact_fixture();allrows=[];new={};arrs={};audits={}
 for pid in IDS:
  for model in ['P65','P9A']:
   m,a=original(model,pid);allrows.append(summarize(m,a,[],model,env))
  for model in ['P9A1','P65_NEW']:
   m,a,rows=fixed(model,pid);summary=summarize(m,a,rows,model,env);allrows.append(summary)
   if model=='P9A1':new[pid]=summary;arrs[pid]=a;audits[pid]=rows
 csvwrite(D/'same12.csv',allrows);dump(D/'same12.json',allrows)
 inlet=[]
 for pid in [3,15]:
  old=lines(OLD/'replays/P9A/trajectories'/f'mb_{pid:06d}.trials.jsonl.gz')[0];r=audits[pid][0]
  # Exact anchor cap triangle normal, oriented by actual FEM. No fitted-plane gate.
  event=read(O/'admission/birth_ledger.json')['events'][IDS.index(pid)]
  surf=env.boundaries['INLET'];tri=np.array(surf.points)[surf.faces.reshape(-1,4)[event['anchor_triangle'],1:]]
  n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n)
  free=np.array(r['free'][:3]);n*=1 if n@free>0 else -1
  wn=np.array(old['wall_normal_xyz']);v=np.array(r['velocity'][:3])
  inlet.append(dict(particle_id=pid,FEM_inward_m_s=float(free@n),old_inward_m_s=float(np.array(old['velocity_xyz'])@n),new_inward_m_s=float(v@n),
   bulk_tangential_m_s=old['bulk_t_speed'],old_target_tangential_m_s=old['target_tangential_speed'],new_solved_tangential_m_s=float(np.linalg.norm(v-wn*(wn@v))),
   topology=r['planar'][0].get('open_rim_fallback'),anchor_triangle=event['anchor_triangle'],normal=n.tolist(),new_outcome=new[pid]['outcome']))
 csvwrite(D/'inlet_3_15.csv',inlet);dump(D/'inlet_3_15.json',inlet)
 handoff=[];paths=[]
 for pid in IDS:
  a=arrs[pid]
  for row in a:paths.append(dict(particle_id=pid,time_s=row[0],x_m=row[1],y_m=row[2],z_m=row[3],outcome=new[pid]['outcome']))
  if pid in [12,14,17]:
   for model,data in [('P9A',original('P9A',pid)[1]),('P9A1',a)]:
    for row in data:handoff.append(dict(particle_id=pid,model=model,time_s=row[0],g_nf_m=row[15],accepted_dt_s=row[18],dt_over_nominal=row[18]/.00025))
 csvwrite(D/'handoff.csv',handoff);csvwrite(D/'real_paths.csv',paths)
 tests={}
 for name in ['local_tests','remote_tests','shared_core_tests']:
  xml=R/'logs'/f'{name}.xml'
  if xml.exists():
   suites=ET.parse(xml).getroot().iter('testsuite');tests[name]=all(int(s.get('failures','0'))==0 and int(s.get('errors','0'))==0 and int(s.get('tests','0'))>0 for s in suites)
  else:tests[name]=False
 handoff_pass=all(not('SUBDIVISION' in str(new[i]['failure_detail']) and any('HANDOFF_ENDPOINT_BELOW_LOWER' in k for k in new[i]['rejection_reasons'])) for i in [12,14,17])
 gates=dict(unit_regressions=all(tests.values()),canonical=c['old_parity_pass'] and c['raw_2d_with_published_reciprocity_projection_pass'],p65_normal=all(tests.values()),
  id3_no_inlet_escape=new[3]['outcome'] not in ['INLET_ESCAPE','OUTSIDE_OTHER'] and inlet[0]['new_inward_m_s']>0,
  id15_no_inlet_escape=new[15]['outcome'] not in ['INLET_ESCAPE','OUTSIDE_OTHER'] and inlet[1]['new_inward_m_s']>0,
  id22_no_old_collapse=new[22]['mean_speed_m_s']>1e3*next(r['mean_speed_m_s'] for r in allrows if r['particle_id']==22 and r['model']=='P9A'),
  id7_contact_rank=contact['contact_redundancy']['retained_constraint_count']==1 and contact['contact_kkt']['condition']<2 and new[7]['contact_max_condition']<1e10,
  old_handoff_chatter_removed=handoff_pass,
  no_wall_penetration=all(r['minimum_gap_m']>=-env.wall.roundoff_m and r['minimum_g_nf_m']>=-env.wall.roundoff_m for r in allrows if r['model'] in ['P9A1','P65_NEW']),
  no_new_solver_pathology=all(not any(x in str(r['failure_detail']) for x in ['CONTACT_PRIMAL','ILL_CONDITIONED','Dual','dual','HANDOFF_EVENT_COUNT_LIMIT']) for r in allrows if r['model'] in ['P9A1','P65_NEW']))
 samepass=all(gates.values());smoke=[]
 if (O/'provenance/smoke30_completed.json').exists():
  ledger=read(O/'admission/birth_ledger.json')
  for event in ledger['events'][:30]:
   m,a,rows=fixed('P9A1',event['particle_id']);smoke.append(summarize(m,a,rows,'P9A1',env))
  csvwrite(D/'smoke30.csv',smoke)
  # Completion is not a user-required numerical gate. A previously established
  # stationary hard-contact state is retained as incomplete, never relabelled
  # completed. New failures block escalation. This policy does not select events.
  stationary=[]
  for r in smoke:
   if r['completed']:continue
   previous=new.get(r['particle_id'])
   _,_,audit=fixed('P9A1',r['particle_id'])
   last=[x for x in audit if x['accepted'] and 'solver' in x][-32:]
   verified=bool(previous and previous['failure_detail']==r['failure_detail'] and
       'EXACT_STATIONARY_COORDINATES_32' in str(r['failure_detail']) and len(last)==32)
   for x in last:
    s=x['solver'];k=s.get('contact_kkt',{});rank=s.get('contact_redundancy',{}).get('rank_after',0)
    verified=verified and rank==3 and bool(np.all(np.array(s.get('multipliers',[]))>0)) and max(abs(np.array(x['velocity'][:3])))<=k.get('velocity_budget_m_s',0)
   stationary.append(dict(particle_id=r['particle_id'],same12_failure_unchanged=previous is not None,
       verified_independent_three_contact_stationarity=bool(verified),last_solve=last[-1] if last else None))
  dump(D/'stationary_stop_audit.json',stationary)
  gates['smoke30_no_systematic_failure']=all(x['verified_independent_three_contact_stationarity'] for x in stationary) and all(
      r['minimum_gap_m']>=-env.wall.roundoff_m and r['minimum_g_nf_m']>=-env.wall.roundoff_m and r['contact_max_condition']<1e10 for r in smoke)
 production_allowed=samepass and len(smoke)==30 and gates.get('smoke30_no_systematic_failure',False)
 result=dict(gates=gates,test_logs=tests,same12_pass=samepass,production_authorized_by_all_gates=production_allowed,
  automated_status='P9A1_AUTOMATED_VALIDATION_PASS' if all(gates.values()) else 'P9A1_AUTOMATED_VALIDATION_FAIL',human_status='PENDING_USER_REVIEW',
  completed_counts={m:sum(r['completed'] for r in allrows if r['model']==m) for m in ['P65','P9A','P65_NEW','P9A1']},
  outcome_counts={m:dict(Counter(r['outcome'] for r in allrows if r['model']==m)) for m in ['P65','P9A','P65_NEW','P9A1']},
  canonical=c,id22_old_mean_speed_m_s=next(r['mean_speed_m_s'] for r in allrows if r['particle_id']==22 and r['model']=='P9A'),id22_new_mean_speed_m_s=new[22]['mean_speed_m_s'],
  smoke30_status='EXECUTED' if smoke else 'NOT_LAUNCHED',smoke30_outcomes=dict(Counter(r['outcome'] for r in smoke)),production_status='EXECUTED' if (O/'provenance/production_completed.json').exists() else ('RUNNING' if (O/'provenance/production_config.json').exists() else 'NOT_LAUNCHED'),
  contact_id7_fixture_condition_old=read(R/'reference/id7_duplicate_contact.json')['original']['contact_kkt']['condition'],contact_id7_fixture_condition_new=contact['contact_kkt']['condition'])
 dump(D/'gates.json',result);print((D/'gates.json').read_text())
if __name__=='__main__':main()
