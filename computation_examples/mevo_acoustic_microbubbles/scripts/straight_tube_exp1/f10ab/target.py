"""Finite preregistered fixed-shape target, with discovery and confirmation separated."""
import copy,numpy as np,h5py
from .common import *
from .plan import task,freeze
from .analyze import numerical

def plan():
 rep=read(BASE/'analysis/REPLICATION.json');g=next(g for g in rep['groups'] if g['shape']=='S2' and g['action']=='A5_X180')
 assert g['full_seed_numerical_confirmation'] and g['signed_xy_direction_repeated']
 protocol=dict(target='T2: increase fixed-ID 3D aspect ratio by 3 percent in 2 ms',
  initial_state='S2_r0',target_aspect_fraction=.03,aspect_tolerance_fraction=.005,
  com_limit_m=20e-9,require_final_single_component=True,require_original_numerical_gates=True,
  cost='((aspect_fraction-0.03)/0.005)^2 + (COM_norm/20nm)^2; select legal completed candidates, report connectivity independently and include in final admission',
  candidate_amplitudes_Pa_equivalent=[600,700,800,900],phase_deg=[0,180],
  selected_action_confirmation='independent rerun S2_r0; both S2_r1 and S2_r2, never selected by response; all three receive dt/cache controls',
  holdout_scope='r1/r2 were mechanism confirmation states, not used to tune target amplitude; not a new population-wide external test',
  frozen_before_new_target_runs=True,identity=identity())
 assert not (BASE/'TARGET_PREREGISTRATION.json').exists();write(BASE/'TARGET_PREREGISTRATION.json',protocol)
 tasks=[]
 for amp in protocol['candidate_amplitudes_Pa_equivalent']:
  rec=copy.deepcopy(read(BASE/'ACTIONS.json')['A5_X180']);rec['segments'][0]['amplitude_Pa_equivalent'][2:]=[amp,amp]
  tasks.append(task('target_search','S2_r0',f'T_A{amp}',.002,schedule=rec))
 freeze('target_search',tasks,150,'Predeclared +3% aspect target, 20nm COM, endpoint connectivity. Four amplitudes only, phase180. No hand-adjusted particle positions.')

def score(r):
 a=r['initial_shape'];m=r['final_shape'];p=read(BASE/'TARGET_PREREGISTRATION.json')
 frac=m['aspect_3d']/a['aspect_3d']-1;drift=np.linalg.norm(np.asarray(m['center_3d_m'])-a['center_3d_m'])
 return dict(run_id=r['run_id'],state_id=r['state_id'],aspect_fraction=frac,aspect_error=frac-p['target_aspect_fraction'],com_m=drift,
  final_components=m['proximity_components'],final_largest=m['largest_component'],similarity_nm=m['similarity_rmsd_3d_m']*1e9,
  cost=((frac-.03)/.005)**2+(drift/20e-9)**2,
  target_metrics_pass=bool(abs(frac-.03)<=.005 and drift<=20e-9 and m['proximity_components']==1))

def confirm_plan():
 rows=read(BASE/'target_search/INDEX.json');completed=[r for r in rows if r['status']=='PASS'];assert completed
 best=min(completed,key=lambda r:score(r)['cost']);selection=dict(selected=best,all_candidates=[score(r) for r in completed],identity=identity(),selection='minimum preregistered cost; no seed cherry-picking')
 write(BASE/'TARGET_SELECTION.json',selection)
 tasks=[task('target_confirmation',state,best['action'],.002,v,best['task']['schedule']) for state in ['S2_r0','S2_r1','S2_r2'] for v in ['DT1','DT05','FINE']]
 freeze('target_confirmation',tasks,250,'Independent recompute of chosen action; both mechanism holdouts; dt/cache for each, same target tolerances.')

def evaluate():
 rows=read(BASE/'target_confirmation/INDEX.json');comparisons=[];cases=[]
 for a in rows:
  if a.get('variant')!='DT1':continue
  cc=[]
  for b in rows:
   if b.get('state_id')==a['state_id'] and b.get('variant') in ['DT05','FINE']:
    if a['status']=='PASS' and b['status']=='PASS':c=numerical(a,b)
    else:c=dict(reference=a['run_id'],comparison=b['run_id'],original=dict(admitted=False),shape_uncertainty={})
    comparisons.append(c);cc.append(c)
  d=score(a);d['numerical_pass']=len(cc)==2 and all(c['original']['admitted'] for c in cc)
  d['admitted']=d['target_metrics_pass'] and d['numerical_pass'];d['trajectory_path']=a.get('trajectory_path');cases.append(d)
 best=read(BASE/'TARGET_SELECTION.json')['selected'];repeat=next(r for r in rows if r['state_id']=='S2_r0' and r['variant']=='DT1')
 independent=numerical(best,repeat)
 d=dict(cases=cases,comparisons=comparisons,independent_recompute=independent,
  status='PASS' if all(r['admitted'] for r in cases) and len(cases)==3 else 'FAIL',
  local_target_pass=any(r['state_id']=='S2_r0' and r['admitted'] for r in cases),
  generalization_pass=len(cases)==3 and all(r['admitted'] for r in cases),
  summary='Actual finite-amplitude search; target and tolerances were frozen before integration.',
  reason='A local target demonstration does not imply all seed configurations retain connectivity or meet the same target.',identity=identity())
 write(BASE/'analysis/TARGET.json',d);print([(r['state_id'],r['aspect_fraction'],r['final_components'],r['admitted']) for r in cases])

def extension_plan():
 rep=read(BASE/'analysis/REPLICATION.json');g=next(g for g in rep['groups'] if g['shape']=='S2' and g['action']=='A5_X180')
 assert g['full_seed_numerical_confirmation'] and g['signed_xy_direction_repeated']
 base=next(r for r in read(BASE/'discovery/INDEX.json') if r['state_id']=='S2_r0' and r['action']=='A5_X180')
 forecast=base['wall_s']*2.5*4;assert forecast<180
 write(BASE/'EXTENSION_PREREGISTRATION.json',dict(reference=base['run_id'],projected_5ms_triple_s=forecast,
  projected_raw_output_bytes=3*2.5*(BASE/base['trajectory_path']).stat().st_size if (BASE/base['trajectory_path']).exists() else 220000000,
  resource_policy='reuse audited <=2 workers; <1GB RAM per task, RTX4090 cache unchanged; require remaining physics budget >300s in runner',
  purpose='Persistence on one admitted initial realization; NOT long-time stability or new-seed generalization',identity=identity()))
 freeze('extension',[task('extension','S2_r0','A5_X180',t,v) for t in [.003,.005] for v in ['DT1','DT05','FINE']],250,
  'Repeated signed anisotropy with original dt/cache admission; bounded 3/5ms persistence on r0 only. No extrapolation to seconds.')

def extension_evaluate():
 rows=read(BASE/'extension/INDEX.json');checks=[]
 for a in rows:
  if a.get('variant')!='DT1':continue
  for b in rows:
   if b.get('variant') in ['DT05','FINE'] and b['task']['duration_s']==a['task']['duration_s']:
    checks.append(numerical(a,b) if a['status']==b['status']=='PASS' else dict(reference=a['run_id'],comparison=b['run_id'],original=dict(admitted=False)))
 write(BASE/'analysis/EXTENSION_NUMERICAL.json',checks)
