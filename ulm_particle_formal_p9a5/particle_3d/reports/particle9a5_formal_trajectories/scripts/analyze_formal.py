"""Read-only CORE500/FULL statistics, sequential-look uncertainty and safety gates."""
from pathlib import Path
import json,csv,sys,xml.etree.ElementTree as ET
import numpy as np
from scipy.stats import binomtest
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *
R=ROOT/REL;D=R/'data'


def stats(values):
 a=np.asarray([v for v in values if v is not None],dtype=float)
 return dict(n=len(a),minimum=float(a.min()),p10=float(np.quantile(a,.1)),median=float(np.median(a)),mean=float(a.mean()),p90=float(np.quantile(a,.9)),maximum=float(a.max())) if len(a) else dict(n=0,minimum=None,p10=None,median=None,mean=None,p90=None,maximum=None)


def tests(path):
 r=ET.parse(path).getroot();ss=[r] if r.tag=='testsuite' else list(r)
 return {k:sum(int(s.get(k,0)) for s in ss) for k in ['tests','failures','errors','skipped']}


def group(rows,adaptive=False):
 n=len(rows);counts={o:sum(r['status']=='COMPLETED' and r['outlet']==o for r in rows) for o in ['O1','O2','O3']}
 counts.update(stationary=sum(r['status']=='SUPPORTED_STATIONARY' for r in rows),censored=sum(r['status']=='LONG_RESIDENCE_CENSORED' for r in rows),failed=sum(r['status']=='SOLVER_FAILURE' for r in rows))
 fractions={}
 for name,k in counts.items():
  # Familywise 95% coverage over 6 outcome categories x at most 13 scheduled
  # looks. This conservative exact interval remains valid at the coverage stop.
  level=1-.05/(6*13) if adaptive else .95
  ci=binomtest(k,n).proportion_ci(confidence_level=level,method='exact' if adaptive else 'wilson')
  fractions[name]=dict(count=k,fraction=k/n,lower=float(ci.low),upper=float(ci.high))
 return dict(N=n,**counts,completed=sum(counts[o] for o in ['O1','O2','O3']),fractions=fractions,
             uncertainty_method='Clopper-Pearson with Bonferroni over 6 categories and 13 predeclared looks; familywise >=95%' if adaptive else 'Fixed CORE500; per-category Wilson 95% interval',
             diameter_stats=stats([r['diameter_um'] for r in rows]),
             transit_stats=stats([r['transit_time_s'] for r in rows if r['status']=='COMPLETED']),
             wall_gap_stats=stats([r['minimum_wall_gap_m'] for r in rows]),nearwall_stats=stats([r['nearwall_exposure_s'] for r in rows]))


def main():
 cohort_path=D/'FINAL_FORMAL_COHORT.json';frozen=json.loads(cohort_path.read_text());n=frozen['count']
 core=json.loads((D/'CORE500_COHORT.json').read_text());require_prefix(core,frozen)
 population=load_births(ROOT);assert frozen['events']==population['events'][:n]
 rows=merge_rows(json.loads((D/'formal_metrics.json').read_text()),range(1,n+1))
 for e,row in zip(frozen['events'],rows):
  folder=R/row['track_relative_path'];marker=json.loads((folder/'COMPLETE.json').read_text())
  if not completion_matches(folder,marker['identity'],e):raise ValueError('Incomplete/corrupt formal trajectory: '+str(e['particle_id']))
  assert row['source_event_id']==e['source_event_id'] and row['diameter_um']==e['diameter_um']
  assert digest(folder/'trajectory.npz')==row['trajectory_file_sha256']
 points=merge_rows(json.loads((D/'point_metrics.json').read_text()),range(1,n+1))
 assert all(p['cohort_sha256']==digest(cohort_path) and p['initial_center_m']==e['birth_center_m'] for e,p in zip(frozen['events'],points))
 core_stats=group(rows[:500]);full_stats=group(rows,True);bench=json.loads((D/'worker_scaling_results.json').read_text());production=json.loads((D/'production_complete.json').read_text());context=json.loads((D/'run_context.json').read_text())
 bench.update(json.loads((D/'RESOURCE_PLAN.json').read_text()))
 first={o:next((r['particle_id'] for r in rows if r['status']=='COMPLETED' and r['outlet']==o),None) for o in ['O1','O2','O3']};coverage=all(v is not None for v in first.values())
 safety={key:sum(r[key] for r in rows) for key in ['penetration_count','handoff_violation_count','inlet_escape_count','nan_inf_count','unclassified_corruption_count']}
 baseline=tests(R/'logs/baseline_tests.xml');new=tests(R/'logs/new_tests.xml');final=tests(R/'logs/final_tests.xml') if (R/'logs/final_tests.xml').exists() else None
 failed_tests=any(t[k] for t in [baseline,new]+([final] if final else []) for k in ['failures','errors','skipped'])
 if failed_tests:status='P9A5_REGRESSION_FAIL'
 elif any(safety.values()) or full_stats['failed']:status='P9A5_FORMAL_TRAJECTORIES_NUMERICAL_SAFETY_FAIL'
 elif full_stats['censored']:status='P9A5_FORMAL_TRAJECTORIES_LONG_RESIDENCE_REVIEW_REQUIRED'
 elif full_stats['stationary']:status='P9A5_FORMAL_TRAJECTORIES_PASS_WITH_SUPPORTED_STATIONARY'
 elif coverage:status='P9A5_FORMAL_TRAJECTORIES_PASS_ALL_OUTLETS_OBSERVED'
 else:status='P9A5_FORMAL_TRAJECTORIES_PASS_OUTLET_COVERAGE_INCOMPLETE'
 differences={k:full_stats['fractions'][k]['fraction']-core_stats['fractions'][k]['fraction'] for k in core_stats['fractions']}
 transitions={};point_counts={k:0 for k in ['O1','O2','O3','NO_EXIT']}
 for row,p in zip(rows,points):
  po='O'+str(int(p['outlet'].split('_')[-1])) if p['outlet'] else 'NO_EXIT';point_counts[po]+=1
  mb=row['outlet'] if row['status']=='COMPLETED' else row['status'];key=po+' → '+mb;transitions[key]=transitions.get(key,0)+1
  row['point_outlet']=po;row['point_end_reason']=p['end_reason']
 horizons=[]
 for h in [1.5,3.,6.,12.]:
  resolved=[r for r in rows if r['trajectory_age_s']<=h+1e-12]
  horizons.append(dict(horizon_s=h,completed=sum(r['status']=='COMPLETED' for r in resolved),stationary=sum(r['status']=='SUPPORTED_STATIONARY' for r in resolved),
    failed=sum(r['status']=='SOLVER_FAILURE' for r in resolved),still_active_or_censored=n-sum(r['status'] in ['COMPLETED','SUPPORTED_STATIONARY','SOLVER_FAILURE'] for r in resolved)))
 diameter_outcomes={k:stats([r['diameter_um'] for r in rows if (r['outlet'] if r['status']=='COMPLETED' else r['status'])==k]) for k in ['O1','O2','O3','SUPPORTED_STATIONARY','LONG_RESIDENCE_CENSORED','SOLVER_FAILURE']}
 summary=dict(baseline_commit=context['baseline_commit'],new_branch=context['branch'],new_worktree=context['root'],server_root=context['remote'],
  flow_sha256=FLOW_SHA,source_contract_sha=SOURCE_CONTRACT_SHA,previous_formal_horizon_s=json.loads((D/'historical_horizon_audit.json').read_text())['previous_formal_horizon_s'],
  initial_horizon_s=3.,maximum_horizon_used_s=max(r['maximum_horizon_used_s'] for r in rows),dt_s=DT,
  horizon_extension_3_to_6_count=sum([3.,6.] in r['horizon_extensions'] for r in rows),horizon_extension_6_to_12_count=sum([6.,12.] in r['horizon_extensions'] for r in rows),
  core_N=500,core_cohort_sha=digest(D/'CORE500_COHORT.json'),final_formal_N=n,final_formal_cohort_sha=digest(cohort_path),
  increment_size_initial=250,increment_size_after_2000=500,first_O1_at_N=first['O1'],first_O2_at_N=first['O2'],first_O3_at_N=first['O3'],
  all_outlets_observed_at_N=max(first.values()) if coverage else None,
  all_outlets_observed_at_checkpoint_N=next((c['N'] for c in production['population_checkpoints'] if c['all_outlets_naturally_observed']),None),
  all_outlets_naturally_observed=coverage,minimum_500_satisfied=n>=500,core500=core_stats,full_formal=full_stats,
  core_full_fraction_difference=differences,**safety,diameter_stats=full_stats['diameter_stats'],diameter_by_outcome=diameter_outcomes,
  transit_stats=full_stats['transit_stats'],wall_gap_stats=full_stats['wall_gap_stats'],nearwall_stats=full_stats['nearwall_stats'],
  point_basin_counts=point_counts,point_mb_transitions=transitions,point_budget_note='Frozen CPU native diagnostic tracer: 30 s / 2 mm budget; not used for MB selection or MB stationary classification',
  production_workers=production['production_workers'],worker_scaling_results=bench,
  runtime_seconds=production['batch_compute_seconds'],benchmark_runtime_seconds=sum(c['wall_seconds'] for c in bench['configurations']),
  point_runtime_seconds=json.loads((D/'point_complete.json').read_text())['wall_seconds'],
  baseline_tests=baseline,new_tests=new,final_tests=final,old_files_changed=json.loads((D/'protected_verification.json').read_text())['old_files_changed'] if (D/'protected_verification.json').exists() else None,dt_revision=json.loads((D/'dt_revision.json').read_text()),
  core500_all_outlets_observed=all(core_stats[o]>0 for o in ['O1','O2','O3']),stop_reason=production['stop_reason'],
  population_checkpoints=production['population_checkpoints'],resource_limited_N_max=bench['resource_limited_N_max'],
  horizon_resolution=horizons,final_status=status,
  scope='Independent finite-size single MB trajectories in fixed NEW Network-H0 FEM field; no MB-MB/RBC coupling; supported stationary is a constrained contact equilibrium, not physiological trapping proof')
 # Analysis summaries may be refreshed after final regression; raw evidence is immutable.
 (D/'final_summary.json').write_bytes(canonical(summary))
 (D/'analysis_metrics.json').write_bytes(canonical(rows))
 with (D/'trajectory_catalog.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
  w.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
 print(json.dumps(dict(final_N=n,core500=core_stats,full_formal=full_stats,safety=safety,final_status=status),indent=2))
if __name__=='__main__':main()
