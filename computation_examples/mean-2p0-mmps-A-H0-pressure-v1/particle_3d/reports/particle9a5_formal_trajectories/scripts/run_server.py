"""One-deployment benchmark and sequential formal production controller."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,json,math,multiprocessing,os,sys,time,threading,traceback
import numpy as np
import psutil
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *
import particle_3d.formal_dynamics_p9a5 as dynamics
from particle_3d.population_inlet_p9a4 import load_new_environment,make_source
R=ROOT/REL


def monitor(stop,rows):
 process=psutil.Process()
 while not stop.is_set():
  processes=[process]+process.children(recursive=True);rss=pss=0
  for p in processes:
   try:
    rss+=p.memory_info().rss;pss+=getattr(p.memory_full_info(),'pss',p.memory_info().rss)
   except (psutil.NoSuchProcess,psutil.AccessDenied):pass
  rows.append(dict(unix_s=time.time(),rss_bytes=rss,pss_bytes=pss,load_average=os.getloadavg()))
  stop.wait(1.)


def pool_run(events,folder,workers):
 folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);start=time.perf_counter();rows=[];errors=[]
 measurements=[];stop=threading.Event()
 # Create fork pool before the monitoring thread; no background thread at fork.
 with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork')) as pool:
  futures={pool.submit(dynamics.job,(e,str(folder/'tracks'/f"mb_{e['particle_id']:06d}"))):e for e in events}
  thread=threading.Thread(target=monitor,args=(stop,measurements),daemon=True);thread.start()
  for f in as_completed(futures):
   e=futures[f]
   try:row=f.result();rows.append(row);print(json.dumps(dict(progress=len(rows),total=len(events),particle_id=e['particle_id'],status=row['status'],outlet=row['outlet'])),flush=True)
   except Exception:
    errors.append(dict(particle_id=e['particle_id'],error=traceback.format_exc()));print('EXECUTION_ERROR',e['particle_id'],flush=True)
  stop.set();thread.join()
 elapsed=time.perf_counter()-start
 if errors:
  write_new(folder/'execution_errors.json',errors)
  raise RuntimeError('Unfiltered trajectory execution errors; complete/partial records retained')
 rows=merge_rows(rows,[e['particle_id'] for e in events]);receipts=[json.loads((folder/'tracks'/f"mb_{e['particle_id']:06d}"/'receipt.json').read_text()) for e in events]
 cpu=sum(r['cpu_seconds'] for r in receipts);hashes=[r['scientific_result_sha256'] for r in receipts]
 result=dict(count=len(events),workers=workers,wall_seconds=elapsed,tracks_per_second=len(events)/elapsed,
             cpu_seconds=cpu,cpu_utilization_core_equivalents=cpu/elapsed,peak_tree_rss_bytes=max((r['rss_bytes'] for r in measurements),default=0),
             peak_tree_pss_bytes=max((r['pss_bytes'] for r in measurements),default=0),load_average_samples=measurements,
             solver_failures=sum(r['status']=='SOLVER_FAILURE' for r in rows),rows=rows,
             scientific_results_sha256=content_sha(hashes),per_particle_scientific_sha256=hashes)
 write_new(folder/'batch_result.json',result)
 return result


def init():
 if not str(ROOT).startswith('/workspace/particle9a5_formal_trajectories_'):raise ValueError('Isolated P9-A.5 server directory required')
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
  if os.environ.get(k)!='1':raise ValueError('Scientific library thread count must equal one: '+k)
 deployment=json.loads((ROOT/'deployment_manifest.json').read_text())
 for name,expected in deployment.items():
  if digest(ROOT/name)!=expected:raise ValueError('Deployment SHA mismatch: '+name)
 env=load_new_environment(ROOT);source=make_source(env,json.loads((ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json').read_text()))
 core=json.loads((R/'data/CORE500_COHORT.json').read_text());population=load_births(ROOT)
 if source.identity!=population['identity'] or core['events']!=population['events'][:500]:raise ValueError('Authoritative population mismatch')
 identity=dict(flow_sha256=FLOW_SHA,source_births_sha256=BIRTHS_SHA,core_cohort_sha256=digest(R/'data/CORE500_COHORT.json'),
  formal_contract_sha256=digest(ROOT/'particle_3d/contracts/P9A5_FORMAL_PRODUCTION_V1.json'),
  scientific_source_snapshot=source.identity,
  adapter_sha256=digest(ROOT/'particle_3d/src/particle_3d/formal_dynamics_p9a5.py'),
  orchestration_sha256=digest(ROOT/'particle_3d/src/particle_3d/formal_cohort_p9a5.py'))
 dynamics.ENV=env;dynamics.IDENTITY=identity
 return population,core,identity


def counts(rows):
 return dict(N=len(rows),**{o:sum(r['status']=='COMPLETED' and r['outlet']==o for r in rows) for o in ['O1','O2','O3']},
             completed=sum(r['status']=='COMPLETED' for r in rows),stationary=sum(r['status']=='SUPPORTED_STATIONARY' for r in rows),
             censored=sum(r['status']=='LONG_RESIDENCE_CENSORED' for r in rows),failed=sum(r['status']=='SOLVER_FAILURE' for r in rows))


def reference_dt_job(event):
 """Original runner at the explicitly requested dt, independent of formal scheduling."""
 import importlib.util
 import particle_3d.particle82a_integration as legacy
 from particle_3d.particle6_stepper import bind_query_dependency
 spec=importlib.util.spec_from_file_location('p9a5_same_dt_reference',ROOT/'particle_3d/reports/network_derived_flow_mb_validation_v1/scripts/runner.py')
 runner=importlib.util.module_from_spec(spec);sys.modules[spec.name]=runner;spec.loader.exec_module(runner)
 sample=dynamics.ENV.field.sample(event['birth_center_m'])
 runner.ENV=dynamics.ENV;runner.ROOT=ROOT;runner.IDENTITY={'role':'SAME_DT_FROZEN_KERNEL_REFERENCE','dt_s':DT}
 runner.INITIAL={str(event['particle_id']):dict(position_m=event['birth_center_m'],q=event['q'],velocity_m_s=sample.velocity_m_s.tolist(),omega_s_inv=(.5*sample.vorticity_s_inv).tolist())}
 original=legacy.integrate_admitted
 reference=bind_query_dependency(original,dict(HORIZON=HORIZONS[-1],MAX_PROVIDER_CALLS=round(16000*HORIZONS[-1]/1.5)))
 reference.__defaults__=(DT,None,False);legacy.integrate_admitted=reference
 try:return runner.mb_job((event,str(R/'outputs/same_dt_reference')))
 finally:legacy.integrate_admitted=original


def benchmark():
 population,core,identity=init();events=core['events'][:24];results=[]
 for workers in [1,2,4,6,8,12,16]:
  folder=R/'outputs/benchmark'/f'workers_{workers:02d}'
  if (folder/'batch_result.json').exists():
   result=json.loads((folder/'batch_result.json').read_text())
   if not all(completion_matches(folder/'tracks'/f"mb_{e['particle_id']:06d}",identity,e) for e in events):raise ValueError('Benchmark resume mismatch')
  else:result=pool_run(events,folder,workers)
  if results and result['scientific_results_sha256']!=results[0]['scientific_results_sha256']:raise ValueError('Worker-dependent science result')
  if result['solver_failures']:raise ValueError('Benchmark contains solver failures')
  safety_keys=['penetration_count','handoff_violation_count','inlet_escape_count','nan_inf_count','unclassified_corruption_count']
  if any(row[k] for row in result['rows'] for k in safety_keys):raise ValueError('Benchmark numerical safety gate failed')
  if workers==1:
   with ProcessPoolExecutor(max_workers=8,mp_context=multiprocessing.get_context('fork')) as pool:
    list(pool.map(reference_dt_job,events))
  for e in events:
   reference=R/'outputs/same_dt_reference/trajectories'/f"mb_{e['particle_id']:06d}.npz"
   new=folder/'tracks'/f"mb_{e['particle_id']:06d}"/'trajectory.npz'
   if not np.array_equal(np.load(reference)['samples'],np.load(new)['samples']):raise ValueError('Same-dt original-kernel reference parity failed')
  results.append(result);print('WORKER_BENCHMARK_COMPLETE',workers,result['wall_seconds'],result['tracks_per_second'],flush=True)
 resources=json.loads((ROOT/'logs/server_resources.json').read_text());q,p=map(int,resources['cgroups']['/sys/fs/cgroup/cpu.max'].split());quota=q/p
 memory=int(resources['cgroups']['/sys/fs/cgroup/memory.max']);eligible=[r for r in results if r['workers']<=math.ceil(quota) and r['peak_tree_pss_bytes']<.7*memory]
 best=max(eligible,key=lambda r:r['tracks_per_second'])
 # Translate observed per-step CPU/bytes to a conservative 12 s all-track envelope.
 sample=results[0]['rows'];trackdir=R/'outputs/benchmark/workers_01/tracks'
 costs=[]
 for row in sample:
  d=trackdir/f"mb_{row['particle_id']:06d}";rec=json.loads((d/'receipt.json').read_text());size=sum(p.stat().st_size for p in d.rglob('*') if p.is_file())
  costs.append(dict(particle_id=row['particle_id'],bytes=size,cpu_seconds=rec['cpu_seconds'],accepted_steps=row['accepted_steps'],provider_calls=row['provider_calls']))
 worst_bytes_per_step=max(r['bytes']/max(r['accepted_steps'],1) for r in costs)
 worst_cpu_per_call=max(r['cpu_seconds']/max(r['provider_calls'],1) for r in costs)
 worst_5000_bytes=math.ceil(worst_bytes_per_step*(round(HORIZONS[-1]/DT)+1)*5000*1.25)
 worst_5000_seconds=worst_cpu_per_call*128000*5000/min(best['workers'],quota)*1.25
 disk=__import__('shutil').disk_usage(ROOT);safe_disk=disk.free*.8
 # Declare a hard bound before formal production. Existing sample membership is never reduced.
 limit=5000 if worst_5000_bytes<=safe_disk else max(500,min(5000,math.floor(safe_disk/(worst_5000_bytes/5000)/250)*250))
 summary=dict(identity=identity,workers_production=best['workers'],CPU_quota_cores=quota,
  GPU_TRAJECTORY_KERNEL_AVAILABLE=False,benchmark_cohort_count=24,benchmark_science_identical=True,exact_same_dt_reference_24_trajectory_parity=True,dt_s=DT,numerical_safety_pass=True,
  selection='Maximum measured throughput among worker counts <= ceil(CPU quota), PSS <70% memory; 12/16 measured as oversubscription controls',
  configurations=[{k:v for k,v in r.items() if k not in ['rows','load_average_samples']} for r in results],
  resource_limited_N_max=limit,preferred_N_max=5000,hard_limit_reason=None if limit==5000 else 'CONSERVATIVE_ALL_12S_TRAJECTORIES_DISK_ENVELOPE',
  current_free_disk_bytes=disk.free,worst_case_5000_storage_bytes=worst_5000_bytes,worst_case_5000_runtime_seconds=worst_5000_seconds,
  cost_method=f'Maximum measured bytes per accepted step x{round(HORIZONS[-1]/DT)+1}x5000x1.25; maximum measured CPU per provider call x128000x5000/usable cores x1.25; conservative extrapolations, not strict analytical bounds',
  expected_5000_runtime_seconds=5000/best['tracks_per_second'],benchmark_cost_rows=costs)
 write_new(R/'data/worker_scaling_results.json',summary)
 print(json.dumps(summary,indent=2),flush=True)


def production():
 start=time.perf_counter();population,core,identity=init()
 gate=json.loads((R/'data/production_gate.json').read_text())
 if not gate['new_tests_pass'] or not gate['baseline_tests_pass'] or gate['identity']!=identity:raise ValueError('Production test gate required')
 bench=json.loads((R/'data/worker_scaling_results.json').read_text());plan=json.loads((R/'data/RESOURCE_PLAN.json').read_text());workers=bench['workers_production'];maximum=plan['resource_limited_N_max']
 if digest(R/'data/RESOURCE_PLAN.json')!=gate['resource_policy_sha256']:raise ValueError('Preproduction resource plan mismatch')
 if bench['identity']!=identity:raise ValueError('Benchmark source identity changed')
 previous=None;rows=[];checkpoints=[];n=500;offset=0
 while True:
  current=cohort(population,n,core['master_seed'])
  if previous:require_prefix(previous,current)
  batch=R/'outputs/formal'/('core500' if offset==0 else f'extension_{offset+1:05d}_{n:05d}')
  cohort_sha=write_new(R/'data/cohorts'/f'FORMAL_{n:05d}_COHORT.json',current)
  write_new(batch/'BATCH_MANIFEST.json',dict(start_id=offset+1,end_id=n,count=n-offset,events=current['events'][offset:n],
             formal_cohort_sha256=cohort_sha,core_cohort_sha256=digest(R/'data/CORE500_COHORT.json'),identity=identity))
  if (batch/'batch_result.json').exists():
   result=json.loads((batch/'batch_result.json').read_text())
   if not all(completion_matches(batch/'tracks'/f"mb_{e['particle_id']:06d}",identity,e) for e in current['events'][offset:n]):raise ValueError('Completed batch hash mismatch')
  else:result=pool_run(current['events'][offset:n],batch,workers)
  for r in result['rows']:r['track_relative_path']=str((batch/'tracks'/f"mb_{r['particle_id']:06d}").relative_to(R))
  rows.extend(result['rows']);rows=merge_rows(rows,range(1,n+1));summary=counts(rows)
  all_observed=all(summary[o]>0 for o in ['O1','O2','O3'])
  checkpoint=dict(**summary,added=n-offset,all_outlets_naturally_observed=all_observed,cohort_sha256=cohort_sha,batch_result_sha256=digest(batch/'batch_result.json'))
  write_new(batch/'BATCH_COMPLETE.json',checkpoint);checkpoints.append(checkpoint)
  print('FORMAL_CHECKPOINT',json.dumps(checkpoint),flush=True)
  next_n=next_size(n,all_observed,maximum)
  if next_n==n:break
  previous=current;offset=n;n=next_n
 final_sha=write_new(R/'data/FINAL_FORMAL_COHORT.json',current)
 write_new(R/'data/formal_metrics.json',rows)
 write_new(R/'data/population_checkpoints.json',checkpoints)
 summary=dict(final_formal_N=n,final_formal_cohort_sha256=final_sha,core500=counts(rows[:500]),full_formal=counts(rows),
              population_checkpoints=checkpoints,all_outlets_naturally_observed=all_observed,
              stop_reason='ALL_OUTLETS_NATURALLY_OBSERVED' if all_observed else 'PREDECLARED_N_MAX_REACHED',
              production_workers=workers,controller_runtime_seconds=time.perf_counter()-start,
              batch_compute_seconds=sum(json.loads((R/r['track_relative_path']).parents[1].joinpath('batch_result.json').read_text())['wall_seconds'] for r in rows if r['particle_id'] in [1]+[c['N']-c['added']+1 for c in checkpoints[1:]]))
 write_new(R/'data/production_complete.json',summary);print('FORMAL_PRODUCTION_COMPLETE',json.dumps(summary),flush=True)


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['benchmark','production'],required=True);a=p.parse_args()
 (benchmark if a.stage=='benchmark' else production)()
