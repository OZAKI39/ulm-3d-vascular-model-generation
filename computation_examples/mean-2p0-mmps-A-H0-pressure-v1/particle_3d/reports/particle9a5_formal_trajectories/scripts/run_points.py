"""Post-selection fluid-path explanation for the already frozen FINAL cohort."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,sys,os,multiprocessing,time
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *
from particle_3d.population_inlet_p9a4 import load_new_environment
from particle_3d.particle82_point_native import NativePointTracer
R=ROOT/REL;ENV=None;TRACER=None;COHORT_SHA=None

def point(e):
 folder=R/'outputs/points'/f"point_{e['particle_id']:06d}";folder.mkdir(parents=True,exist_ok=True)
 marker=folder/'COMPLETE.json'
 if marker.exists():
  m=json.loads(marker.read_text())
  if m['cohort_sha256']!=COHORT_SHA or m['event_sha256']!=content_sha(e) or any(digest(folder/n)!=h for n,h in m['files'].items()):raise ValueError('Point resume mismatch')
  return json.loads((folder/'point.json').read_text())
 result=TRACER.trace(e['birth_center_m'],step_m=.2e-6,error=1e-11,horizon_m=2e-3)
 path=result.pop('path');np.savez_compressed(folder/'path.npz',path=path)
 result.update(particle_id=e['particle_id'],source_event_id=e['source_event_id'],initial_center_m=e['birth_center_m'],
  cohort_sha256=COHORT_SHA,source_event_sha256=content_sha(e),cohort_preexists_point_execution=True,
  final_time_s=float(path[-1,0]),path_sha256=digest(folder/'path.npz'),sample_count=len(path),
  diagnostic_time_limit_s=30.,diagnostic_path_length_limit_m=2e-3,point_kernel_source_sha256=TRACER.source_sha256,
  scope='Post-selection fluid path only; point budget is the frozen 30-second diagnostic budget, not MB residence classification')
 write_new(folder/'point.json',result)
 write_new(marker,dict(cohort_sha256=COHORT_SHA,event_sha256=content_sha(e),files={n:digest(folder/n) for n in ['point.json','path.npz']}))
 return result


def main():
 global ENV,TRACER,COHORT_SHA
 if not str(ROOT).startswith('/workspace/particle9a5_formal_trajectories_'):raise ValueError('Server production directory required')
 final=R/'data/FINAL_FORMAL_COHORT.json';COHORT_SHA=digest(final);cohort=json.loads(final.read_text())
 production=json.loads((R/'data/production_complete.json').read_text())
 if production['final_formal_cohort_sha256']!=COHORT_SHA or production['final_formal_N']!=cohort['count']:raise ValueError('Final MB membership must preexist points')
 ENV=load_new_environment(ROOT);TRACER=NativePointTracer(ENV);start=time.perf_counter();rows=[]
 with ProcessPoolExecutor(max_workers=production['production_workers'],mp_context=multiprocessing.get_context('fork')) as pool:
  futures=[pool.submit(point,e) for e in cohort['events']]
  for f in as_completed(futures):
   rows.append(f.result())
   if len(rows)%50==0:print('POINT_PROGRESS',len(rows),cohort['count'],flush=True)
 rows=merge_rows(rows,range(1,cohort['count']+1))
 if digest(final)!=COHORT_SHA:raise ValueError('Point execution changed MB cohort')
 write_new(R/'data/point_metrics.json',rows)
 write_new(R/'data/point_complete.json',dict(count=len(rows),cohort_sha256=COHORT_SHA,cohort_unchanged=True,wall_seconds=time.perf_counter()-start,
  outlet_counts={(o or 'NO_EXIT'):sum(r['outlet']==o for r in rows) for o in ['OUTLET_01','OUTLET_02','OUTLET_03',None]}))
 print('POINT_EXPLANATION_COMPLETE',len(rows),flush=True)
if __name__=='__main__':main()
