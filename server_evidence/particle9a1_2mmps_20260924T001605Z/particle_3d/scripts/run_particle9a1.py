#!/usr/bin/env python3
"""P9-A.1 same-cohort validation and strictly gated production; fresh cache."""
import argparse,gzip,hashlib,json,os,resource,socket,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from copy import deepcopy
import numpy as np
from particle_3d.particle9a_provenance import require_current_flow,sha256,REPO
from particle_3d.particle81_simulation import SavedTrajectoryStepper,environment,dump,DT,HORIZON
from particle_3d.particle9a_motion import Particle9AStepper
from particle_3d.particle9a1_audit import instrument
from particle_3d.particle82a_integration import integrate_admitted
from particle_3d.particle6_stepper import bind_query_dependency
from particle_3d.particle8_replay import canonical_hash


def identity():
 return dict(flow=require_current_flow(),source_sha256={p.name:sha256(p) for p in sorted((REPO/'particle_3d/src/particle_3d').glob('*.py'))},runner_sha256=sha256(__file__),dt_s=DT,horizon_s=HORIZON)


def job(spec):
 event,model,output,logging,provenance=spec
 start=time.time();cpu=time.process_time();env=environment();rows=[];created=[]
 folder=Path(output)/model;meta=folder/'trajectories'/f'mb_{event["particle_id"]:06d}.json'
 if meta.exists():
  previous=json.loads(meta.read_text())
  if previous.get('p9a1_identity')!=provenance or previous.get('model')!=model:raise ValueError('P9A1 cache identity mismatch')
 def factory(*args,**kwargs):
  stepper=(SavedTrajectoryStepper(*args,**kwargs) if model=='P65_NEW' else Particle9AStepper(*args,gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv,**kwargs))
  created.append(stepper);return instrument(stepper,enabled=logging,rows=rows)
 integrate=bind_query_dependency(integrate_admitted,{'SavedTrajectoryStepper':factory})
 before=canonical_hash(event);result=integrate(deepcopy(event),output=folder)
 assert canonical_hash(event)==before
 if created:
  result.update(model=model,dynamics=model,p9a1_identity=provenance,diagnostics_enabled=logging,
   planar_statistics=getattr(created[0],'planar_statistics',None),
   receipt=dict(hostname=socket.gethostname(),pid=os.getpid(),started_unix_s=start,ended_unix_s=time.time(),
    cpu_seconds=time.process_time()-cpu,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,remote_server_compute=str(REPO).startswith('/workspace/particle9a1_2mmps_')))
  dump(meta,result)
  if logging:
   with gzip.open(meta.with_suffix('.audit.jsonl.gz'),'wt') as f:
    for row in rows:f.write(json.dumps(row,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())+'\n')
 return dict(particle_id=event['particle_id'],model=model,completed=result['completed'],end_reason=result['end_reason'],failure_detail=result['failure_detail'],
  minimum_gap_m=result.get('minimum_original_wall_gap_m'),last_time_s=result.get('last_elapsed_time_s'),path_length_m=result.get('path_length_m'),outlet=result['exit_outlet'],
  accepted_steps=result.get('accepted_steps'),rejected_trials=result.get('rejected_trials'),wall_seconds=result['wall_seconds'],reused=not bool(created))


def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['same12','smoke30','production'],required=True);p.add_argument('--workers',type=int,default=6)
 p.add_argument('--ledger',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--diagnostics',action='store_true');p.add_argument('--gate',type=Path)
 a=p.parse_args()
 if not str(REPO).startswith('/workspace/particle9a1_2mmps_'):raise ValueError('Use isolated authorized server directory')
 if a.stage!='same12':
  gate=json.loads(a.gate.read_text())
  required='same12_pass' if a.stage=='smoke30' else 'production_authorized_by_all_gates'
  if gate.get(required) is not True:raise ValueError('Escalation gate not passed: '+required)
 ledger=json.loads(a.ledger.read_text());provenance=identity();events=ledger['events'][:{'same12':12,'smoke30':30,'production':500}[a.stage]]
 if a.stage=='same12':assert [e['particle_id'] for e in events]==[3,5,7,10,11,12,14,15,17,18,19,22]
 a.output.mkdir(parents=True,exist_ok=True);start=time.time()
 dump(a.output/'admission/birth_ledger.json',dict(ledger,schema='P9A1_CURRENT_2MMPS_IDENTICAL_ADMISSION_EVENTS_NEW_STAGE_LEDGER',p9a1_identity=provenance,source_ledger_sha256=sha256(a.ledger)))
 env=environment();dump(a.output/'provenance/open_rim_topology.json',dict(env.wall.open_boundary_topology['provenance'],edges=[dict(global_node_ids_zero_based=k,roles=v) for k,v in sorted(env.wall.open_boundary_topology['edges'].items())]))
 config=dict(stage=a.stage,workers=a.workers,diagnostics=a.diagnostics,server_root=str(REPO),hostname=socket.gethostname(),started_unix_s=start,p9a1_identity=provenance,event_sha256={e['particle_id']:canonical_hash(e) for e in events},Q_in_m3_s=env.sampler.Q_m3_s)
 dump(a.output/'provenance'/f'{a.stage}_config.json',config)
 models=['P9A1','P65_NEW'] if a.stage=='same12' else ['P9A1'];results=[]
 with ProcessPoolExecutor(max_workers=a.workers) as pool:
  futures=[pool.submit(job,(e,m,str(a.output),a.diagnostics,provenance)) for m in models for e in events]
  for f in as_completed(futures):
   r=f.result();results.append(r);print(json.dumps(r),flush=True)
   dump(a.output/'provenance'/f'{a.stage}_progress.json',results)
 dump(a.output/'provenance'/f'{a.stage}_completed.json',dict(config,results=results,wall_seconds=time.time()-start))
if __name__=='__main__':main()
