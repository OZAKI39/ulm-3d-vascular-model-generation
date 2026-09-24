#!/usr/bin/env python3
"""Isolated server execution of audit calculations; protected inputs are read only."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,time,socket
from particle_3d.routing_stationary_audit import *
from particle_3d.particle8_replay import REPO

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['p65','point']);p.add_argument('--workers',type=int,default=6);a=p.parse_args()
    assert str(REPO).startswith('/workspace/particle9a1_routing_stationary_audit_')
    R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';D=R/'data';D.mkdir(parents=True,exist_ok=True)
    ref=REPO/'reference';ledger=read(ref/'birth_ledger.json');events=ledger['events'];candidates=read(ref/'all_scheduled_admission_records.json')
    dump(D/'event_identity.json',verify_events(candidates,events))
    original=read(ref/'production_source_manifest.json')
    for name,digest in original.items():assert sha(REPO/'particle_3d/src/particle_3d'/name)==digest,name
    flow=REPO/'formal_3D_flow_solver/FEM_SimVascular/frozen_reference/flow/steady_flow_mean_2p0_mmps.vtu';assert sha(flow)==FLOW_SHA
    identity=dict(flow_sha256=FLOW_SHA,production_source_sha256=original,audit_source_sha256=sha(Path(__file__).parents[1]/'src/particle_3d/routing_stationary_audit.py'),runner_sha256=sha(__file__),dt_s=DT,horizon_s=HORIZON)
    start=time.time();config=dict(stage=a.stage,workers=a.workers,host=socket.gethostname(),root=str(REPO),started_unix_s=start,identity=identity,blas_threads=1)
    dump(D/(a.stage+'_config.json'),config);environment()
    if a.stage=='p65':tasks=[(e,str(R/'diagnostic_outputs/P65_NEW'),identity) for e in events];fn=p65_job
    else:
        from particle_3d.particle82_point_native import NativePointTracer
        environment().native_point=NativePointTracer(environment())
        tasks=[(candidates[i:i+32],str(R/'diagnostic_outputs/point'),dict(step_m=.2e-6,error=1e-11,horizon_m=2e-3)) for i in range(0,len(candidates),32)];fn=point_job
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(fn,t) for t in tasks]):
            v=f.result();results.extend(v if isinstance(v,list) else [v]);print(a.stage,len(results),'/',len(events) if a.stage=='p65' else len(candidates),flush=True)
            dump(D/(a.stage+'_progress.json'),results)
    results.sort(key=lambda r:r['particle_id'])
    dump(D/(a.stage+'_completed.json'),dict(config,results=results,wall_seconds=time.time()-start))
    csvwrite(D/('p65_new_500_outcomes.csv' if a.stage=='p65' else 'all_candidate_point_routing.csv'),results)
    if a.stage=='point':csvwrite(D/'accepted500_point_routing.csv',[r for r in results if r['accepted']])
if __name__=='__main__':main()
