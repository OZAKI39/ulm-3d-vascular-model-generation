#!/usr/bin/env python3
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import time,socket,argparse
from particle_3d.routing_stationary_audit import read,sha
from particle_3d.stationary_radius_audit import radius_job,environment,dump
from particle_3d.particle8_replay import REPO

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=6);a=p.parse_args()
    assert str(REPO).startswith('/workspace/particle9a1_routing_stationary_audit_')
    R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';old=Path('/workspace/particle9a1_2mmps_20260924T001605Z/particle_3d/outputs/particle9a1_2mmps/P9A1/trajectories')
    metas=[p for p in sorted(old.glob('*.json')) if not read(p)['completed']]
    assert all('STATIONARY' in read(p)['failure_detail'] for p in metas)
    start=time.time();environment();results=[]
    config=dict(workers=a.workers,host=socket.gethostname(),started_unix_s=start,diagnostic_source_sha256=sha(REPO/'particle_3d/src/particle_3d/stationary_radius_audit.py'),reference_hashes={str(p):sha(p) for p in metas})
    dump(R/'data/radius_config.json',config)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(radius_job,(str(p),str(R/'diagnostic_outputs/virtual_radius'))) for p in metas]):
            result=f.result();results.append(result);dump(R/'data/radius_progress.json',results);print(result['particle_id'],result['classification'],[(t['radius_ratio'],t['status']) for t in result['trials']],flush=True)
    dump(R/'data/stationary_audit.json',sorted(results,key=lambda r:r['particle_id']));dump(R/'data/radius_completed.json',dict(config,wall_seconds=time.time()-start,count=len(results)))
if __name__=='__main__':main()
