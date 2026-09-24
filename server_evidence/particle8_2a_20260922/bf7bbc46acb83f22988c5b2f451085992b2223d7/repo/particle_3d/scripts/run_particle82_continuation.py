#!/usr/bin/env python3
"""Stage D: all original 720 stopped IDs, immutable accepted-prefix audit."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from particle_3d.particle82_batch import run_batch
from particle_3d.particle82_provenance import atomic_json,require_remote


def main(a):
    read=lambda p:json.loads(Path(p).read_text());host=read(a.host_provenance);require_remote(host,host['hostname'])
    previous=Path(a.previous);out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    admission=read(a.admission_audit);scale=read(a.scaling_benchmark)
    original=[read(p) for p in sorted((previous/'trajectories').glob('*.json'))]
    stopped=[r for r in original if r['end_reason']=='INTEGRATION_SAFETY_STOP']
    if len(stopped)!=720:raise ValueError('Expected all original 720 stops')
    ledger=read(previous/'data/birth_ledger.json');ledger['events']=[r['birth_metadata'] for r in stopped]
    summaries=[]
    for factor in [4,16]:
        config=dict(dt_s=.00025,guard_factor=factor,dataset_role='ORIGINAL_720_SAFETY_STOPS_EXTENDED_GUARD_DIAGNOSTIC',
            checkpoint_directory=str(previous/'trajectories'))
        catalog,perf=run_batch(ledger,out/f'guard_x{factor}',scale['chosen_workers'],1,2200,True,1,config,host)
        rows=[]
        for record in catalog['records']:
            pid=record['particle_id'];before=np.load(previous/'trajectories'/f'mb_{pid:06d}.npz')['samples']
            after=np.load(out/f'guard_x{factor}/trajectories/mb_{pid:06d}.npz')['samples']
            prefix=len(after)>=len(before) and before.tobytes()==after[:len(before)].tobytes()
            rows.append(dict(particle_id=pid,point_basin=admission['per_id'][str(pid)]['admitted_basin'],
                exact_original_prefix=bool(prefix),original_states=len(before),extended_states=len(after),
                completed=record['completed'],outlet=record['exit_outlet'],end_reason=record['end_reason'],
                failure_detail=record['failure_detail'],last_elapsed_s=float(after[-1,0]) if len(after) else None))
        summary=dict(guard_factor=factor,baseline_completed=1249,baseline_stopped=720,rows=rows,
            completed_from_stops=sum(r['completed'] for r in rows),still_stopped=sum(not r['completed'] for r in rows),
            all_original_prefixes_exact=all(r['exact_original_prefix'] for r in rows),performance=perf)
        atomic_json(out/f'EXTENDED_GUARD_X{factor}.json',summary);summaries.append(summary)
        if not summary['all_original_prefixes_exact']:raise ValueError('Continuation changed a previously accepted prefix')
    atomic_json(out/'EXTENDED_GUARD_AUDIT.json',dict(factors=summaries,
        same_frozen_physics_radius_admission_and_dt=True,no_physiological_trapping_claim=True))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['previous','output-dir','host-provenance','admission-audit','scaling-benchmark']:p.add_argument('--'+name,required=True)
    main(p.parse_args())
