#!/usr/bin/env python3
"""Fixed 64-ID full-trajectory CPU scaling and exact scientific parity."""
from pathlib import Path
import argparse,json,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from particle_3d.particle82_batch import run_batch,environment
from particle_3d.particle82_provenance import atomic_json,require_remote,sha256


def main(args):
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    read=lambda p:json.loads(Path(p).read_text())
    host=read(args.host_provenance);require_remote(host,host['hostname'])
    full=read(args.input_ledger);ledger=dict(full,events=full['events'][:64])
    config=dict(dt_s=.00025,guard_factor=1,dataset_role='FIXED_64_NATURAL_IDS_SCALING_DIAGNOSTIC_NOT_FORMAL_DATASET')
    # Common warm environment: do not charge cold mesh loading only to workers=1.
    cold=time.perf_counter();environment();cold=time.perf_counter()-cold
    results=[];parities=[]
    memory_limit=host['cgroup'].get('memory.max','max')
    memory_limit=int(memory_limit) if memory_limit!='max' else None
    for workers in [1,2,4,8,16]:
        if workers>host['nproc']:break
        folder=out/f'workers_{workers}'
        if args.resume and (folder/'performance.json').exists() and (folder/'catalog.json').exists():
            catalog=read(folder/'catalog.json');perf=read(folder/'performance.json')
            if len(catalog['records'])!=64 or perf['git_commit']!=host['source_git_commit']:
                raise ValueError('Saved scaling run source/count mismatch')
            for row in catalog['records']:
                pid=row['particle_id'];meta=folder/'trajectories'/f'mb_{pid:06d}.json';receipt=read(folder/'receipts'/f'mb_{pid:06d}.json')
                if receipt['hostname']!=host['hostname'] or receipt['metadata_sha256']!=sha256(meta) or row['samples_sha256']!=sha256(meta.with_suffix('.npz')):
                    raise ValueError('Stale benchmark receipt')
        else:
            # Preserve incomplete attempts; a fresh complete trial is required
            # for comparable throughput, never time only the resumed tail.
            if folder.exists() and any((folder/'receipts').glob('*.json')):
                archive=out/f'workers_{workers}_interrupted_{time.time_ns()}'
                folder.rename(archive)
            catalog,perf=run_batch(ledger,folder,workers,1,64,False,1,config,host)
        if not perf['comparable_full_fresh_batch']:raise ValueError('Cannot benchmark resume lookup throughput')
        peak=max(m['cgroup_memory_current_bytes'] or 0 for m in perf['memory_samples'])
        row={k:v for k,v in perf.items() if k!='memory_samples'}
        row['peak_cgroup_memory_bytes']=peak;row['stable_memory']=memory_limit is None or peak<.90*memory_limit
        results.append(row)
        if workers!=1:
            base=out/'workers_1';reference=read(base/'catalog.json')
            ids=[]
            for r in catalog['records']:
                pid=r['particle_id'];a=np.load(base/'trajectories'/f'mb_{pid:06d}.npz')['samples'];b=np.load(folder/'trajectories'/f'mb_{pid:06d}.npz')['samples']
                same=a.shape==b.shape and a.dtype==b.dtype and a.tobytes()==b.tobytes()
                same &= r==reference['records'][pid-1]
                ids.append(dict(particle_id=pid,exact_scientific_record_and_sample_bytes=bool(same)))
            parities.append(dict(workers=workers,ids=ids,all_exact=all(r['exact_scientific_record_and_sample_bytes'] for r in ids)))
            if not parities[-1]['all_exact']:raise ValueError('Parallelism changed scientific results')
        if not row['stable_memory']:break
    chosen=max((r for r in results if r['stable_memory']),key=lambda r:r['tracks_per_hour'])
    serial=results[0]['tracks_per_hour'];speedup=chosen['tracks_per_hour']/serial
    result=dict(fixed_ids=list(range(1,65)),selection='First 64 original natural P81 IDs, fixed before scaling; all terminal types retained',
        cold_environment_load_s=cold,worker_scaling_results=results,chosen_workers=chosen['workers'],
        serial_tracks_hour=serial,parallel_tracks_hour=chosen['tracks_per_hour'],parallel_speedup=speedup,
        parallel_efficiency=speedup/chosen['workers'],selection_rule='Maximum measured throughput among completed runs below 90% cgroup memory limit; no concurrent heavy stages',
        parallel_vs_serial=parities,all_exact=all(p['all_exact'] for p in parities),
        host=host['hostname'],source_commit=host['source_git_commit'])
    atomic_json(out/'SERVER_SCALING_BENCHMARK.json',result);print(json.dumps({k:result[k] for k in ['chosen_workers','serial_tracks_hour','parallel_tracks_hour','parallel_speedup','all_exact']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input-ledger',required=True);p.add_argument('--output-dir',required=True)
    p.add_argument('--host-provenance',required=True);p.add_argument('--resume',action='store_true');main(p.parse_args())
