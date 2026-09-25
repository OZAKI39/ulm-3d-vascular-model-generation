"""Remote-only, atomic, resume-safe independent trajectory multiprocessing.

Linux fork shares the preloaded immutable Frozen mesh/flow via copy-on-write.
Scientific merge order is stable-ID order; execution receipts intentionally
remain separate because PID, timing and worker assignment are not deterministic.
"""
from pathlib import Path
from datetime import datetime,timezone
import concurrent.futures
import json,multiprocessing,os,socket,time
import numpy as np
from .particle82_provenance import atomic_json,sha256,require_remote
from .particle8_replay import canonical_hash
from .particle82_integration import integrate_one,environment


def utc():return datetime.now(timezone.utc).isoformat()


def sorted_unique(records):
    records=sorted(records,key=lambda r:r['particle_id'])
    ids=[r['particle_id'] for r in records]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate stable IDs in merged catalog')
    return records


def physics_record(meta):
    fields=['particle_id','birth_metadata','birth_metadata_sha256','birth_time_s','radius_m','diameter_um',
        'initial_q','admission_candidates','integration_config','end_reason','failure_detail','completed',
        'exit_outlet','exit_time_s','residence_time_s','path_length_m','sample_count','samples_sha256',
        'inlet_position_m','inlet_face','inlet_barycentric','provider_calls','accepted_steps','rejected_trials']
    return {k:meta[k] for k in fields if k in meta}


def task_run(task):
    events,folder,config,host,resume=task
    require_remote(host,host['hostname'])
    folder=Path(folder);result=[]
    for event in events:
        pid=event['particle_id'];receipt_path=folder/'receipts'/f'mb_{pid:06d}.json'
        extra={} if not config.get('checkpoint_directory') else dict(checkpoint=Path(config['checkpoint_directory'])/f'mb_{pid:06d}.json')
        mp=folder/'trajectories'/f'mb_{pid:06d}.json'
        if mp.exists() and not resume:raise ValueError('Existing ID requires --resume')
        if receipt_path.exists():
            receipt=json.loads(receipt_path.read_text())
            if receipt['config_sha256']!=canonical_hash(config) or receipt['git_commit']!=host['source_git_commit']:
                raise ValueError('Resume configuration/source mismatch')
            if receipt['hostname']!=host['hostname'] or receipt['metadata_sha256']!=sha256(mp):
                raise ValueError('Resume host/metadata mismatch')
            if receipt['frozen_input_sha256']!=host['frozen_input_sha256']:
                raise ValueError('Resume Frozen input mismatch')
            # Native resume also checks original birth and accepted sample digest.
            meta=integrate_one(event,dt=config['dt_s'],guard_factor=config['guard_factor'],output=folder,**extra)
            result.append(physics_record(meta));continue
        start=utc();t=time.perf_counter()
        # A crash between atomic NPZ/metadata and receipt is recoverable: rerun
        # the deterministic ID if its provenance transaction was not committed.
        meta=integrate_one(event,dt=config['dt_s'],guard_factor=config['guard_factor'],output=folder,force=True,**extra)
        receipt=dict(hostname=socket.gethostname(),pid=os.getpid(),worker_id=multiprocessing.current_process().name,
            git_commit=host['source_git_commit'],config_sha256=canonical_hash(config),
            frozen_input_sha256=host['frozen_input_sha256'],start_time=start,end_time=utc(),
            wall_seconds=time.perf_counter()-t,particle_id=pid,metadata_sha256=sha256(mp),
            sample_sha256=meta['samples_sha256'],role='REMOTE_SERVER_HOST',
            ssh_endpoint=host['ssh_endpoint'],host_provenance_sha256=canonical_hash(host))
        atomic_json(receipt_path,receipt);result.append(physics_record(meta))
    ids=[e['particle_id'] for e in events]
    atomic_json(folder/'shards'/f'shard_{ids[0]:06d}_{ids[-1]:06d}.json',dict(
        hostname=socket.gethostname(),pid=os.getpid(),worker_id=multiprocessing.current_process().name,
        git_commit=host['source_git_commit'],config_sha256=canonical_hash(config),
        frozen_input_sha256=host['frozen_input_sha256'],particle_ids=ids,
        individual_receipts=[f'receipts/mb_{i:06d}.json' for i in ids],end_time=utc()))
    return result


def run_batch(ledger,output,workers,id_start,id_end,resume,shard_size,config,host):
    import psutil
    require_remote(host,host['hostname'])
    if not np.isfinite(config['dt_s']) or config['dt_s']<=0:raise ValueError('Positive finite dt required')
    if config['guard_factor'] not in [1,4,16]:raise ValueError('Unregistered guard factor')
    if not 1<=id_start<=id_end:raise ValueError('Invalid stable-ID interval')
    if workers<1 or workers>len(os.sched_getaffinity(0)):raise ValueError('Workers exceed available CPUs')
    if shard_size<1:raise ValueError('Positive shard size required')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    all_events=ledger['events'];events=[e for e in all_events if id_start<=e['particle_id']<=id_end]
    if not events:raise ValueError('No IDs selected')
    if len({e['particle_id'] for e in events})!=len(events):raise ValueError('Duplicate birth IDs')
    source=Path(__file__).resolve().parents[3]/'SOURCE_COMMIT'
    if not source.exists() or source.read_text().strip()!=host['source_git_commit']:
        raise ValueError('Deployed source snapshot commit mismatch')
    for relative,digest in host['frozen_input_sha256'].items():
        if sha256(source.parent/relative)!=digest:raise ValueError('Frozen input changed: '+relative)
    for relative,digest in host['source_sha256'].items():
        if sha256(source.parent/relative)!=digest:raise ValueError('Deployed source changed: '+relative)
    for path,digest in host.get('external_scientific_input_sha256',{}).items():
        if sha256(path)!=digest:raise ValueError('External scientific input changed: '+path)
    existing=output/'RUN_CONFIG.json'
    manifest=dict(config=config,ledger_sha256=canonical_hash(ledger),source_commit=host['source_git_commit'],
                  id_start=id_start,id_end=id_end,formal_role=config['dataset_role'])
    if existing.exists() and json.loads(existing.read_text())!=manifest:raise ValueError('Resume run contract mismatch')
    atomic_json(existing,manifest)
    reused_ids=[e['particle_id'] for e in events if (output/'receipts'/f'mb_{e["particle_id"]:06d}.json').exists()]
    start=time.perf_counter();started=utc();env=environment()
    # TetraGeometry, nodal field arrays and wall triangles are read-only. VTK
    # locator caches are private process COW pages; peak PSS/RSS is measured.
    readonly=all(not a.flags.writeable for a in [env.field.points_m,env.field.tetra,env.field.velocity_nodes_m_s])
    if not readonly:raise ValueError('Frozen field must be read-only before fork')
    tasks=[(events[i:i+shard_size],str(output),config,host,resume) for i in range(0,len(events),shard_size)]
    records=[];metrics=[];proc=psutil.Process();psutil.cpu_percent()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork')) as pool:
        pending={pool.submit(task_run,t) for t in tasks}
        while pending:
            done,pending=concurrent.futures.wait(pending,timeout=1,return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:records.extend(future.result())
            processes=[proc,*proc.children(recursive=True)];rss=pss=cpu_s=0.
            for child in processes:
                try:
                    m=child.memory_full_info();rss+=m.rss;pss+=getattr(m,'pss',0)
                    ct=child.cpu_times();cpu_s+=ct.user+ct.system
                except psutil.NoSuchProcess:pass
            vm=psutil.virtual_memory();swap=psutil.swap_memory()
            metrics.append(dict(elapsed_s=time.perf_counter()-start,process_rss_sum_bytes=rss,process_pss_sum_bytes=pss,
                process_cpu_seconds=cpu_s,system_cpu_percent=psutil.cpu_percent(),load_average=list(os.getloadavg()),
                system_available_bytes=vm.available,system_total_bytes=vm.total,swap_used_bytes=swap.used,
                cgroup_memory_current_bytes=int(Path('/sys/fs/cgroup/memory.current').read_text()) if Path('/sys/fs/cgroup/memory.current').exists() else None))
            if done:print('MB_BATCH',len(records),'/',len(events),'workers',workers,flush=True)
    records=sorted_unique(records);elapsed=time.perf_counter()-start
    catalog=dict(schema='PARTICLE82_SCIENTIFIC_CATALOG_V1',dataset_role=config['dataset_role'],records=records,
        scientific_records_sha256=canonical_hash(records),scheduled=len(records),admitted=sum(bool(r['sample_count']) for r in records),
        completed=sum(r['completed'] for r in records),outlet_counts={o:sum(r['exit_outlet']==o for r in records) for o in ['OUTLET_01','OUTLET_02','OUTLET_03']},
        config=config,source_commit=host['source_git_commit'],host_provenance=host,server_only_formal_dataset=True)
    atomic_json(output/'catalog.json',catalog)
    computed=len(records)-len(reused_ids)
    summary=dict(hostname=socket.gethostname(),workers=workers,tracks=len(records),computed_tracks=computed,
        resumed_tracks=len(reused_ids),wall_seconds=elapsed,comparable_full_fresh_batch=not reused_ids,
        tracks_per_second=computed/elapsed,tracks_per_hour=computed/elapsed*3600,start_time=started,end_time=utc(),
        peak_process_rss_sum_bytes=max(m['process_rss_sum_bytes'] for m in metrics),
        peak_process_pss_sum_bytes=max(m['process_pss_sum_bytes'] for m in metrics),
        sharing='LINUX_FORK_COPY_ON_WRITE_READONLY_FROZEN_ARRAYS',memory_samples=metrics,
        config_sha256=canonical_hash(config),git_commit=host['source_git_commit'])
    atomic_json(output/'performance.json',summary)
    return catalog,summary
