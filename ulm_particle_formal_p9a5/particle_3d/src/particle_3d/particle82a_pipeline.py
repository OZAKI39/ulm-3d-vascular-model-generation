"""Remote-only sharded admission audit, deterministic scheduling and receipts."""
from pathlib import Path
from collections import Counter
import argparse, gzip, hashlib, json, multiprocessing as mp, os, resource, socket, time
import numpy as np
from .particle82a_admission import context, common_event, audit_event, MASTER_SEED
from .particle82a_geometry import FAILURES, maximum_handoff_radius
from .particle82_provenance import atomic_json, sha256, require_remote
from .particle8_replay import canonical_hash

_RUN = None
CAUSES = ['ACCEPTED', *FAILURES]
BASINS = ['OUTLET_01', 'OUTLET_02', 'OUTLET_03', 'UNRESOLVED_POINT_PATH']


def atomic_npz(target, **arrays):
    path = Path(target); temp = path.with_name(path.stem+f'.{os.getpid()}.tmp.npz')
    np.savez_compressed(temp, **arrays); os.replace(temp, path)


def atomic_gzip(path, value):
    path = Path(path); temp = path.with_name(path.name+f'.{os.getpid()}.tmp')
    with gzip.open(temp, 'wt') as stream:
        json.dump(value, stream, ensure_ascii=False, allow_nan=False)
    os.replace(temp, path)


def initialize(output, provenance, config):
    global _RUN
    require_remote(provenance, provenance['hostname'])
    _RUN = dict(output=str(output), provenance=provenance, config=config,
                config_hash=canonical_hash(config))
    context()  # Immutable FEM/VTK/BVH inherited by Linux fork, no spawn copies.


def receipt(kind, first, count, start, **extra):
    p = _RUN['provenance']
    return dict(kind=kind, first_id=first, count=count, hostname=socket.gethostname(),
        pid=os.getpid(), git_commit=p['source_git_commit'], config_hash=_RUN['config_hash'],
        frozen_input_hash=canonical_hash(p['frozen_input_sha256']), start_time=start, end_time=time.time(),
        REMOTE_SERVER_COMPUTE=True, worker=mp.current_process().name,
        peak_worker_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, **extra)


def valid_resume(meta_path):
    if not meta_path.exists():
        return None
    meta = json.loads(meta_path.read_text())
    if meta['config_hash'] != _RUN['config_hash'] or meta['git_commit'] != _RUN['provenance']['source_git_commit']:
        raise ValueError('Resume config/source mismatch: '+str(meta_path))
    if meta['hostname'] != socket.gethostname() or not meta['REMOTE_SERVER_COMPUTE']:
        raise ValueError('Resume compute-host mismatch')
    for file, digest in meta['files'].items():
        if sha256(meta_path.parent/file) != digest:
            raise ValueError('Corrupt shard '+file)
    return meta


def event_shard(task):
    first, count, subdir = task
    folder = Path(_RUN['output'])/subdir; folder.mkdir(parents=True, exist_ok=True)
    stem = f'events_{first:07d}_{first+count-1:07d}'; meta_path = folder/(stem+'.receipt.json')
    old = valid_resume(meta_path)
    if old is not None:
        return old
    start = time.time(); rows = []; attempts = []; paths = []; offsets = [0]; c = context()
    original_parity = []; accepted = Counter()
    from .injection_admission import FiniteSizeAdmission
    original = FiniteSizeAdmission(None, None, wall=c.env.wall, field=c.env.field)
    for event_id in range(first, first+count):
        e = common_event(event_id)
        e, methods, path = audit_event(e, point_step=_RUN['config']['point_step_m'],
            search_tolerance=_RUN['config']['search_tolerance_m'],
            horizon_factor=_RUN['config'].get('horizon_factor', 1.))
        for mi, result in enumerate(methods):
            accepted[result['method']] += int(result['accepted'])
            candidate_rows = result.pop('attempts')
            result['failure_counts'] = dict(Counter(x['cause'] for x in candidate_rows if x['cause'] != 'ACCEPTED'))
            for attempt in candidate_rows:
                attempts.append([mi, event_id, attempt['draw'], attempt['diameter_um'],
                    *attempt['position_m'], attempt['triangle'], CAUSES.index(attempt['cause']),
                    attempt.get('wall_distance_m', -1), attempt.get('perimeter_distance_m', -1)])
            # Independent original-check parity for all accepted A/B and all
            # first candidates plus deterministic rejection samples, not just easy passes.
            if mi < 2:
                indices = sorted(set([0, len(candidate_rows)-1]+list(range(0, len(candidate_rows), 64))))
                for k in indices:
                    row = candidate_rows[k]
                    p, status, detail = original.check(dict(species='MB', particle_id=event_id,
                        radius_m=row['diameter_um']*.5e-6, q=[1, 0, 0, 0]), np.asarray(row['position_m']), {})
                    equal = (p is not None) == (row['cause'] == 'ACCEPTED')
                    original_parity.append(dict(event_id=event_id, method=result['method'], draw=k,
                        original_status=status, classified_cause=row['cause'], equal=equal))
                    if not equal:
                        raise ValueError('Original admission checker parity failed: '+str(original_parity[-1]))
        rows.append(dict(event=e, methods=dict(zip('ABC', methods))))
        paths.append(path); offsets.append(offsets[-1]+len(path))
    raw_path = folder/(stem+'.npz'); rows_path = folder/(stem+'.json.gz')
    atomic_npz(raw_path, attempts=np.asarray(attempts, float), paths=np.concatenate(paths), offsets=np.asarray(offsets))
    atomic_gzip(rows_path, dict(rows=rows, original_checker_parity=original_parity,
        attempt_columns=['method_index','event_id','draw','diameter_um','x_m','y_m','z_m','triangle','cause_index','wall_distance_m','perimeter_distance_m'],
        cause_labels=CAUSES))
    meta = receipt('COMMON_LEDGER_PAIRED_ABC', first, count, start,
        accepted=dict(accepted), original_checker_parity_count=len(original_parity), original_checker_parity=True,
        files={p.name:sha256(p) for p in [raw_path, rows_path]})
    atomic_json(meta_path, meta)
    return meta


def map_shard(task):
    first, points, faces, weights, subdir = task
    folder = Path(_RUN['output'])/subdir; folder.mkdir(parents=True, exist_ok=True)
    stem = f'map_{first:07d}_{first+len(points)-1:07d}'; meta_path = folder/(stem+'.receipt.json')
    old = valid_resume(meta_path)
    if old is not None:
        return old
    start = time.time(); c = context(); rows = []
    for point, face, weight in zip(points, faces, weights):
        trace = c.native.trace(point, step_m=_RUN['config']['point_step_m'])
        wall, perimeter, _ = c.geometry.distances(point)
        radius = min(wall, perimeter)
        rows.append([*point, face, weight, BASINS.index(trace['outlet'] or 'UNRESOLVED_POINT_PATH'),
            wall, perimeter, 2*radius*1e6, c.distribution.cdf(2*radius*1e6),
            c.distribution.cdf(2*maximum_handoff_radius(wall)*1e6)])
    path = folder/(stem+'.npz'); atomic_npz(path, rows=np.asarray(rows))
    meta = receipt('INLET_QUADRATURE_MAP', first, len(points), start, files={path.name:sha256(path)})
    atomic_json(meta_path, meta); return meta


def triangle_quadrature(n, c=None):
    c = c or context(); bary = []
    # Equal-area n^2 subtriangles; centroids integrate linear flux exactly.
    for i in range(n):
        for j in range(n-i):
            vertices = np.array([[i,j],[i+1,j],[i,j+1]], float)/n
            uv = vertices.mean(axis=0); bary.append([1-uv.sum(), *uv])
            if i+j < n-1:
                vertices = np.array([[i+1,j],[i+1,j+1],[i,j+1]], float)/n
                uv = vertices.mean(axis=0); bary.append([1-uv.sum(), *uv])
    bary = np.asarray(bary); points = []; faces = []; weights = []
    assert len(bary) == n*n
    # Positive P1 pieces preserve clipped reverse-flux regions exactly.
    for face, xyz, q, flux in c.env.sampler.pieces:
        area = np.linalg.norm(np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0]))/2
        points.extend(bary@xyz); weights.extend(area/len(bary)*(bary@q)); faces.extend([face]*len(bary))
    return np.asarray(points), np.asarray(faces), np.asarray(weights)


def run_pool(fn, tasks, workers, label):
    start = time.time(); result = []
    with mp.get_context('fork').Pool(workers) as pool:
        for item in pool.imap_unordered(fn, tasks, chunksize=1):
            result.append(item)
            if len(result) % 10 == 0 or len(result) == len(tasks):
                print(label, len(result), '/', len(tasks), 'seconds', round(time.time()-start, 2), flush=True)
    return sorted(result, key=lambda x:x['first_id'])


def run_map(n, workers):
    p, f, w = triangle_quadrature(n); directory = f'map_n{n}'
    tasks = [(i, p[i:i+256], f[i:i+256], w[i:i+256], directory) for i in range(0, len(p), 256)]
    run_pool(map_shard, tasks, workers, directory)
    data = np.concatenate([np.load(x)['rows'] for x in sorted((Path(_RUN['output'])/directory).glob('map_*.npz'))])
    atomic_npz(Path(_RUN['output'])/(directory+'.npz'), rows=data)
    stats = {}
    for k, basin in enumerate(BASINS):
        m = data[:,5] == k; weight = data[m,4]; diameter = data[m,8]
        order = np.argsort(diameter); cumulative = np.cumsum(weight[order]); cumulative /= cumulative[-1] if len(cumulative) else 1
        stats[basin] = dict(count=int(m.sum()), flux_fraction=float(weight.sum()/w.sum()),
            geometry_pass_probability=float(np.average(data[m,9], weights=weight)) if len(weight) else None,
            center_on_plane_probability=float(np.average(data[m,10], weights=weight)) if len(weight) else None,
            Dmax_um_quantiles=dict(zip(['min','p10','median','p90','max'], np.interp([0,.1,.5,.9,1], cumulative, diameter[order]).tolist())) if len(weight) else None)
    atomic_json(Path(_RUN['output'])/(directory+'_summary.json'), dict(n=n, count=len(data),
        sum_quadrature_flux=float(w.sum()), official_flux=context().env.sampler.Q_m3_s, by_basin=stats,
        columns=['x_m','y_m','z_m','inlet_triangle','flux_weight_m3_s','basin_index','wall_distance_m','perimeter_distance_m','Dmax_um','geometry_pass_probability','center_on_plane_probability']))


def benchmark():
    results = []
    for workers in [1,2,4,8,16]:
        start = time.time();usage0=resource.getrusage(resource.RUSAGE_CHILDREN)
        metas = run_pool(event_shard, [(i,1,f'benchmark/w{workers}') for i in range(1,33)], workers, 'BENCHMARK')
        elapsed = time.time()-start
        usage1=resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu=usage1.ru_utime+usage1.ru_stime-usage0.ru_utime-usage0.ru_stime
        results.append(dict(workers=workers, events=32, seconds=elapsed, events_per_s=32/elapsed,
            worker_peak_rss_kib=max(x['peak_worker_rss_kib'] for x in metas),
            worker_cpu_seconds=cpu,mean_busy_cpu_cores=cpu/elapsed,
            worker_capacity_utilization_percent=100*cpu/(elapsed*workers)))
    fingerprints = []
    for workers in [1,2,4,8,16]:
        all_rows = []
        for p in sorted((Path(_RUN['output'])/f'benchmark/w{workers}').glob('*.json.gz')):
            with gzip.open(p, 'rt') as f:all_rows.extend(json.load(f)['rows'])
        fingerprints.append(canonical_hash(all_rows))
    if len(set(fingerprints)) != 1:raise ValueError('Parallel vs serial scientific parity failed')
    chosen = max(results, key=lambda r:r['events_per_s'])['workers']
    atomic_json(Path(_RUN['output'])/'ADMISSION_SCALING.json', dict(results=results, chosen_workers=chosen,
        scientific_parity=True, fingerprints=fingerprints, selection='MAXIMUM_OBSERVED_COMPLETED_EVENT_THROUGHPUT',
        shared_input_policy='FORK_COPY_ON_WRITE_READ_ONLY_FEM_AND_WALL_BVH'))
    return chosen


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['benchmark','map','events'])
    p.add_argument('--output',required=True);p.add_argument('--provenance',required=True)
    p.add_argument('--count',type=int,default=100000);p.add_argument('--workers',type=int,default=8)
    p.add_argument('--map-n',type=int,default=12);p.add_argument('--first',type=int,default=1)
    p.add_argument('--point-step',type=float,default=.2e-6);p.add_argument('--search-tolerance',type=float,default=1e-9)
    p.add_argument('--horizon-factor',type=float,default=1.)
    a=p.parse_args();provenance=json.loads(Path(a.provenance).read_text())
    config=dict(master_seed=MASTER_SEED, A_guard=512, B_guard=512, point_step_m=a.point_step,
        search_tolerance_m=a.search_tolerance, horizon_factor=a.horizon_factor,
        physics='UNCHANGED_P65_SPHERE', C_fallback='DISABLED', pair_interactions=False)
    initialize(a.output,provenance,config);out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    atomic_json(out/'RUN_CONFIG.json',config)
    atomic_json(out/'INLET_CAP_SEMANTICS.json',context().geometry.semantics())
    if a.action=='benchmark':benchmark()
    elif a.action=='map':run_map(a.map_n,a.workers)
    else:
        tasks=[(i,min(64,a.first+a.count-i),'events') for i in range(a.first,a.first+a.count,64)]
        start=time.time();metas=run_pool(event_shard,tasks,a.workers,'ADMISSION')
        atomic_json(out/'ADMISSION_RUN_RECEIPT.json',receipt('FORMAL_ADMISSION_COMPLETE',a.first,a.count,start,
            workers=a.workers,shards=len(metas),events_per_s=a.count/(time.time()-start)))


if __name__=='__main__':main()
