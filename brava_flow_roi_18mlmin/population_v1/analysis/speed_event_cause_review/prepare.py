"""Read original server evidence and only the four frozen windows. No dynamics."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
POP = HERE.parents[1]
BASE = POP.parent
RUN = POP / 'runs/realization_0000_contactfix'
TRACKS = POP / 'analysis/realization_0000_contactfix/tracks'
SHA = {}


def register(path):
    with path.open('rb') as stream:
        SHA[str(path.relative_to(BASE))] = hashlib.file_digest(stream, 'sha256').hexdigest()
    return path


def read(path):
    return json.loads(register(path).read_text())


def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def main():
    plan = json.loads((HERE / 'FROZEN_PLAN.json').read_text())
    identity = read(RUN / 'IDENTITY.json')
    manifest = read(RUN / 'MANIFEST.json')
    source = POP / 'source/particle_3d/src/particle_3d'
    historical = json.loads((HERE / 'HISTORICAL_SOURCE_SHA.json').read_text())
    differences = []
    for name, expected in historical['files'].items():
        path = register(source / name)
        if SHA[str(path.relative_to(BASE))] != expected:
            differences.append(name)
    if differences != ['population_v1/run.py']:
        raise ValueError(('HISTORICAL_SOURCE_IDENTITY_MISMATCH', differences))
    for name, expected in identity['code_files_sha256'].items():
        path = register(source / 'population_v1' / name)
        if name != 'run.py' and SHA[str(path.relative_to(BASE))] != expected:
            raise ValueError(('ORIGINAL_RUN_SOURCE_MISMATCH', name))
    for key, path in dict(mesh=identity['flow_provenance']['mesh'],
                          flow=identity['flow_provenance']['flow'],
                          **identity['flow_provenance']['boundaries']).items():
        p = register(Path(path))
        assert SHA[str(p.relative_to(BASE))] == identity['flow_provenance']['sha256'][key]
    historical_results = {}
    for name in ['SUMMARY', 'largest_step_increments', 'largest_increment_context', 'display_review']:
        historical_results[name] = read(POP / 'analysis/failure_speed_review_20260929' / (name + '.json'))
    assert historical_results['largest_step_increments']['top_per_track'][0]['particle_id'] == 822
    for name in ['SCAN_SUMMARY', 'DIAGNOSIS_SUMMARY', 'CPU_STATE_REPLAY', 'WALL_RELEASE_AND_FEM']:
        read(POP / 'analysis/acceleration_diagnosis' / (name + '.json'))
    read(BASE / 'reports/ACTIVE_FLOW.json')
    read(BASE / 'production_workflow/PARTICLE_CFD_INPUT_CONTRACT.json')
    catalog = {int(r['particle_id']): r for r in csv.DictReader(register(
        POP / 'analysis/realization_0000_contactfix/trajectory_catalog.csv').open())}
    events = {}
    for event in plan['events']:
        pid = event['particle_id']
        with np.load(register(TRACKS / f'mb_{pid:06d}.npz')) as z:
            t = z['physical_time_s']
            a, b = np.searchsorted(t, event['window_s'], side='left')
            data = {k: z[k][a:b] for k in z.files}
        assert len(data['physical_time_s']) > 2
        np.savez_compressed(HERE / f'original_{pid}.npz', **data)
        events[pid] = dict(**event, radius_m=float(catalog[pid]['radius_m']),
                           birth_time_s=float(catalog[pid]['birth_time_s']),
                           track_start_index=int(a), samples=len(data['physical_time_s']),
                           first_recorded_s=float(data['physical_time_s'][0]),
                           last_recorded_s=float(data['physical_time_s'][-1]),
                           pair_rows=[], rejected_rows=[], handoff_rows=[], cpu_component_rows=[],
                           accepted_audit_rows=0, initial_state=None,
                           minimum_neighbor_surface_gap_m=None,
                           minimum_swept_margin_above_pair_shell_m=None)
    births = {}
    with gzip.open(register(RUN / 'source.jsonl.gz'), 'rt') as f:
        for line in f:
            row = json.loads(line)
            if row.get('particle_id') in events:
                births[row['particle_id']] = row
    with gzip.open(register(RUN / 'interaction_audit.jsonl.gz'), 'rt') as f:
        for line in f:
            row = json.loads(line)
            t = row.get('evaluated_at_time_s', row['physical_time_s'])
            for pid, e in events.items():
                if not e['window_s'][0] <= t <= e['window_s'][1]:
                    continue
                if row.get('event') == 'REJECTED_TRIAL':
                    e['rejected_rows'].append(row)
                    continue
                e['accepted_audit_rows'] += 1
                if any(pid in pair for pair in row.get('pair_ids', [])):
                    e['pair_rows'].append(row)
                handoff = row.get('handoff_event')
                if handoff and pid in [handoff.get('particle_id'), handoff.get('particle_j_id')]:
                    e['handoff_rows'].append(row)
                cs = row.get('cpu_components')
                if cs and any(pid in c['particle_ids'] for c in cs['components']):
                    e['cpu_component_rows'].append(row)
    chunks = [c for c in manifest['chunks'] if any(
        c['last_time_s'] >= e['window_s'][0] and c['first_time_s'] <= e['window_s'][1]
        for e in events.values())]
    def load_chunk(c):
        path = register(RUN / c['path'])
        assert SHA[str(path.relative_to(BASE))] == c['file_sha256']
        with np.load(path) as z:
            return {k: z[k] for k in z.files}
    previous = {}
    neighbors = {pid: [] for pid in events}
    # Two independent small-file readers; bounded queue (no full-run load).
    with ThreadPoolExecutor(2) as pool:
        for start in range(0, len(chunks), 2):
            for z in pool.map(load_chunk, chunks[start:start + 2]):
                for index, t in enumerate(z['state_time_s']):
                    for pid, e in events.items():
                        if not e['window_s'][0] <= t <= e['window_s'][1]:
                            continue
                        a, b = z['state_offsets'][index:index + 2]
                        ids = z['particle_id'][a:b]
                        ix = np.flatnonzero(ids == pid)
                        if not len(ix):
                            continue
                        k = int(ix[0])
                        x, r = z['x'][a:b], z['radius'][a:b]
                        if e['initial_state'] is None:
                            e['initial_state'] = {name: z[name][a+k].tolist() for name in z
                                if name not in ['state_time_s', 'state_offsets']}
                            e['initial_state']['physical_time_s'] = float(t)
                        separation = np.linalg.norm(x - x[k], axis=1) - r - r[k]
                        separation[k] = np.inf
                        minimum = float(separation.min())
                        e['minimum_neighbor_surface_gap_m'] = min(e['minimum_neighbor_surface_gap_m'] or np.inf, minimum)
                        sweep = None
                        if pid in previous:
                            oldt, oldids, oldx, oldr, oldk = previous[pid]
                            if t > oldt:
                                common, ia, ib = np.intersect1d(oldids, ids, return_indices=True)
                                mask = common != pid
                                ia, ib = ia[mask], ib[mask]
                                relative = oldx[ia] - oldx[oldk]
                                delta = (x[ib] - x[k]) - relative
                                den = np.einsum('ij,ij->i', delta, delta)
                                fraction = np.divide(-np.einsum('ij,ij->i', relative, delta), den,
                                    out=np.zeros_like(den), where=den > 0).clip(0, 1)
                                reference = np.minimum(r[k], r[ib])
                                shell = np.maximum(.05*reference, np.maximum(2e-9, .001*reference))
                                sweep = float((np.linalg.norm(relative + fraction[:, None]*delta, axis=1)
                                               - r[k] - r[ib] - shell).min())
                                e['minimum_swept_margin_above_pair_shell_m'] = min(
                                    e['minimum_swept_margin_above_pair_shell_m'] or np.inf, sweep)
                        neighbors[pid].append((float(t), minimum, sweep))
                        previous[pid] = (float(t), ids.copy(), x.copy(), r.copy(), k)
    for pid, e in events.items():
        e['birth'] = births[pid]
        e['isolation_confirmed'] = (not e['pair_rows'] and
            e['minimum_swept_margin_above_pair_shell_m'] > 0)
        e['state_restart_limit'] = None if pid == 822 else (
            'Saved x/v/omega/q and event times exist, but safety provider-call counters, '
            'no-acceptance counters and nominal/history buffers are not checkpointed. '
            'No invented safety state; RHS reevaluation only in this review.')
        save(f'ledger_{pid}.json', e)
        save(f'neighbors_{pid}.json', neighbors[pid])
        print('WINDOW', pid, 'samples', e['samples'], 'pair_rows', len(e['pair_rows']),
              'swept_margin_um', e['minimum_swept_margin_above_pair_shell_m']*1e6, flush=True)
    save('SERVER_EVIDENCE.json', dict(source_differences=differences, historical_results=historical_results,
         chunks_read=len(chunks), read_scope='Only overlapping original state chunks, four tracks and ledgers',
         held_velocity_semantics=manifest['velocity']))
    save('INPUT_SHA256.json', SHA)


if __name__ == '__main__':
    main()
