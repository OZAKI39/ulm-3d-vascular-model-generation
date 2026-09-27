"""Outcome-blind P9-A.4 prefix cohorts and deterministic production bookkeeping."""
from pathlib import Path
import gzip, hashlib, json, os

CORE_N = 500
MAX_N = 5000
DT = 0.001  # Explicit user revision: nominal time step = 1.0 ms.
HORIZONS = (3.0, 6.0, 12.0)
FLOW_SHA = '064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4'
BIRTHS_SHA = 'a36cbd03dcae3e4f98692d75f2e4ccf6639702a4eeb6d2a4613396742cc12156'
SOURCE_CONTRACT_SHA = '345d529f216df9aefc867ee17b557202ef8ee74c21e2b82636ce0ab5e702e48f'
REL = 'particle_3d/reports/particle9a5_formal_trajectories'


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)+'\n').encode()


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def content_sha(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def write_new(path, value):
    """Exclusive immutable record; only an identical existing record is reusable."""
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    data = canonical(value)
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError('Refuse replacement of frozen record: '+str(path))
        return digest(path)
    # Publish a complete record only, even if a process is interrupted.
    temp = path.with_name(path.name+'.pending.'+str(os.getpid()))
    with temp.open('xb') as f:
        f.write(data); f.flush(); os.fsync(f.fileno())
    try:
        os.link(temp, path)
    finally:
        temp.unlink()
    return hashlib.sha256(data).hexdigest()


def load_births(root):
    root = Path(root)
    p = root/'particle_3d/reports/particle9a4_population_inlet/data/inlet100k/accepted_births.json.gz'
    raw = gzip.decompress(p.read_bytes())
    if hashlib.sha256(raw).hexdigest() != BIRTHS_SHA:
        raise ValueError('P9-A.4 reviewed accepted ledger changed')
    population = json.loads(raw)
    events = population['events']
    if [e['particle_id'] for e in events] != list(range(1, len(events)+1)):
        raise ValueError('Accepted ordering is not contiguous')
    if any(b['source_event_id'] <= a['source_event_id'] for a, b in zip(events, events[1:])):
        raise ValueError('Source event ordering changed')
    if any(e['flow_sha256'] != FLOW_SHA for e in events):
        raise ValueError('OLD flow is prohibited')
    return population


def cohort(population, n, master_seed):
    if type(n) is not int or not CORE_N <= n <= min(MAX_N, len(population['events'])):
        raise ValueError('Formal cohort must contain a prefix of 500..5000 accepted births')
    events = population['events'][:n]
    return dict(schema='P9A5_FORMAL_ACCEPTED_PREFIX_V1', count=n,
                selection='FIRST_N_ACCEPTED_ORDER_NO_OUTCOME_INPUT', source_births_sha256=BIRTHS_SHA,
                source_identity=population['identity'], master_seed=master_seed,
                events=events, event_sha256=[content_sha(e) for e in events])


def require_prefix(old, new):
    if new['count'] < old['count'] or new['events'][:old['count']] != old['events']:
        raise ValueError('A formal extension may not replace or remove existing births')
    if [e['particle_id'] for e in new['events']] != list(range(1, new['count']+1)):
        raise ValueError('A formal extension may not skip IDs')


def next_size(n, all_observed, maximum=MAX_N):
    if maximum < CORE_N or maximum > MAX_N or n < CORE_N or n > maximum:
        raise ValueError('Invalid formal sample/resource bound')
    return n if all_observed or n == maximum else min(maximum, n+(250 if n < 2000 else 500))


def merge_rows(rows, expected_ids):
    rows = sorted(rows, key=lambda r: r['particle_id'])
    if [r['particle_id'] for r in rows] != list(expected_ids):
        raise ValueError('Missing, duplicate or unexpected formal result IDs')
    return rows


def completion_matches(folder, identity, event):
    folder = Path(folder); marker = folder/'COMPLETE.json'
    if not marker.is_file():
        return False
    try:
        m = json.loads(marker.read_text())
        if m['identity'] != identity or m['event_sha256'] != content_sha(event):
            return False
        required = {'trajectory.npz', 'trajectory.json', 'audit.jsonl.gz', 'metrics.json', 'support.json'}
        if not required <= m['files'].keys():
            return False
        for name, expected in m['files'].items():
            path = (folder/name).resolve()
            if not path.is_relative_to(folder.resolve()) or not path.is_file() or digest(path) != expected:
                return False
        return True
    except (KeyError, ValueError, OSError, TypeError):
        return False


class HorizonSchedule:
    """One iterator for one living stepper; only expand after reaching a boundary."""
    def __init__(self, horizons=HORIZONS, dt=DT, checkpoint=None):
        self.horizons = tuple(horizons); self.dt = dt; self.checkpoint = checkpoint
        if not self.horizons or any(b <= a for a, b in zip(self.horizons, self.horizons[1:])):
            raise ValueError('Horizons must increase strictly')
        if self.horizons[0] <= 0 or dt <= 0:
            raise ValueError('Positive horizon and dt required')
        self.used = self.horizons[0]; self.extensions = []

    def steps(self):
        previous = 0
        for j, horizon in enumerate(self.horizons):
            if j:
                old = self.horizons[j-1]
                if self.checkpoint:
                    self.checkpoint(old, horizon)
                self.extensions.append([old, horizon]); self.used = horizon
            end = round(horizon/self.dt)
            if abs(end*self.dt-horizon) > 1e-12:
                raise ValueError('Horizon must align with the unchanged nominal dt')
            for k in range(previous+1, end+1):
                yield k
            previous = end


def terminal_status(meta, stationary_supported, maximum_horizon):
    if meta.get('completed'):
        return 'COMPLETED'
    if meta.get('end_reason') == 'PHYSICAL_RESIDENCE_HORIZON_REACHED':
        return 'LONG_RESIDENCE_CENSORED' if meta.get('last_elapsed_time_s', 0) >= maximum_horizon-1e-12 else 'HORIZON_REACHED'
    candidate = any(s in str(meta.get('failure_detail')) for s in ['STATIONARY', 'ROUNDOFF_SCALE_STAGNATION'])
    return 'SUPPORTED_STATIONARY' if candidate and stationary_supported else 'SOLVER_FAILURE'
