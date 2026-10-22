"""P9-A.4 counter-addressed Poisson proposals and single-shot source thinning.

Legacy schedulers, samplers, admission and dynamics are imported, never patched.
The deterministic reducer owns physical time and accepted IDs. Workers only
produce independent marks; they cannot condition a proposal on earlier outcomes.
"""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import math
import os
import numpy as np

MODEL = 'STEADY_CONTINUOUS_INFUSION_FINITE_SIZE_FLUX_INLET_V1'
NEW_FLOW_SHA256 = '064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4'
P9A4_ARRIVAL, P9A4_DIAMETER, P9A4_POSITION, P9A4_ORIENTATION = 940, 941, 942, 943
ROLES = (P9A4_ARRIVAL, P9A4_DIAMETER, P9A4_POSITION, P9A4_ORIENTATION)


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(name+' must be a positive integer')
    return int(value)


def counter_rng(master_seed, source_event_id, role):
    if isinstance(master_seed, (bool, np.bool_)) or not isinstance(master_seed, (int, np.integer)) or master_seed < 0:
        raise ValueError('master_seed must be a nonnegative integer')
    source_event_id = positive_integer(source_event_id, 'source_event_id')
    if role not in ROLES:
        raise ValueError('Unknown P9-A.4 RNG role')
    return np.random.default_rng([int(master_seed), source_event_id, int(role)])


def poisson_interval(uniform, rate):
    """Inverse CDF; zero rate has no next event, never an infinite ledger row."""
    uniform, rate = float(uniform), float(rate)
    if not math.isfinite(rate) or rate < 0:
        raise ValueError('Finite nonnegative source rate required')
    if not math.isfinite(uniform) or not 0 <= uniform < 1:
        raise ValueError('Uniform must lie in [0,1)')
    if rate == 0:
        raise StopIteration('ZERO_SOURCE_RATE')
    interval = -math.log1p(-uniform)/rate
    if not math.isfinite(interval):
        raise FloatingPointError('Source interval exceeds finite float64 time range')
    return interval


class ContinuousInfusionSource:
    """One size draw, one full-flux position draw and one admission.check call.

    The synthetic switch permits analytical benchmark fixtures only. Real
    production setup independently verifies both file SHA and mesh provenance.
    """
    def __init__(self, *, sampler, distribution, checker, concentration_m3,
                 flow_sha256, master_seed, source_contract_sha256, synthetic=False, scientific_identity=None):
        if not synthetic and flow_sha256 != NEW_FLOW_SHA256:
            raise ValueError('P9A4_REQUIRES_AUTHORITATIVE_NEW_NETWORK_H0_FLOW')
        concentration_m3 = float(concentration_m3)
        q = float(sampler.Q_m3_s)
        if not math.isfinite(concentration_m3) or concentration_m3 < 0:
            raise ValueError('Finite nonnegative modeled concentration required')
        if not math.isfinite(q) or q <= 0:
            raise ValueError('Finite positive inlet flux required')
        rate = concentration_m3*q
        if not math.isfinite(rate) or (concentration_m3 > 0 and rate == 0):
            raise ValueError('Source rate overflow or underflow')
        counter_rng(master_seed, 1, P9A4_ARRIVAL)
        self.sampler, self.distribution, self.checker = sampler, distribution, checker
        self.rate, self.master_seed = rate, int(master_seed)
        self.identity = dict(model=MODEL, flow_sha256=flow_sha256,
            source_contract_sha256=source_contract_sha256, concentration_m3=concentration_m3,
            Q_in_m3_s=q, lambda_source_s_inv=rate, master_seed=self.master_seed,
            rng='PCG64_SEED_SEQUENCE_MASTER_SOURCE_EVENT_ROLE', roles=list(ROLES),
            active_particles='EMPTY_INDEPENDENT_SINGLE_MB_POPULATION', synthetic=bool(synthetic),
            scientific_identity=deepcopy(scientific_identity or {}))
        self.identity_sha256 = hashlib.sha256(canonical_bytes(self.identity)).hexdigest()

    def proposal(self, source_event_id):
        sid = positive_integer(source_event_id, 'source_event_id')
        arrival_u = float(counter_rng(self.master_seed, sid, P9A4_ARRIVAL).random())
        dt = poisson_interval(arrival_u, self.rate)
        diameter, draw = self.distribution.sample(counter_rng(self.master_seed, sid, P9A4_DIAMETER))
        diameter = float(diameter)
        if not math.isfinite(diameter) or not 0 < diameter <= 4e-6:
            raise ValueError('P9A4 source diameter must be in (0,4 um]')
        positions, triangles = self.sampler.sample(counter_rng(self.master_seed, sid, P9A4_POSITION))
        position = np.asarray(positions[0], dtype=np.float64)
        if position.shape != (3,) or not np.isfinite(position).all():
            raise ValueError('Finite 3D source position required')
        orientation = counter_rng(self.master_seed, sid, P9A4_ORIENTATION).standard_normal(4)
        orientation /= np.linalg.norm(orientation)
        # Temporary tag used by the pure check; active={} means no particle-pair
        # interaction depends on it. Only the parent reducer assigns particle_id.
        event = dict(species='MB', particle_id=sid, radius_m=diameter/2, q=orientation.tolist())
        particle, admission_status, detail = self.checker.check(event, position, {})
        accepted = particle is not None
        if accepted != (admission_status == 'ACCEPTED'):
            raise ValueError('Inconsistent authoritative admission result')
        wall = getattr(self.checker, 'wall', None)
        clearance = None if wall is None else float(wall.nearest_center_triangle(position)[1]-diameter/2)
        row = dict(source_event_id=sid, particle_id=None, arrival_uniform=arrival_u,
            arrival_delta_t_s=dt, diameter_m=diameter, diameter_um=diameter*1e6,
            radius_m=diameter/2, position_m=position.tolist(), inlet_triangle_id=int(triangles[0]),
            q=orientation.tolist(), status='ACCEPTED_BIRTH' if accepted else 'REJECTED_SOURCE_EVENT',
            rejection_reason=None if accepted else admission_status, admission_status=admission_status,
            diameter_source_draw=draw, clearance_m=clearance, admission_detail=detail,
            position_draw_count=1, diameter_draw_count=1, admission_check_count=1,
            identity_sha256=self.identity_sha256)
        canonical_bytes(row)  # fail explicitly on NaN / Inf before publication
        return row


class PopulationLedger:
    """Ordered reduction, with rejected source events retained permanently."""
    def __init__(self, identity):
        self.identity = deepcopy(identity)
        self.identity_sha256 = hashlib.sha256(canonical_bytes(identity)).hexdigest()
        self.rows = []
        self.time_s = 0.
        self.accepted = 0

    @property
    def next_source_event_id(self):
        return len(self.rows)+1

    def append(self, proposal):
        row = deepcopy(proposal)
        if row['source_event_id'] != self.next_source_event_id or row['identity_sha256'] != self.identity_sha256:
            raise ValueError('Missing, duplicate, unordered source ID or changed source contract')
        if row['status'] not in ('ACCEPTED_BIRTH', 'REJECTED_SOURCE_EVENT'):
            raise ValueError('Unknown proposal status')
        dt = float(row['arrival_delta_t_s'])
        if not math.isfinite(dt) or dt < 0:
            raise ValueError('Invalid interarrival time')
        expected_u = float(counter_rng(self.identity['master_seed'], row['source_event_id'], P9A4_ARRIVAL).random())
        if row['arrival_uniform'] != expected_u or dt != poisson_interval(expected_u, self.identity['lambda_source_s_inv']):
            raise ValueError('Arrival key or inverse-transform mismatch')
        if row['particle_id'] is not None:
            raise ValueError('Workers cannot allocate particle IDs')
        if any(row[k] != 1 for k in ('position_draw_count', 'diameter_draw_count', 'admission_check_count')):
            raise ValueError('Physical source proposals cannot retry')
        time_s = self.time_s+dt  # fixed source ID order, independent of worker partition
        if not math.isfinite(time_s) or (dt > 0 and time_s <= self.time_s):
            raise FloatingPointError('Physical source clock has exhausted float64 resolution')
        accepted = row['status'] == 'ACCEPTED_BIRTH'
        if accepted != (row['admission_status'] == 'ACCEPTED') or accepted != (row['rejection_reason'] is None):
            raise ValueError('Invalid accepted/rejected event semantics')
        row['proposal_time_s'] = time_s
        row['particle_id'] = self.accepted+1 if accepted else None
        canonical_bytes(row)
        self.time_s = time_s
        self.accepted += int(accepted)
        self.rows.append(row)
        return row

    def merge(self, proposals):
        proposals = sorted(proposals, key=lambda r:r['source_event_id'])
        ids = [r['source_event_id'] for r in proposals]
        if ids != list(range(self.next_source_event_id, self.next_source_event_id+len(ids))):
            raise ValueError('Worker merge contains missing or duplicate IDs')
        for proposal in proposals:
            self.append(proposal)

    def birth_events(self, count=None):
        if count is not None:
            positive_integer(count, 'count')
        births = []
        for row in self.rows:
            if row['status'] != 'ACCEPTED_BIRTH':
                continue
            d = row['diameter_m']
            births.append(dict(particle_id=row['particle_id'], source_event_id=row['source_event_id'],
                species='MB', birth_time_s=row['proposal_time_s'], scheduled_time_s=row['proposal_time_s'],
                diameter_m=d, diameter_um=row['diameter_um'], radius_m=d/2,
                volume_m3=4*math.pi*(d/2)**3/3, q=row['q'], birth_center_m=row['position_m'],
                anchor_m=row['position_m'], anchor_triangle=row['inlet_triangle_id'],
                inlet_triangle_id=row['inlet_triangle_id'], clearance_m=row['clearance_m'],
                entry_transition_time_s=0., admission_strategy=MODEL, admission_status='ACCEPTED',
                source_identity_sha256=self.identity_sha256, flow_sha256=self.identity['flow_sha256'],
                source_distribution_contract_sha256=self.identity['source_contract_sha256']))
            if count is not None and len(births) == count:
                break
        if count is not None and len(births) != count:
            raise ValueError('Insufficient accepted population; never draw extra accepted-only proposals')
        return births

    def state(self):
        return dict(schema='P9A4_POPULATION_CHECKPOINT_V1', identity=self.identity,
            next_source_event_id=self.next_source_event_id, next_particle_id=self.accepted+1,
            physical_source_time_s=self.time_s, ledger_count=len(self.rows),
            ledger_sha256=hashlib.sha256(b''.join(canonical_bytes(r) for r in self.rows)).hexdigest())

    def checkpoint(self, folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=False)
        with (folder/'proposal_ledger.jsonl').open('xb') as stream:
            for row in self.rows:
                stream.write(canonical_bytes(row))
            stream.flush(); os.fsync(stream.fileno())
        with (folder/'state.json').open('xb') as stream:
            stream.write(canonical_bytes(self.state())); stream.flush(); os.fsync(stream.fileno())
        manifest = {name:sha256(folder/name) for name in ('proposal_ledger.jsonl','state.json')}
        (folder/'manifest.tmp').write_bytes(canonical_bytes(manifest))
        os.replace(folder/'manifest.tmp', folder/'manifest.json')  # completion marker last

    @classmethod
    def restore(cls, folder, expected_identity):
        folder = Path(folder)
        manifest = json.loads((folder/'manifest.json').read_text())
        if set(manifest) != {'proposal_ledger.jsonl','state.json'}:
            raise ValueError('Incomplete checkpoint schema')
        for name, expected in manifest.items():
            if sha256(folder/name) != expected:
                raise ValueError('Checkpoint bytes changed: '+name)
        state = json.loads((folder/'state.json').read_text())
        if state['identity'] != expected_identity:
            raise ValueError('Checkpoint source/flow/RNG identity mismatch')
        ledger = cls(expected_identity)
        with (folder/'proposal_ledger.jsonl').open() as stream:
            for line in stream:
                expected = json.loads(line)
                proposal = deepcopy(expected)
                proposal.pop('proposal_time_s'); proposal['particle_id'] = None
                actual = ledger.append(proposal)
                if canonical_bytes(actual) != canonical_bytes(expected):
                    raise ValueError('Checkpoint event reduction mismatch')
        if ledger.state() != state:
            raise ValueError('Checkpoint counters or rejected events changed')
        return ledger
