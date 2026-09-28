"""One authoritative clock; all active particles share each accepted interval."""
from collections import Counter
import numpy as np
from ..rbc_orientation import advance_orientation
from .state import ActiveParticlePool
from .reference import ReferenceDynamics


class PopulationEngine:
    def __init__(self, environment, *, dt=.0005, pairs=False, realization_id=0, recorder=None):
        if not np.isfinite(dt) or dt <= 0: raise ValueError('Positive finite nominal dt required')
        self.env = environment; self.dt = float(dt); self.time = 0.
        self.realization_id = int(realization_id)
        self.pool = ActiveParticlePool()
        self.dynamics = ReferenceDynamics(environment, pairs=pairs)
        self.recorder = recorder
        self.snapshots = []; self.terminals = []; self.timeline = []
        self.births = {}; self.source_ledger = []; self.counts = Counter()
        self.maximum_active = 0
        self._simultaneous_outlets = {}
        self._previous_pair_nearfield = set(); self._previous_pair_contacts = set()

    def insert_accepted(self, particle, event):
        if event.get('realization_id', self.realization_id) != self.realization_id:
            raise ValueError('Particle belongs to a different realization')
        self.pool.insert(particle, event, self.time)
        self.births[particle.particle_id] = dict(event)
        self.counts['accepted_birth_count'] += 1
        self.maximum_active = max(self.maximum_active, self.pool.n)
        self.record()

    def record(self):
        row = self.pool.snapshot(self.time)
        if self.recorder is None: self.snapshots.append(row)
        else: self.recorder.snapshot(row)

    def _accepted(self, world):
        interval = world.time_s - self.time
        if interval <= 0: raise ValueError('Accepted state must advance global physical time')
        for k, pid in enumerate(self.pool.particle_id[:self.pool.n]):
            pid = int(pid)
            if self.env.classifier is not None and world.boundary_events:
                hit=self.env.classifier.first_event(self.pool.x[k],world.shapes[pid].center_m)
                # All particles reached the same physical end time. Endpoint
                # events belong together; do not advance a tied outlet a second time.
                if hit is not None and hit.role.startswith('OUTLET_') and abs(hit.segment_fraction-1.)<=64*np.finfo(float).eps:
                    self._simultaneous_outlets[pid]=hit.role
            self.pool.x[k] = world.shapes[pid].center_m
            self.pool.v[k] = world.velocities[pid]
            self.pool.omega[k] = world.omegas[pid]
            self.pool.q[k] = advance_orientation(self.pool.q[k], self.pool.omega[k], interval)
            self.pool.age[k] = world.time_s - self.pool.birth_time[k]
            if pid in world.wall_gaps: self.pool.wall_gap[k] = world.wall_gaps[pid].gap_m
            self.pool.tetra_id[k] = self.env.field.locate(self.pool.x[k])[0]
        self.time = world.time_s
        near={tuple(sorted((r['particle_i_id'],r['particle_j_id'])))
              for r in world.projection.get('potential_blocks',[]) if r['kind']=='PAIR' and r['active']}
        contact={tuple(sorted(r['canonical_id'][-2:])) for r in world.projection.get('handoff_constraints',[])
                 if r.get('canonical_id',('', ''))[1]=='PAIR'}
        self.counts['pair_nearfield_event_count']+=len(near-self._previous_pair_nearfield)
        self.counts['pair_collision_event_count']+=len(contact-self._previous_pair_contacts)
        self._previous_pair_nearfield=near;self._previous_pair_contacts=contact
        self.counts['maximum_simultaneous_pair_contacts']=max(self.counts['maximum_simultaneous_pair_contacts'],len(contact))
        self.record()
        if self.recorder is not None and hasattr(self.recorder, 'audit'):
            self.recorder.audit(world)

    def _remove_outlets(self, world):
        events=dict(world.boundary_events,**{})
        events.update(self._simultaneous_outlets);self._simultaneous_outlets.clear()
        for pid, role in sorted(events.items()):
            if not role.startswith('OUTLET_'): raise ValueError('Only official open outlets terminate as exits')
            birth = self.births[pid]
            row = dict(particle_id=pid, source_event_id=birth['source_event_id'],
                realization_id=self.realization_id, birth_time_s=birth['birth_time_s'],
                physical_time_s=self.time, residence_time_s=self.time-birth['birth_time_s'],
                status='EXITED', terminal_reason='FIRST_CENTER_CROSSING_OF_OFFICIAL_OPEN_OUTLET',
                outlet=role, exit_time_s=self.time)
            self.terminals.append(row)
            self.counts['exited_' + role] += 1
            if self.recorder is not None: self.recorder.terminal(row)
        self.pool.remove(events)

    def advance(self, target):
        if not np.isfinite(target) or target < self.time: raise ValueError('Global time cannot decrease')
        while self.time < target:
            if not self.pool.n:
                self.time = float(target); break
            old = self.time
            world = self.dynamics.advance(self.pool, self.time, float(target), self._accepted)
            self._remove_outlets(world)
            if self.time <= old: raise RuntimeError('Population made no physical progress')
        self.accounting()

    def step_to(self, target):
        while self.time < target:
            # The global grid, not particle age, defines the nominal time boundary.
            ratio=self.time/self.dt
            tick = (int(np.floor(ratio + 8*abs(np.spacing(ratio)))) + 1)*self.dt
            if tick <= self.time: tick = self.time + self.dt
            self.advance(min(float(target), tick))

    def accounting(self):
        accepted = self.counts['accepted_birth_count']
        exited = sum(v for k, v in self.counts.items() if k.startswith('exited_'))
        removed = sum(self.counts[k] for k in ['supported_stationary_terminated_count', 'solver_failed_count', 'censored_count'])
        if accepted != self.pool.n + exited + removed:
            raise AssertionError('POPULATION_ACCOUNTING_FAILED')
        if self.counts['source_event_count']:
            rejected=sum(self.counts[k] for k in ['wall_rejected_source_count','pair_rejected_inlet_count','other_rejected_source_count'])
            if self.counts['source_event_count']!=accepted+rejected:
                raise AssertionError('SOURCE_ACCOUNTING_FAILED')
        row = dict(self.counts, physical_time_s=self.time, active_count=self.pool.n,
                   maximum_simultaneous_active_particles=self.maximum_active,
                   conservation_pass=True)
        self.timeline.append(row)
        return row
