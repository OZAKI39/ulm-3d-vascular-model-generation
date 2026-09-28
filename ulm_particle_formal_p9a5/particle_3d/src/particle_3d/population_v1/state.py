"""Structure-of-arrays active state, stable IDs and explicit terminal removal."""
import numpy as np
from ..particle_shapes import Sphere

ACTIVE, EXITED, STATIONARY, FAILED, CENSORED = range(5)


class ActiveParticlePool:
    def __init__(self, capacity=64):
        self.n = 0
        self.used_ids = set()
        self.capacity = max(1, int(capacity))
        for name in ['particle_id', 'source_event_id', 'realization_id', 'tetra_id', 'status']:
            setattr(self, name, np.zeros(self.capacity, dtype=np.int64))
        for name in ['birth_time', 'age', 'radius', 'diameter', 'wall_gap']:
            setattr(self, name, np.zeros(self.capacity, dtype=np.float64))
        for name in ['x', 'v', 'omega']:
            setattr(self, name, np.zeros((self.capacity, 3), dtype=np.float64))
        self.q = np.zeros((self.capacity, 4), dtype=np.float64)

    def arrays(self):
        return {k: v for k, v in vars(self).items() if isinstance(v, np.ndarray)}

    def insert(self, particle, event, physical_time):
        pid = int(particle.particle_id)
        if pid in self.used_ids:
            raise ValueError('STABLE_ID_NEVER_REUSED')
        if physical_time != event['birth_time_s']:
            raise ValueError('BIRTH_MUST_OCCUR_AT_EXACT_PHYSICAL_TIME')
        if self.n == self.capacity:
            self.capacity *= 2
            for name, array in self.arrays().items():
                enlarged = np.zeros((self.capacity, *array.shape[1:]), dtype=array.dtype)
                enlarged[:self.n] = array[:self.n]
                setattr(self, name, enlarged)
        k = self.n; self.n += 1; self.used_ids.add(pid)
        values = dict(particle_id=pid, source_event_id=event['source_event_id'],
                      realization_id=event.get('realization_id', 0), tetra_id=-1,
                      status=ACTIVE, birth_time=physical_time, age=0.,
                      radius=particle.radius, diameter=2*particle.radius,
                      x=particle.position, v=particle.velocity, omega=particle.omega,
                      q=particle.q, wall_gap=event.get('clearance_m', 0.))
        for name, value in values.items(): getattr(self, name)[k] = value

    def remove(self, ids):
        keep = ~np.isin(self.particle_id[:self.n], list(ids))
        count = int(keep.sum())
        for value in self.arrays().values(): value[:count] = value[:self.n][keep]
        self.n = count

    def shapes(self):
        return {int(self.particle_id[k]): Sphere(self.x[k], self.radius[k]) for k in range(self.n)}

    def snapshot(self, time):
        result = {k: v[:self.n].copy() for k, v in self.arrays().items()}
        result['physical_time_s'] = float(time)
        return result
