"""Bounded diagnostic observer; inherited production kernels are unchanged.

Host readbacks are intentional for <=50us diagnosis, never a timing reference.
The independent physics configuration permits 24 rounds without changing F2.
"""

from dataclasses import dataclass, asdict
import numpy as np
from scipy.spatial.distance import pdist
from .nearfield import F2Physics
from .f6_engine import F6Simulator
from . import kernels as k


@dataclass(frozen=True)
class ContactDiagnosticPhysics(F2Physics):
    def __post_init__(self):
        if self.contact_iterations not in (12, 24):
            raise ValueError("Only preregistered12/24 rounds")
        F2Physics(**(asdict(self) | {"contact_iterations": 12}))


def distances(x, radius):
    d = pdist(x)
    return dict(min_center_distance_m=float(d.min()), max_overlap_m=float(max(0, 2 * radius - d.min())))


class ContactTrace(F6Simulator):
    def __init__(self, positions, physics, *args, **kwargs):
        self.trace = None
        # F3's constructor freezes twelve rounds. Allocate its unchanged buffers
        # under that configuration, then select the explicitly registered local
        # diagnostic loop count. No historical class or acceptance test changes.
        allocated = F2Physics(**(asdict(physics) | {"contact_iterations": 12}))
        super().__init__(positions, allocated, *args, **kwargs)
        self.p = physics

    def fields(self):
        super().fields()
        if self.trace is not None:
            self.trace["primary_force_N"] = self.force.numpy().copy()

    def step(self, dt):
        self.trace = dict(
            initial_positions_m=self.x.numpy().copy(), dt_s=dt, rounds=[], wall_passes=[], lubrication=[]
        )
        return super().step(dt)

    def record_pair(self, x):
        after = x.numpy().copy()
        dx = after - self._pair_before
        self.trace["rounds"].append(
            dict(
                round=len(self.trace["rounds"]) + 1,
                positions_m=after,
                correction_m=dx,
                mean_correction_m=dx.mean(0),
                max_correction_m=float(np.linalg.norm(dx, axis=1).max()),
                mean_correction_path_m=float(np.linalg.norm(dx, axis=1).mean()),
                **distances(after, self.p.radius_m),
            )
        )
        self._recorded = True

    def wall_pass(self, x):
        if self.trace is None:
            return super().wall_pass(x)
        # F4 calls this before launch(contact) returns: record pre-wall PBD.
        if self._stage == "pair" and not self._recorded:
            self.record_pair(x)
        if self._stage == "predict":
            self.trace["predicted_positions_m"] = x.numpy().copy()
        before = x.numpy().copy()
        result = super().wall_pass(x)
        self.trace["wall_passes"].append(
            dict(after_pair_round=self._pair_contact_pass, correction_m=x.numpy() - before)
        )
        return result

    def launch(self, kernel, dim, args):
        if self.trace is None:
            return super().launch(kernel, dim, args)
        if kernel is k.predict:
            self._stage = "predict"
            self.trace.update(
                total_force_N=self.force.numpy().copy(),
                free_velocity_m_s=self.free.numpy().copy(),
                flow_velocity_m_s=self.flow_query.velocity.numpy().copy(),
                local_velocity_m_s=args[1].numpy().copy(),
            )
        elif kernel is k.contact:
            self._stage = "pair"
            self._pair_before = args[0].numpy().copy()
            self._recorded = False
        elif kernel is k.accept:
            self._stage = "accept"
            self.trace["proposal_positions_m"] = args[0].numpy().copy()
        result = super().launch(kernel, dim, args)
        if kernel is k.contact and not self._recorded:
            self.record_pair(args[1])
        elif kernel is k.local_mobility:
            self.trace["lubrication"].append(args[3].numpy().copy())
        return result
