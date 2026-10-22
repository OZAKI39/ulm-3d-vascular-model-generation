"""Opt-in SQP adapter retaining original fields, forces and cavity wall queries."""

import time
from .water_b2 import WaterB2
from .accelerated_contact import resistance
from .accelerated_sequential import sequential_step


class WaterAccelerated(WaterB2):
    def __init__(self, *args, policy, backend="cpu", **kwargs):
        super().__init__(*args, **kwargs)
        self.policy = policy
        self.backend = backend

    def step(self, dt):
        begin = time.perf_counter()
        f = self.fields(self.x)
        t = time.perf_counter()
        op = resistance(self.x, self.radius)
        assembly = time.perf_counter() - t
        r = sequential_step(self.x, f["total_force_N"], self.radius, self.zeta, dt, op["A"], self.policy, backend=self.backend)
        walls = [self.wall_check(r[k]) for k in ["prediction_m", "positions_m"]]
        gates = r["gates"].copy()
        gates["WALL_PASS"] = all(w["error_bits"] == 0 and w["correction_m"] == 0 for w in walls)
        reasons = [k for k, v in gates.items() if not v]
        gates["FULL_STEP_ACCEPTED"] = not reasons
        if not reasons:
            self.x = r["positions_m"].copy()
        return dict(
            accepted=not reasons,
            reasons=reasons,
            gates=gates,
            result=r,
            fields=f,
            wall_correction_m=max(w["correction_m"] for w in walls),
            wall_error_bits=max(w["error_bits"] for w in walls),
            resistance_s=assembly,
            full_step_s=time.perf_counter() - begin,
        )
