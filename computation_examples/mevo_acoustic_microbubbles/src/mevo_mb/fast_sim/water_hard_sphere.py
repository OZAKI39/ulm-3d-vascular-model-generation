"""F10P opt-in adapter; inherits unchanged acoustic and real-cavity interfaces."""

import time
from .water_b2 import WaterB2
from .f10g_contact import resistance
from .hard_sphere import hard_step


class WaterHardSphere(WaterB2):
    def __init__(self, *args, policy, **kwargs):
        super().__init__(*args, **kwargs)
        self.policy = policy

    def step(self, dt):
        start = time.perf_counter()
        f = self.fields(self.x)
        t = time.perf_counter()
        op = resistance(self.x, self.radius)
        assembly = time.perf_counter() - t
        result = hard_step(self.x, f["total_force_N"], self.radius, self.zeta, dt, op["A"], self.policy)
        walls = [self.wall_check(result[k]) for k in ["prediction_m", "positions_m"]]
        gates = result["gates"].copy()
        gates["WALL_QUALITY"] = all(w["error_bits"] == 0 and w["correction_m"] == 0 for w in walls)
        reasons = [k for k, v in gates.items() if not v]
        gates["FULL_STEP_ACCEPTANCE"] = not reasons
        if not reasons:
            self.x = result["positions_m"].copy()
        return dict(
            accepted=not reasons,
            reasons=reasons,
            gates=gates,
            result=result,
            fields=f,
            wall_correction_m=max(w["correction_m"] for w in walls),
            wall_error_bits=max(w["error_bits"] for w in walls),
            resistance_s=assembly,
            full_step_s=time.perf_counter() - start,
        )
