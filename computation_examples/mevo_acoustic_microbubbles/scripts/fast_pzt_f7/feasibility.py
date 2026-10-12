"""One exact failed-state feasibility replay; no expanded search or integration."""

import argparse
import numpy as np
from mevo_mb.fast_sim.f7_navigation import MobilityFeedback
from .common import ROOT, BASE, js, write, identity, objects


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f7/protocol.json")
    env = proto["environments"][1]
    failure = js(BASE / "slow_v2/run/FAILURE_PROPOSAL.json")
    p, b, r, _, _ = objects()
    center = np.array(env["initial_centroid_m"])
    tensor = np.array(env["flow_gradient_s_inv"])
    u0 = np.array(env["flow_velocity_m_s"])
    def flow(x):
        return u0 + (x - center) @ tensor.T
    policy = dict(proto["f7_controller"], flow_enabled=True, analytic_pressure_feasibility_repair=True)
    ctl = MobilityFeedback(b, r, p, env["target_m"], flow, policy)
    ctl.previous = np.array(failure["controller"]["previous"])
    try:
        event = ctl.update(np.array(failure["positions_m"]), failure["time_s"])
        row = ctl.records[-1]
        assert row["min_pressure_Pa"] >= 100 - 1e-6
        assert max(event["amplitude_Pa"]) <= policy["amplitude_cap_Pa"]
        assert row["feasibility_repair"]["optimizer_iterations_added"] == 0
        write(
            out / "CHECKS.json",
            dict(
                status="PASS",
                identity=identity(),
                event=event,
                optimization=row,
                source="slow_v2 actual failed control observation; same optimizer warm start",
                policy=policy,
            ),
        )
    except RuntimeError as exc:
        write(
            out / "CHECKS.json",
            dict(status="FAIL", identity=identity(), reason=str(exc), failure=ctl.failed_update),
        )


if __name__ == "__main__":
    main()
