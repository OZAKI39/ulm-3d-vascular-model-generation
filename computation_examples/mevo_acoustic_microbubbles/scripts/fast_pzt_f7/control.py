"""Admit cached GPU predictor, then bounded RPY-aware feedback and dt check."""

import argparse
import time
import numpy as np
from mevo_mb.fast_sim.f7_predictor import CollectivePredictor
from .common import ROOT, BASE, write, js, identity, objects
from .fixed import probe
from .runner import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--phase", choices=["check", "pilot", "formal", "half", "slow"], required=True)
    ap.add_argument("--dt", type=float, default=1e-5)
    ap.add_argument("--repair", action="store_true")
    a = ap.parse_args()
    assert js(BASE / "replay_v1/REPLAY.json")["status"] == "PASS"
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f7/protocol.json")
    obj = objects()
    p, b, r, _, _ = obj
    proto["feedback_policy"] = proto["f7_controller"]
    if a.repair:
        assert a.phase == "slow"
        assert js(BASE / "feasibility_v1/CHECKS.json")["status"] == "PASS"
        proto["feedback_policy"]["analytic_pressure_feasibility_repair"] = True
    write(out / "REGISTERED_PROTOCOL.json", dict(proto, phase=a.phase, dt_s=a.dt))
    write(out / "IDENTITY.json", identity())
    start = time.perf_counter()
    if a.phase == "check":
        refs = js(BASE / "fixed_v1/CHECKS.json")["states"]
        data = np.load(BASE / "fixed_v1/fixed_states.npz")
        predictor = CollectivePredictor(b, r, p)
        rows = {}
        for name, ref in refs.items():
            x = data[name + "_positions_m"]
            uf = data[name + "_M2_flow_m_s"]
            ev = ref["event"]
            predictor.prepare(x, uf)
            q = predictor.evaluate(ev["amplitude_Pa"], ev["phase_rad"], ev["enabled"], details=True)
            for key in ["primary_N", "secondary_N", "velocity_m_s"]:
                np.testing.assert_allclose(
                    q[key], data[name + "_M2_" + key], rtol=1e-9, atol=1e-23 if key.endswith("_N") else 1e-16
                )
            # Changed amplitude/phase candidate must match full production sampling too.
            alt = dict(
                ev,
                amplitude_Pa=(np.array(ev["amplitude_Pa"]) * 0.7).tolist(),
                phase_rad=(np.array(ev["phase_rad"]) + np.array([0.0, 0.1, -0.1, 0.05])).tolist(),
            )
            actual = probe(x, alt, obj, "M2", name != "mature", name != "mature")
            test = predictor.evaluate(alt["amplitude_Pa"], alt["phase_rad"], details=True)
            np.testing.assert_allclose(test["velocity_m_s"], actual["velocity_m_s"], rtol=1e-9, atol=1e-16)
            rows[name] = dict(
                max_velocity_difference_m_s=np.abs(test["velocity_m_s"] - actual["velocity_m_s"]).max()
            )
        write(
            out / "CHECKS.json",
            dict(
                status="PASS",
                rows=rows,
                wall_s=time.perf_counter() - start,
                scope="Cached linear basis + table + primary + F2 secondary + full RPY agrees with production at three real states and two controls each. No PBD prediction.",
            ),
        )
        return
    assert js(BASE / "predictor_v1/CHECKS.json")["status"] == "PASS"
    if a.phase in ["formal", "half", "slow"]:
        assert js(BASE / "controller_pilot_v2/run/SUMMARY.json")["numerical_quality"] == "PASS"
    if a.phase == "slow":
        main = js(BASE / "controller_v2/run/SUMMARY.json")
        assert main["numerical_quality"] == "PASS" and main["navigation"]["arrived_and_held"]
    env = proto["environments"][1 if a.phase == "slow" else 0]
    s = run_case(
        out / "run",
        env,
        proto["frozen_F6_events"][0],
        proto,
        obj,
        feedback="rpy",
        mobility="rpy_free_space",
        duration=0.03 if a.phase == "pilot" else 1.0,
        dt=a.dt,
    )
    write(
        out / "CAMPAIGN.json",
        dict(status="COMPLETE", phase=a.phase, runs={"run": s}, process_wall_s=time.perf_counter() - start),
    )


if __name__ == "__main__":
    main()
