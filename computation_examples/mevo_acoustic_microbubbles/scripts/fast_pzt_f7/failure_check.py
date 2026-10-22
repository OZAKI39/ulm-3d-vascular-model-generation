"""Inject one observer failure; verify true prefix and unmatched-diagnostic rollback."""

import argparse
import numpy as np
from mevo_mb.fast_sim.f7_navigation import RecordedFeedbackSimulator
from .common import ROOT, BASE, objects, js, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f7/protocol.json")
    env = proto["environments"][0]
    p, b, r, g, _ = objects()

    class Injected:
        policy = dict(interval_s=0.0005)
        records = []

        def update(self, x, t):
            if t > 0:
                raise RuntimeError("INJECTED_OBSERVER_FAILURE")
            return dict(proto["frozen_F6_events"][0])

    sim = RecordedFeedbackSimulator(
        np.array(env["positions_m"]),
        p,
        r,
        b,
        vessel=g,
        vessel_enabled=False,
        flow_enabled=False,
        mobility_mode="rpy_free_space",
        controller=Injected(),
    )
    result = sim.run(0.001, 1e-5, [proto["frozen_F6_events"][0]], sample_interval_s=0.0005)
    assert result["summary"]["status"] == "FAILED"
    assert result["summary"]["reason"] == "INJECTED_OBSERVER_FAILURE"
    assert len(result["time_s"]) == len(result["quality_samples"]) == len(result["flow"]["velocity_m_s"]) == 1
    assert abs(result["failure_proposal"]["time_s"] - 0.0005) < 1e-12
    assert not np.array_equal(result["failure_proposal"]["positions_m"], np.array(env["positions_m"]))
    write(
        out / "CHECKS.json",
        dict(
            status="PASS",
            identity=identity(),
            saved_times=result["time_s"],
            failure=result["failure_proposal"],
            meaning="Synthetic controller exception after real500us integration; no scientific failure suppressed or quality budget changed.",
        ),
    )
    print("Failure-prefix preservation PASS", flush=True)


if __name__ == "__main__":
    main()
