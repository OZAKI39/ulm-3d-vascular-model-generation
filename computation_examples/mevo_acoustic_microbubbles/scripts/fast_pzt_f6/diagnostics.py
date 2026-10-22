"""Minimal halved-step diagnosis triggered by failed 10 us pilot budgets."""

import argparse
from .common import ROOT, BASE, write, js, identity, objects
from .science import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--refine", action="store_true")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f6/protocol.json")
    obj = objects()
    env0, env1, _ = proto["environments"]
    rows = {}
    if a.refine:
        for name, env, feedback in [("phase", env0, False), ("slow_feedback", env1, True)]:
            rows[name] = run_case(
                out / name,
                env,
                env["events"]["target_phase"],
                proto,
                obj,
                dt=2.5e-6,
                duration=0.03,
                feedback=feedback,
            )
        write(
            out / "DIAGNOSTICS.json",
            dict(
                rows=rows,
                identity=identity(),
                reason="5 us pilot also exceeded 0.01R; final bounded 2.5 us diagnostic, unchanged control and thresholds",
            ),
        )
        return
    # One target-phase test; the affine version differs only by the small stated flow.
    rows["phase_dt5"] = run_case(
        out / "phase_dt5", env0, env0["events"]["target_phase"], proto, obj, dt=5e-6, duration=0.005
    )
    rows["slow_feedback_dt5"] = run_case(
        out / "slow_feedback_dt5",
        env1,
        env1["events"]["target_phase"],
        proto,
        obj,
        dt=5e-6,
        duration=0.03,
        feedback=True,
    )
    rows["zero_feedback_event"] = run_case(
        out / "zero_feedback_event",
        env0,
        env0["events"]["target_phase"],
        proto,
        obj,
        duration=0.03,
        feedback=True,
    )
    rows["zero_RPY"] = run_case(
        out / "zero_RPY",
        env0,
        env0["events"]["target_phase"],
        proto,
        obj,
        duration=0.03,
        feedback=True,
        mobility="rpy_free_space",
    )
    write(
        out / "DIAGNOSTICS.json",
        dict(
            rows=rows,
            identity=identity(),
            reason="10 us target-phase and slow-affine feedback exceeded unchanged 0.01R contact budget in short_v2; no target/control/threshold change",
        ),
    )


if __name__ == "__main__":
    main()
