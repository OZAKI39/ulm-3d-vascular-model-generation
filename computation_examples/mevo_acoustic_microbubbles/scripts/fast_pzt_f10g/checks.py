"""Supplemental adapter tests: two-bubble dual oracle, control switches and raw rejection."""

import argparse
import gzip
import json
import numpy as np
from mevo_mb.fast_sim.f10g_diagnostic import SoftContactDiagnostic
from mevo_mb.fast_sim.f10g_contact import resistance, enumerated_reference
from mevo_mb.fast_sim.flow_gpu import AffineFlow
from scripts.fast_pzt_f9.common import objects
from .common import BASE, ROOT, CASES, state, js, write, identity, clean


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    p, b, r, g, _ = objects()
    k = 2 * js(ROOT / "configs/fast_pzt/f10g/stiffness.json")["kc_N_m"]
    center = np.array(state(CASES[0])[1]["initial_positions_m"]).mean(0)
    cfg = js(ROOT / "results/fast_pzt/f10e/contact_v1/O_1000_1102_C0/CONFIG.json")
    ev = cfg["events"][0] | {"time_s": 0, "flow_enabled": True}
    events = [ev, ev | {"time_s": 2e-5, "amplitude_Pa": [0] * 4}, ev | {"time_s": 4e-5}]
    x = np.array([center + [-p.radius_m - 1e-9, 0, 0], center + [p.radius_m + 1e-9, 0, 0]])
    s = SoftContactDiagnostic(x, p, r, b, vessel=g, flow=AffineFlow(np.array([0.0002, 0, 0])), kc_N_m=k)
    s.configure_schedule(events, 1e-5)
    s.control(events[0])
    rows = []
    for step in range(6):
        before = s.x.numpy().copy()
        s.advance(1e-5, 1, graphs=False)
        free = np.array(s.trace["free_velocity_m_s"])
        A = resistance(before, p.radius_m)["A"]
        ref = enumerated_reference(
            before, free, p.radius_m, 6 * np.pi * p.viscosity_Pa_s * p.radius_m, 1e-5, k, A
        )
        err = float(np.max(np.abs(s.x.numpy() - ref["positions_m"])))
        assert err < 1e-12
        if step in [2, 3]:
            assert np.all(s.trace["total_force_N"] == 0)
        rows.append(
            dict(step=step + 1, reference_position_error_m=err, reference=ref, step_evidence=s.soft_steps[-1])
        )
    np.testing.assert_allclose(
        [v["time_s"] for v in s.executed_controls], [0, 2e-5, 4e-5], rtol=0, atol=1e-15
    )
    # Save complete deliberately rejected state, candidate, quality and clock.
    x = np.array([center + [-p.radius_m + 200e-9, 0, 0], center + [p.radius_m - 200e-9, 0, 0]])
    s = SoftContactDiagnostic(x, p, r, b, vessel=g, flow=AffineFlow(np.zeros(3)), kc_N_m=k)
    s.configure_schedule([ev], 1e-5)
    s.control(ev)
    try:
        s.advance(1e-5, 1, graphs=False)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Invalid overlap accepted")
    assert s.first_rejection and s.accepted_steps == 0 and np.array_equal(s.x.numpy(), x)
    try:
        s.advance(1e-5, 1, graphs=False)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Continuation after failure")
    with gzip.open(out / "ENGINEERING_EVIDENCE.json.gz", "wt") as stream:
        json.dump(
            clean(
                dict(
                    pair_steps=rows,
                    rejected_step=s.soft_steps[-1],
                    trace=s.trace,
                    first_rejection=s.first_rejection,
                    accepted_clock=int(s.step_count.numpy()[0]),
                )
            ),
            stream,
            separators=(",", ":"),
        )
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            two_bubble_adapter_steps=6,
            deliberate_rejections=1,
            failure_continuation_blocked=True,
            max_independent_reference_error_m=max(v["reference_position_error_m"] for v in rows),
            control_switches=[0, 2e-5, 4e-5],
            training_ready=False,
        ),
    )
    print("F10G diagnostic adapter checks PASS")


if __name__ == "__main__":
    main()
