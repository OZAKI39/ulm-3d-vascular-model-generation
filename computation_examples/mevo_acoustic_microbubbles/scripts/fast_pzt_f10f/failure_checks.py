"""Targeted odd-round graph rejection and non-vacuous edge-time checks."""

import argparse
from dataclasses import replace
import h5py
import numpy as np
from mevo_mb.fast_sim.f10f_smallsteps import SmallStepsBatch, MODES
from mevo_mb.fast_sim.flow_gpu import AffineFlow
from .common import BASE, HISTORY, ROOT, objects, identity, js, write


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    p, b, r, g, f = objects()
    chosen = js(HISTORY / "SELECTION.json")["branch"]
    cfg = js(ROOT / chosen["relative_path"] / "CONFIG.json")
    x = np.array(cfg["initial_positions_m"])
    failures = []
    edges = []
    for mode in ["SS2", "SS4"]:
        dt, it = MODES[mode]
        pp = replace(p, contact_iterations=it)
        s = SmallStepsBatch(x, pp, r, b, flow=f, vessel=g, numerical_mode=mode, diagnostic=False)
        s.configure_schedule(cfg["events"], dt)
        s.control(cfg["events"][0])
        try:
            s.advance(dt, 20, graphs=True)
        except RuntimeError:
            pass
        else:
            raise AssertionError("Archived first-substep residual must fail")
        assert s.physical_steps == 1 and s.accepted_steps == 0
        with h5py.File(BASE / "short_v1" / ("O_1000_1102_" + mode) / "trajectory.h5") as h:
            np.testing.assert_array_equal(s.final_proposal.numpy(), h["rejected_interval/proposal_m"][0])
            np.testing.assert_array_equal(s.x.numpy(), h["last_legal_state/positions_m"][0])
        failures.append(
            dict(
                mode=mode,
                first_rejection=s.first_rejection,
                graph_rejected_proposal_equals_diagnostic=True,
                accepted_steps=0,
            )
        )
    center = x[0].mean(0)
    for mode, (dt, it) in MODES.items():
        pair = np.repeat(center[None], 2, axis=0)
        pair[:, 0] += np.array([-1, 1]) * (2.1 * p.radius_m + 0.2e-9) / 2
        ev = cfg["events"][0] | {"amplitude_Pa": (np.array(cfg["events"][0]["amplitude_Pa"]) * 0.1).tolist()}
        pp = replace(p, contact_iterations=it)
        s = SmallStepsBatch(
            pair, pp, r, b, flow=AffineFlow([0, 0, 0]), vessel=g, numerical_mode=mode, diagnostic=False
        )
        s.configure_schedule([ev], dt)
        s.control(ev)
        s.advance(dt, 8, graphs=True)
        count = int(s.event_count.numpy()[0])
        assert count >= 1, "Edge-time check must exercise an actual crossing"
        times = s.event_time.numpy()[:count]
        ids = s.event_ids.numpy()[:count]
        assert np.all(np.abs(times / dt - np.round(times / dt)) < 1e-10)
        assert np.all(times > 0) and np.all(times <= 8 * dt + 1e-14)
        np.testing.assert_array_equal(ids[:, 1:4], [[0, 1, 1]])
        edges.append(dict(mode=mode, ids=ids, times_s=times, dt_s=dt))
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            graph_first_error=failures,
            actual_edge_crossings=edges,
            scope="SOFTWARE_REGRESSION_ONLY; no additional independent scientific samples",
        ),
    )
    print("PASS failure and event checks", flush=True)


if __name__ == "__main__":
    main()
