"""Replay accepted RK steps from a stable diagnostic copy; never resume a run."""
import argparse
import copy
import json
import time

import h5py
import numpy as np

from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.coupling.bifurcation_experiment import BifurcationExperimentTracker
from scripts.stage4.common import ROOT, new_run, require, finish
from scripts.stage4.experiments import setup
from scripts.stage4.geometry_certificate import install


def replay(config, field, flow, positions, radii, start_time, certificate, per_particle=False):
    tracker = BifurcationExperimentTracker(copy.deepcopy(config), field, flow)
    counts = install(tracker.cylinder, per_particle=per_particle) if certificate else None
    x, t = positions.copy(), start_time
    nominal = config["particles"]["dt_s"]
    start = time.perf_counter()
    state = tracker.evaluate(x, radii, t)
    initial = {key: value.copy() for key, value in state.items()}
    rows = []
    for _ in range(64):
        target = (np.floor(t/nominal+1e-10)+1)*nominal
        speed = np.linalg.norm(state["total_velocity_m_s"], axis=1)
        dt = min(target-t, float(np.min(config["particles"]["max_displacement_per_radius"]*radii/np.maximum(speed, 1e-30))))
        x, wall, contacts, correction, iterations = tracker.step(x, radii, t, dt, state)
        t = min(t+dt, target)
        state = tracker.evaluate(x, radii, t)
        tracker.monitor.record_accepted_interval(dt)
        rows.append(dict(x=x.copy(), wall=wall.copy(), contacts=contacts, correction=correction,
                         iterations=iterations, dt=dt, t=t, state={k: v.copy() for k, v in state.items()}))
    return dict(initial=initial, rows=rows, runtime_s=time.perf_counter()-start,
                certificate=counts, timings=tracker.timings, validity=tracker.monitor.summary())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-particle", action="store_true")
    args = parser.parse_args()
    snapshots = require("performance/secondary_contact_snapshot")
    _, _, _, flow, _, _ = setup()
    out = new_run("performance/later_particle_geometry_certificate" if args.per_particle else "performance/later_geometry_certificate",
                  dict(accepted_steps=64, per_particle=args.per_particle, production_deployed=False), ["performance/secondary_contact_snapshot"])
    results = []
    for source in snapshots["rows"]:
        run = (ROOT/source["source"]).parents[1]
        config = json.loads((run/"RUN_IDENTITY.json").read_text())["config"]
        field = ComplexAcousticField.load(ROOT/config["effective_acoustic_field"]["path"])
        with h5py.File(ROOT/source["snapshot"]) as h:
            x, r, t = h["positions_m"][-1], h["radii_m"][:], float(h["time_s"][-1])
            saved = {k: h[k][-1] for k in ["flow_velocity_m_s", "free_total_velocity_m_s", "total_velocity_m_s"]}
        baseline = replay(config, field, flow, x, r, t, False)
        candidate = replay(config, field, flow, x, r, t, True, args.per_particle)
        initial_equal = all(np.array_equal(baseline["initial"][k], v) for k, v in saved.items())
        equality = []
        for a, b in zip(baseline["rows"], candidate["rows"], strict=True):
            equality.append(bool(all(np.array_equal(a[k], b[k]) for k in ["x", "wall"])
                and all(a[k] == b[k] for k in ["contacts", "correction", "iterations", "dt", "t"])
                and all(np.array_equal(a["state"][k], b["state"][k], equal_nan=True) for k in a["state"])))
        validity_equal = baseline["validity"] == candidate["validity"]
        results.append(dict(source_snapshot=source, replay_initial_velocity_bitwise_matches_saved=initial_equal,
                            accepted_step_bitwise_equal=equality, validity_equal=validity_equal,
                            baseline_runtime_s=baseline["runtime_s"], candidate_runtime_s=candidate["runtime_s"],
                            start_time_s=t, end_time_s=baseline["rows"][-1]["t"], certificate=candidate["certificate"]))
    passed = all(row["replay_initial_velocity_bitwise_matches_saved"] and all(row["accepted_step_bitwise_equal"])
                 and row["validity_equal"] for row in results)
    finish(out, dict(status="PASS" if passed else "FAIL", rows=results, production_deployed=False,
                     scope="64-step deterministic replays from copied partial frames, not resumed production or completed branch outcomes. Original contact NNLS retained."))


if __name__ == "__main__":
    main()
