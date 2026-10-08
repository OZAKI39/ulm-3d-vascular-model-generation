"""One matched N200/5ms regression per execution mode, plus bounded GPU timings."""

import argparse
import h5py
import numpy as np
from mevo_mb.fast_sim.f3_engine import F3Simulator
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.models import ResponseSection
from scripts.fast_pzt.performance import event_measure
from .common import BASE, F2, write, js, identity
from .experiment import run_case


def main(label, repeats):
    assert js(BASE / repeats / "RESULT.json")["status"] == "PASS_SAME_DEVICE_BITWISE"
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    old = js(F2 / "performance_v1/RESULT.json")
    cfg = old["run"]["configuration"]
    rows = {}
    timings = {}
    with h5py.File(F2 / "performance_v1/N200_5ms/trajectory.h5") as h:
        initial = h["positions_m"][0]
    for mode in ["fast_atomic", "stable_reduction"]:
        rows[mode] = run_case(out / mode, dict(cfg, execution_mode=mode))
        with h5py.File(out / mode / "trajectory.h5") as h:
            np.testing.assert_array_equal(h["positions_m"][0], initial)
            x = h["positions_m"][-1]
        sim = F3Simulator(
            x, F2Physics(), ResponseSection(BASE / "inputs/response.h5", F2Physics()), execution_mode=mode
        )
        timings[mode] = dict(
            pair_force_stage=event_measure(sim, sim.pairs),
            full_step=event_measure(sim, lambda: sim.step(1e-5)),
        )
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            identity=identity(),
            rows=rows,
            GPU_timings=timings,
            F2_wall_s=old["run"]["wall_s"],
            stable_to_atomic_wall_ratio=rows["stable_reduction"]["wall_s"] / rows["fast_atomic"]["wall_s"],
            same_initial_positions_as_F2=True,
            scope="One N200/5ms trajectory per execution mode at 10 us, same GPU/input specs; CUDA-event timings compiled/warmed and not additive; no new long/batch campaign",
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--repeats", required=True)
    a = ap.parse_args()
    main(a.label, a.repeats)
