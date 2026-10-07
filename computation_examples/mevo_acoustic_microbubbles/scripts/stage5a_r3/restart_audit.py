"""Exact full-system restart and accepted-cache rollback; bounded local only."""

import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, write
from scripts.stage5a_r3.run import run, restore


def main():
    args = dict(
        case="B2",
        duration=20e-6,
        max_dt=5e-6,
        tol_factor=0.25,
        initial_dt=50e-9,
        budget=300,
        nonlinear_strategy="quasi_newton",
        nonsmooth_recovery=True,
        fixture=BASE / "diagnosis/reproduction_v1/nominal/terminal.npz",
        optimized=True,
    )
    full, a = run(**args, label="restart_v2_uninterrupted")
    part, b = run(**args, label="restart_v2_intentional_stop", stop_after_macros=3)
    final, c = run(**args, label="restart_v2_resumed", resume=part / "CHECKPOINT.npz")

    def read(path):
        with h5py.File(path / "trajectories.h5") as h:
            return h["time_s"][:], h["positions_m"][:]

    ta, xa = read(full)
    tb, xb = read(part)
    tc, xc = read(final)
    merged = np.concatenate([xb, xc[1:]])
    times = np.r_[tb, tc[1:]]
    same = np.array_equal(xa, merged, equal_nan=True) and np.array_equal(ta, times)
    sa, ma, ca = restore(full / "CHECKPOINT.npz")
    sc, mc, cc = restore(final / "CHECKPOINT.npz")
    same_cache = (
        ca.particle_ids == cc.particle_ids
        and ca.dual == cc.dual
        and ca.active == cc.active
        and all(np.array_equal(ca.primal[k], cc.primal[k]) for k in ca.primal)
    )
    report = dict(
        status="PASS" if same and same_cache else "FAIL",
        bitwise_actual_positions_and_time=same,
        identical_accepted_cache=same_cache,
        full=str(full),
        part=str(part),
        resumed=str(final),
        original_start=False,
        scope="Exact continuation of unchanged R3 fixture/controller, no output clock truncation introduced by interruption",
    )
    write(BASE / "diagnosis/RESTART_V2.json", report)
    print(report, flush=True)
    assert report["status"] == "PASS"


if __name__ == "__main__":
    main()
