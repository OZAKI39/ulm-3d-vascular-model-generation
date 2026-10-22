"""One final-backend N=200, 5 ms regression and a small cost decomposition."""

import argparse
import json
import h5py
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.f2_engine import F2Simulator
from mevo_mb.fast_sim import kernels as k
from scripts.fast_pzt.performance import event_measure
from .common import BASE, F1, write
from .experiment import run_case


def main(label, adoption):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    selected = json.loads((BASE / adoption / "SELECTED.json").read_text())
    cfg = dict(
        count=200,
        seed=20261008,
        cloud_radius_m=50e-6,
        duration_s=0.005,
        dt_s=selected["dt_s"],
        sample_interval_s=0.0005,
        physics=selected["physics"],
        events=[{"time_s": 0, "pressure_Pa": 1000}],
    )
    row = run_case(out / "N200_5ms", cfg)
    with h5py.File(out / "N200_5ms/trajectory.h5") as h:
        x = h["positions_m"][-1]
    with h5py.File(F1 / "balanced_performance_v1/N200_5ms/trajectory.h5") as h:
        with h5py.File(out / "N200_5ms/trajectory.h5") as f:
            import numpy as np

            np.testing.assert_array_equal(f["positions_m"][0], h["positions_m"][0])
    p = F2Physics(**selected["physics"])
    r = ResponseSection(BASE / "inputs/response.h5", p)
    sim = F2Simulator(x, p, r)
    components = {}
    components["field_response"] = event_measure(sim, sim.fields)
    components["pair_force"] = event_measure(sim, sim.pairs)
    components["one_local_mobility_pass"] = event_measure(
        sim,
        lambda: sim.launch(
            k.local_mobility,
            sim.total,
            [
                sim.x,
                sim.free,
                sim.free,
                sim.vel1,
                sim.n,
                p.radius_m,
                p.lubrication_activation,
                p.lubrication_floor,
            ],
        ),
    )
    sim.tmp1.assign(x.reshape(-1, 3))
    components["one_contact_pass"] = event_measure(
        sim,
        lambda: sim.launch(
            k.contact,
            sim.total,
            [sim.tmp1, sim.tmp2, sim.n, p.radius_m, p.contact_relaxation, sim.diag, sim.error],
        ),
    )
    components["full_step"] = event_measure(sim, lambda: sim.step(selected["dt_s"]))
    f1 = json.loads((F1 / "balanced_performance_v1/RESULT.json").read_text())
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            selected=selected,
            run=row,
            components=components,
            F1_simulation_wall_s=f1["short_run"]["wall_s"],
            wall_ratio_to_F1=row["wall_s"] / f1["short_run"]["wall_s"],
            macro_step_ratio=row["macro_step_count"] / f1["short_run"]["macro_step_count"],
            same_initial_positions=True,
            status="PASS_10S_TARGET" if row["numerical_quality"] == "PASS" and row["wall_s"] < 10 else "FAIL",
            runtime_scope="Simulation wall includes sample synchronization; cold compilation and setup reported separately; N=200,5ms only",
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--adoption", required=True)
    a = ap.parse_args()
    main(a.label, a.adoption)
