"""One 200-bubble/5-ms task and bounded CUDA-event component timings."""

import argparse
import json
from time import perf_counter
import h5py
import warp as wp
from mevo_mb.fast_sim.models import Physics, ResponseSection
from mevo_mb.fast_sim.engine import FastSimulator
from mevo_mb.fast_sim import kernels as k
from .common import BASE, inputs, write
from .experiment import run_case


def choose(studies):
    candidates = []
    for name in studies:
        d = json.loads((BASE / name / "RESULT.json").read_text())
        for row in d["rows"]:
            if row["numerical_quality"] == "PASS" and row.get("step_statistics_pass"):
                candidates.append(
                    dict(
                        study=name,
                        name=row["name"],
                        dt_s=row["nominal_dt_s"],
                        contact_iterations=d["contact_iterations"],
                        wall_s=row["wall_s"],
                    )
                )
    if not candidates:
        raise RuntimeError("No eligible step")
    return min(candidates, key=lambda r: r["wall_s"])


def event_measure(sim, fn, count=100):
    sim.reset()
    sim.control({"time_s": 0, "pressure_Pa": 1000})
    sim.fields()
    sim.pairs()
    wp.synchronize()
    with wp.ScopedCapture(device=sim.device) as capture:
        for _ in range(count):
            fn()
    start = wp.Event(device=sim.device, enable_timing=True)
    end = wp.Event(device=sim.device, enable_timing=True)
    wp.record_event(start)
    t = perf_counter()
    wp.capture_launch(capture.graph)
    enqueue = perf_counter() - t
    wp.record_event(end)
    gpu_ms = wp.get_event_elapsed_time(start, end)
    return dict(
        repetitions=count,
        gpu_mean_us=gpu_ms * 1000 / count,
        host_graph_enqueue_s=enqueue,
        scope="CUDA event elapsed within a captured repeated-kernel graph, compiled/warmed; no cross-component additivity assumed",
    )


def main(label, studies):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    selected = choose(studies)
    write(out / "SELECTED_STEP.json", selected)
    cfg = dict(
        count=200,
        seed=20261008,
        cloud_radius_m=50e-6,
        duration_s=0.005,
        dt_s=selected["dt_s"],
        sample_interval_s=0.0005,
        physics={"contact_iterations": selected["contact_iterations"]},
        events=[{"time_s": 0, "pressure_Pa": 1000}],
    )
    row = run_case(out / "N200_5ms", cfg)
    with h5py.File(out / "N200_5ms/trajectory.h5") as h:
        x = h["positions_m"][-1]
    p = Physics(**cfg["physics"])
    response = ResponseSection(inputs()[0], p)
    sim = FastSimulator(x, p, response)
    components = {}
    components["field_and_response"] = event_measure(sim, sim.fields)
    components["full_pair_force"] = event_measure(sim, sim.pairs)
    sim.tmp1.assign(x.reshape(-1, 3))
    components["one_contact_pass"] = event_measure(
        sim,
        lambda: sim.launch(
            k.contact,
            sim.total,
            [sim.tmp1, sim.tmp2, sim.n, p.radius_m, p.contact_relaxation, sim.diag, sim.error],
        ),
    )
    components["full_step_graph"] = event_measure(sim, lambda: sim.step(selected["dt_s"]))
    sim.reset()
    sim.control({"time_s": 0, "pressure_Pa": 1000})
    wp.synchronize()
    t = perf_counter()
    for _ in range(100):
        sim.step(selected["dt_s"])
    enqueue = perf_counter() - t
    wp.synchronize()
    raw_wall = perf_counter() - t
    t = perf_counter()
    for _ in range(100):
        wp.synchronize()
    empty_sync = (perf_counter() - t) / 100
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            selected=selected,
            short_run=row,
            components=components,
            raw_python_launch_100_steps_enqueue_s=enqueue,
            raw_python_launch_100_steps_wall_s=raw_wall,
            empty_synchronize_mean_s=empty_sync,
            short_target_10s="PASS" if row["numerical_quality"] == "PASS" and row["wall_s"] <= 10 else "FAIL",
            predicted_1s_N32_wall_s=selected["wall_s"] / 0.05,
            prediction_is_completion=False,
        ),
    )


if __name__ == "__main__":
    a = argparse.ArgumentParser()
    a.add_argument("--label", required=True)
    a.add_argument("--studies", nargs="+", required=True)
    q = a.parse_args()
    main(q.label, q.studies)
