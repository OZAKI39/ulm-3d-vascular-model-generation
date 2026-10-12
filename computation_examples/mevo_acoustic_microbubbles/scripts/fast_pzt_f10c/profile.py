"""Bounded existing-engine profiling; never used as throughput timing."""

import argparse
from collections import defaultdict
from time import perf_counter
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f6_engine import F6Simulator
from .common import BASE, objects, case, identity, write


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    begin = perf_counter()
    p, b, r, g, f = objects()
    static = perf_counter() - begin
    cfg = case("O_1000")
    x = np.array(cfg["environment"]["positions_m"])
    ev = cfg["events"][0]
    begin = perf_counter()
    s = F6Simulator(x, p, r, b, flow=f, vessel=g, graph_chunk=5)
    s.control(ev)
    model = perf_counter() - begin
    begin = perf_counter()
    s.graph(1e-5, 5)
    wp.synchronize()
    capture = perf_counter() - begin
    np.testing.assert_array_equal(s.state()[0], x)
    rows = []
    states = []
    for graph in [True, False]:
        s.reset()
        s.control(ev)
        wp.synchronize()
        start = wp.Event(enable_timing=True)
        end = wp.Event(enable_timing=True)
        begin = perf_counter()
        wp.record_event(start)
        for _ in range(20):
            s.advance(1e-5, 5, graphs=graph)
        wp.record_event(end)
        wp.synchronize()
        wall = perf_counter() - begin
        event = wp.get_event_elapsed_time(start, end) / 1000
        states.append(s.state().copy())
        rows.append(
            dict(
                graph_enabled=graph,
                physical_steps=100,
                wall_s=wall,
                cuda_event_span_s=event,
                quality=s.quality(),
            )
        )
    np.testing.assert_array_equal(*states)
    s.reset()
    s.control(ev)
    wp.synchronize()
    wp.timing_begin(cuda_filter=wp.TIMING_ALL)
    begin = perf_counter()
    for _ in range(20):
        s.step(1e-5)
    activities = wp.timing_end()
    instrumented = perf_counter() - begin
    grouped = defaultdict(lambda: dict(count=0, event_span_ms=0.0))
    raw = []
    for v in activities:
        grouped[v.name]["count"] += 1
        grouped[v.name]["event_span_ms"] += v.elapsed
        raw.append(dict(name=v.name, filter=v.filter, elapsed_ms=v.elapsed))
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            static_prepare_s=static,
            model_prepare_s=model,
            module_load_compile_observed_s=s.compile_s,
            graph_capture_s=capture,
            comparison=rows,
            graph_bitwise_equal=True,
            profile_steps=20,
            instrumented_wall_s=instrumented,
            activities=dict(grouped),
            raw=raw,
            pure_kernel_gpu_time_s=None,
            pure_kernel_measurement="NOT_MEASURED: no CUPTI/Nsight installed; Warp inserts events per API activity. These spans can contain launch latency and instrumentation overhead.",
            scope="Existing F6 M0, no algorithm changes; bounded 1 ms comparison and 0.2 ms profiling; output not included; warm disk compilation cache retained.",
            documentation="https://nvidia.github.io/warp/latest/user_guide/execution_and_performance/profiling.html",
        ),
    )
    print(dict(static=static, model=model, capture=capture, rows=rows, activities=dict(grouped)), flush=True)


if __name__ == "__main__":
    main()
