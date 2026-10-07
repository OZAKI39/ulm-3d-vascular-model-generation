"""Changing actual-state substeps; component workers and force device measured."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from time import perf_counter, process_time
import numpy as np
from scripts.stage5a_r3.common import BASE, load, write, identity
from scripts.stage5a_r3.run import restore
from mevo_mb.bubbles.motion_r3 import MotionR3, SemiImplicitR3
from mevo_mb.bubbles.motion_backend_r3 import MotionR3Optimized
from mevo_mb.bubbles.motion_sweep_r3 import WallSweepR3


def benchmark(task):
    case, window, layout, optimized, phase, *requested_dt = task
    tracker, sweep0, cfg, state, inputs = load(case)
    fixture = None
    if window == "event":
        fixture = BASE / "diagnosis/reproduction_v1/nominal/terminal.npz"
    elif window == "late":
        fixture = BASE / "fixtures/B3_actual_3.00490000.npz"
    elif window == "stable":
        fixture = BASE / "fixtures/B2_actual_3.00300000.npz"
    if fixture:
        state, meta, cache = restore(fixture)
    x = state["positions_m"][state["active"]]
    r = state["radii_m"][state["active"]]
    ids = state["particle_id"][state["active"]]
    t = float(state["time_s"])
    workers = {"serial": 1, "parallel2": 2, "parallel4": 4, "global": 1}[layout]
    cls = MotionR3Optimized if optimized else MotionR3
    backend = cls(cfg["medium"]["viscosity_Pa_s"], components=layout != "global", component_workers=workers)
    sweep = WallSweepR3(tracker.cylinder, sweep0.path)
    step = SemiImplicitR3(
        backend, tracker.cylinder, tracker.evaluate, sweep, nonlinear_strategy="quasi_newton"
    )
    dt = requested_dt[0] if requested_dt else 2e-9 if window == "initial" else 5e-8
    out = BASE / "benchmarks/layout_v2" / phase
    out.mkdir(parents=True, exist_ok=True)
    label = f"{case}_{window}_{layout}_{'optimized' if optimized else 'v1'}"
    rows, allx = [], [x.copy()]
    for repeat in range(5):
        begin, cpu = perf_counter(), process_time()
        try:
            y, u, cache, diag = step.advance(x, r, ids, t, dt, "changing_state")
            status, error = "PASS", None
        except Exception as exc:
            status, error = "FAIL", repr(exc)
            diag = getattr(exc, "diagnostics", {})
        wall, cpu_cost = perf_counter() - begin, process_time() - cpu
        row = dict(
            case=case,
            window=window,
            layout=layout,
            optimized=optimized,
            repeat=repeat,
            status=status,
            error=error,
            wall_s=wall,
            process_cpu_s=cpu_cost,
            dt_s=dt,
            time_s=t,
            diagnostics=diag,
            persistent_executor=backend.executor is not None,
        )
        rows.append(row)
        if status != "PASS":
            break
        backend.cache = cache
        x, t = y, t + dt
        allx.append(x.copy())
    if backend.executor:
        backend.executor.shutdown()
    np.savez_compressed(out / (label + ".npz"), positions_m=np.array(allx), particle_id=ids, radii_m=r)
    write(out / (label + ".json"), rows)
    return rows


def gpu():
    from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
    import cupy as cp

    rows = []
    for case in ["B2", "B3"]:
        tr, _, cfg, state, _ = load(case)
        mask = state["active"]
        x, r = state["positions_m"][mask], state["radii_m"][mask]
        p = tr.field.interpolate_complex_pressure(x)
        vol = tr.response.harmonic_response(r, tr.field.frequency_Hz, p)["harmonics_m3"]
        cpu = ClassicalFarFieldBjerknes("numpy", enforce_far_field=False)
        gpu = ClassicalFarFieldBjerknes("cupy", enforce_far_field=False)
        f, rho = tr.field.frequency_Hz, cfg["medium"]["density_kg_m3"]
        for repeat in range(11):
            begin = perf_counter()
            ref = cpu.forces(x, r, vol, f, rho)
            cpu_s = perf_counter() - begin
            begin = perf_counter()
            answer = gpu.forces(x, r, vol, f, rho).get()
            gpu_s = perf_counter() - begin
            begin = perf_counter()
            dev = [cp.asarray(v) for v in [x, r, vol]]
            cp.cuda.Stream.null.synchronize()
            upload = perf_counter() - begin
            e0, e1 = cp.cuda.Event(), cp.cuda.Event()
            begin = perf_counter()
            e0.record()
            result = gpu.forces(*dev, f, rho)
            e1.record()
            kernel_api = perf_counter() - begin
            begin = perf_counter()
            e1.synchronize()
            sync_s = perf_counter() - begin
            kernel_event = cp.cuda.get_elapsed_time(e0, e1) / 1000
            begin = perf_counter()
            result.get()
            download = perf_counter() - begin
            np.testing.assert_allclose(answer, ref, rtol=2e-12, atol=1e-26)
            rows.append(
                dict(
                    case=case,
                    repeat=repeat,
                    CPU_s=cpu_s,
                    GPU_end_to_end_s=gpu_s,
                    H2D_including_sync_s=upload,
                    GPU_API_s=kernel_api,
                    remaining_sync_s=sync_s,
                    CUDA_event_inclusive_s=kernel_event,
                    D2H_s=download,
                    max_force_error_N=float(np.max(abs(answer - ref))),
                )
            )
    return rows


def main(isolated=False):
    tasks = [
        (c, w, layout, o, "isolated" if isolated else "development")
        for c, w in [("B2", "initial"), ("B2", "event"), ("B3", "initial"), ("B3", "late")]
        for layout, o in [
            ("serial", False),
            ("serial", True),
            ("parallel2", True),
            ("parallel4", True),
            ("global", True),
        ]
    ]
    begin = perf_counter()
    # Each single-trajectory layout is measured sequentially. If unrelated
    # development probes are active, label it diagnostic and repeat when isolated.
    rows = [row for task in tasks for row in benchmark(task)]
    gpu_rows = gpu()
    output = "LAYOUT_ISOLATED.json" if isolated else "LAYOUT_DEVELOPMENT.json"
    write(
        BASE / "benchmarks" / output,
        dict(
            rows=rows,
            gpu=gpu_rows,
            identity=identity(),
            wall_s=perf_counter() - begin,
            resource_scope="isolated"
            if isolated
            else "development; other declared diagnostic workers may compete",
            GIL_scope="Native solver calls and SciPy release behavior measured via per-component thread CPU and wall time; no assumption of Python parallelism",
        ),
    )
    # Throughput experiment uses two independent jobs with single-worker QPs.
    subset = [("B2", "event", "serial", True), ("B3", "late", "serial", True)]
    begin = perf_counter()
    serial = [benchmark((*t, "throughput_serial")) for t in subset]
    serial_s = perf_counter() - begin
    begin = perf_counter()
    with ProcessPoolExecutor(max_workers=2) as pool:
        parallel = list(pool.map(benchmark, [(*t, "throughput_parallel") for t in subset]))
    parallel_s = perf_counter() - begin
    write(
        BASE / "benchmarks/INDEPENDENT_THROUGHPUT.json",
        dict(
            serial_wall_s=serial_s,
            parallel_wall_s=parallel_s,
            speedup=serial_s / parallel_s,
            workers=2,
            label="Independent task throughput; NOT single-trajectory speedup",
            serial=serial,
            parallel=parallel,
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--isolated", action="store_true")
    main(p.parse_args().isolated)
