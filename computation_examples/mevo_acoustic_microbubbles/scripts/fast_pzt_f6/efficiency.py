"""Bounded event/adaptive and two-sample throughput checks, no stress campaign."""

import argparse
import time
import h5py
import numpy as np
from mevo_mb.fast_sim.f6_engine import F6Simulator
from mevo_mb.fast_sim.flow_gpu import AffineFlow
from mevo_mb.fast_sim.storage import save_run
from scripts.fast_pzt_f2.common import quality_pass
from .common import BASE, F5, write, js, identity, objects


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    p, basis, response, geometry, realflow = objects()
    # Exact affine advection and noncommensurate switch; imposed flow is a test.
    initial = np.array([[1.1e-6, 80e-6, 40e-6]])
    flow = AffineFlow([1e-4, 0, 0])
    on = dict(
        time_s=0.0,
        amplitude_Pa=[0.0] * 4,
        phase_rad=[0.0] * 4,
        enabled=[False] * 4,
        flow_enabled=True,
        secondary_enabled=False,
    )
    off = dict(on, time_s=17e-6, flow_enabled=False)
    sim = F6Simulator(initial, p, response, basis, flow=flow, acoustic_enabled=False, vessel_enabled=False)
    result = sim.run(48e-6, 1e-5, [on, off], 1e-5, adaptive=True)
    assert result["summary"]["status"] == "COMPLETE"
    np.testing.assert_allclose(
        result["positions_m"][-1, 0], initial + np.array([17e-10, 0, 0]), atol=1e-18, rtol=0
    )
    times = np.cumsum(result["physical_dt_s"])
    assert np.any(np.isclose(times, 17e-6, atol=1e-15, rtol=0))
    np.testing.assert_allclose(times[-1], 48e-6, rtol=1e-12)
    np.savez_compressed(
        out / "ADAPTIVE_EVENT.npz",
        time_s=result["time_s"],
        positions_m=result["positions_m"],
        physical_dt_s=result["physical_dt_s"],
        attempts=result["adaptive_attempts"],
    )
    write(out / "ADAPTIVE_EVENT.json", result["summary"])
    # Exercise rejection/rollback, not just the zero-error constant-flow path.
    nonlinear_flow = AffineFlow([1e-3, 0, 0], np.diag([1000.0, 0, 0]), initial[0])
    check = F6Simulator(
        initial, p, response, basis, flow=nonlinear_flow, acoustic_enabled=False, vessel_enabled=False
    )
    rr = check.run(
        50e-6,
        1e-5,
        [on],
        1e-5,
        adaptive=True,
        adaptive_options=dict(absolute_error_m=5e-12, relative_error=1e-10, minimum_step_s=1.25e-6),
    )
    assert rr["summary"]["status"] == "COMPLETE" and rr["summary"]["rejected_attempts"] > 0
    expected = initial + np.array([1e-6 * np.expm1(0.05), 0, 0])
    np.testing.assert_allclose(rr["positions_m"][-1, 0], expected, rtol=0, atol=1e-10)
    np.testing.assert_allclose(rr["physical_dt_s"].sum(), 50e-6, rtol=1e-12)
    np.savez_compressed(
        out / "ADAPTIVE_REJECTION.npz",
        time_s=rr["time_s"],
        positions_m=rr["positions_m"],
        physical_dt_s=rr["physical_dt_s"],
        attempts=rr["adaptive_attempts"],
    )
    write(out / "ADAPTIVE_REJECTION.json", rr["summary"])
    with h5py.File(F5 / "science_v1/main/trajectory.h5") as h:
        x = h["positions_m"][0, 0, :4]
    initial = np.stack([x, x + np.array([1e-6, 0, 0])])
    event = dict(js(F5 / "ACOUSTIC_RELEASE_PROTOCOL.json")["events"][0], flow_enabled=False)
    event.pop("mode")
    rows = {}
    results = {}
    for mode in ["serial", "batch"]:
        t0 = time.perf_counter()
        valid = 0
        rr = []
        inputs = [initial[0], initial[1]] if mode == "serial" else [initial]
        for i, state in enumerate(inputs):
            sim = F6Simulator(state, p, response, basis, vessel_enabled=False)
            r = sim.run(0.005, 1e-5, [event], 0.0005)
            valid += (
                sum(quality_pass(q, p) for q in r["summary"]["quality_per_sample"])
                if r["summary"]["status"] == "COMPLETE"
                else 0
            )
            save_run(
                out / (mode + str(i) + ".h5"), r, p, identity(), dict(scope="FOUR_BUBBLE_5ms_THROUGHPUT_ONLY")
            )
            rr.append(r)
        elapsed = time.perf_counter() - t0
        rows[mode] = dict(
            valid_samples=valid,
            physical_duration_s=0.005,
            N=4,
            wall_including_setup_metrics_HDF5_s=elapsed,
            valid_samples_per_gpu_hour=valid * 3600 / elapsed,
            scope="warm/cold module load as observed, model setup, graphs, integration, metrics,HDF5 included; shared basis/geometry preparation excluded separately; no claim for 32-bubble 1s navigation",
        )
        results[mode] = rr
    for i in range(2):
        np.testing.assert_allclose(
            results["serial"][i]["positions_m"][:, 0],
            results["batch"][0]["positions_m"][:, i],
            rtol=0,
            atol=1e-15,
        )
    write(
        out / "EFFICIENCY.json",
        dict(
            status="PASS",
            adaptive_event="PASS",
            throughput=rows,
            process_wall_s=time.perf_counter() - start,
            identity=identity(),
        ),
    )
    print("EFFICIENCY CHECKS PASS", flush=True)


if __name__ == "__main__":
    main()
