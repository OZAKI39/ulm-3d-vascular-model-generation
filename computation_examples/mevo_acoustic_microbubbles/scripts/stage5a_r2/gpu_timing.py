"""Explicit transfer/synchronization diagnostics, separate from production timing."""

from time import perf_counter
import numpy as np
from scripts.stage5a_r2.run import load
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes


def main():
    start = perf_counter()
    import cupy as cp

    cp.cuda.Stream.null.synchronize()
    context = perf_counter() - start
    rows = []
    for case in ["B2", "B3"]:
        tr, _, cfg, state, _ = load(case, "S1", "numpy", True)
        mask = state["active"]
        x = state["positions_m"][mask]
        r = state["radii_m"][mask]
        pressure = tr.field.interpolate_complex_pressure(x)
        volumes = tr.response.harmonic_response(r, tr.field.frequency_Hz, pressure)["harmonics_m3"]
        frequency = tr.field.frequency_Hz
        density = cfg["medium"]["density_kg_m3"]
        engine = ClassicalFarFieldBjerknes("cupy", enforce_far_field=False)
        reference = ClassicalFarFieldBjerknes("numpy", enforce_far_field=False).forces(
            x, r, volumes, frequency, density
        )
        for repeat in range(21):
            begin = perf_counter()
            devices = [cp.asarray(v) for v in [x, r, volumes]]
            enqueue = perf_counter() - begin
            begin = perf_counter()
            cp.cuda.Stream.null.synchronize()
            input_sync = perf_counter() - begin
            event0, event1 = cp.cuda.Event(), cp.cuda.Event()
            event0.record()
            begin = perf_counter()
            output = engine.forces(*devices, frequency, density)
            api = perf_counter() - begin
            event1.record()
            begin = perf_counter()
            event1.synchronize()
            kernel_sync = perf_counter() - begin
            event_s = float(cp.cuda.get_elapsed_time(event0, event1) / 1000)
            begin = perf_counter()
            host = output.get()
            download = perf_counter() - begin
            np.testing.assert_allclose(host, reference, rtol=2e-12, atol=1e-26)
            begin = perf_counter()
            actual = engine.forces(x, r, volumes, frequency, density).get()
            end_to_end = perf_counter() - begin
            np.testing.assert_array_equal(actual, host)
            rows.append(
                dict(
                    case=case,
                    repeat=repeat,
                    H2D_API_s=enqueue,
                    H2D_explicit_sync_wait_s=input_sync,
                    force_API_including_internal_scalar_sync_s=api,
                    remaining_kernel_sync_wait_s=kernel_sync,
                    D2H_s=download,
                    CUDA_event_inclusive_s=event_s,
                    production_end_to_end_s=end_to_end,
                    input_bytes=sum(v.nbytes for v in [x, r, volumes]),
                    output_bytes=host.nbytes,
                )
            )
    write(
        BASE / "benchmarks/GPU_TRANSFER_SYNC.json",
        dict(
            context_initialization_s=context,
            rows=rows,
            scope="Separate diagnostic. Repeat zero includes first allocation/JIT. Forced phase synchronizations change scheduling; production end-to-end uses unchanged original API. CUDA event time is inclusive and must not be added to host phase times.",
            device="RTX4090",
            dtype="float64/complex128",
            pair_force="Frozen Stage4L A0, unchanged operation order",
        ),
    )
    print(
        {
            case: float(
                np.median([r["production_end_to_end_s"] for r in rows if r["case"] == case and r["repeat"]])
            )
            for case in ["B2", "B3"]
        }
    )


if __name__ == "__main__":
    main()
