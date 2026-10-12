"""One real fresh-process cold run per admitted B; not a summed cost estimate."""

import argparse
import os
import subprocess
import sys
import threading
from time import perf_counter
import numpy as np
import warp as wp
from .common import BASE, objects, case, identity, write, js
from .runner import run
from .checks import memory_used
from .batch import compare


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--size", type=int)
    ap.add_argument("--reference-label", help="Serial reference campaign; original campaigns by default")
    a = ap.parse_args()
    out = BASE / a.label
    if a.size:
        stop = threading.Event()
        samples = []

        def monitor():
            while not stop.is_set():
                samples.append(memory_used())
                stop.wait(0.01)

        thread = threading.Thread(target=monitor, daemon=True)
        thread.start()
        start = perf_counter()
        obj = objects()
        static = perf_counter() - start
        cfg = case("O_1000")
        x = np.array(cfg["environment"]["positions_m"])
        variants = [x + np.array([0, dy, 0]) for dy in [0, 2e-7, -2e-7, 4e-7]]
        t = perf_counter()
        timing = run(
            out / f"B{a.size}", np.stack([variants[i % 4] for i in range(a.size)]), cfg["events"], obj
        )
        outer = perf_counter() - t
        stop.set()
        thread.join()
        write(
            out / f"B{a.size}" / "COLD.json",
            dict(
                B=a.size,
                static_prepare_s=static,
                outer_simulation_output_s=outer,
                body_with_static_s=perf_counter() - start,
                peak_observed_gpu_MiB=max(samples),
                memory_samples=len(samples),
                memory_poll_s=0.01,
                warp_allocator_high_water_MiB=wp.get_mempool_used_mem_high("cuda:0") / 2**20,
                scope="Fresh Python process/CUDA context/static input; original on-disk compiler cache retained. NVML peak is sampled; Warp allocator high-water excludes CUDA context/driver allocations.",
                timing=timing,
            ),
        )
        return
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    rows = []
    for size in [1, 4, 16, 64]:
        start = perf_counter()
        with open(out / f"B{size}.log", "w") as log:
            subprocess.run(
                [sys.executable, "-m", "scripts.fast_pzt_f10c.cold", "--label", a.label, "--size", str(size)],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                env=os.environ.copy(),
            )
        row = js(out / f"B{size}" / "COLD.json")
        row["full_subprocess_wall_s"] = perf_counter() - start
        prior = BASE / (a.reference_label or ("batch_v1" if size < 64 else "batch64_v1"))
        checks = [compare(prior / f"B{size}_serial_{i:03d}", out / f"B{size}", i) for i in range(size)]
        assert all(r["bitwise"] for r in checks)
        row["bitwise_all_samples"] = True
        row["valid_samples_per_gpu_hour"] = 3600 * size / row["full_subprocess_wall_s"]
        rows.append(row)
        write(out / "RESULTS.json", dict(status="PASS", rows=rows, training_ready=False))
        print(
            dict(B=size, cold_program_s=row["full_subprocess_wall_s"], peak_MiB=row["peak_observed_gpu_MiB"]),
            flush=True,
        )


if __name__ == "__main__":
    main()
