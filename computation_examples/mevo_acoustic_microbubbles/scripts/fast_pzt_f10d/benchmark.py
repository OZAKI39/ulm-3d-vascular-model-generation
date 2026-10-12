"""One controlled legacy/optimized pair per size, with real fresh subprocesses."""

import argparse
import threading
from time import perf_counter
from .common import BASE, objects, identity, case, write, js
from .pipeline import execute, benchmark_records, compare_hdf
from scripts.fast_pzt_f10c.checks import memory_used


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--size", type=int)
    ap.add_argument("--mode", choices=["legacy", "optimized"])
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
        records = benchmark_records(a.size, "O_1000", case("O_1000"))
        timing = execute(out / f"cold_B{a.size}_{a.mode}", records, obj, a.mode)
        stop.set()
        thread.join()
        import warp as wp

        timing.update(
            static_prepare_s=static,
            nvml_sampled_peak_MiB=max(samples),
            warp_allocator_peak_MiB=wp.get_mempool_used_mem_high("cuda:0") / 2**20,
        )
        write(out / f"cold_B{a.size}_{a.mode}" / "COLD.json", timing)
        return
    assert js(BASE / "validate_v1/RESULTS.json")["status"] == "PASS"
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    start = perf_counter()
    obj = objects()
    static = perf_counter() - start
    rows = []
    for size in [16, 64]:
        records = benchmark_records(size, "O_1000", case("O_1000"))
        timing = {}
        for mode in ["legacy", "optimized"]:
            timing[mode] = execute(out / f"resident_B{size}_{mode}", records, obj, mode)
        check = compare_hdf(out / f"resident_B{size}_legacy", out / f"resident_B{size}_optimized")
        assert check["status"] == "PASS"
        row = dict(
            B=size,
            startup="resident",
            static_prepare_s=static,
            modes=timing,
            comparison=check,
            speedup=timing["legacy"]["complete_run_wall_s"] / timing["optimized"]["complete_run_wall_s"],
        )
        rows.append(row)
        write(out / "RESULTS.json", dict(status="RUNNING", rows=rows))
        print(dict(B=size, startup="resident", speedup=row["speedup"]), flush=True)
    # Destroy the resident CUDA context before launching fresh-process tests.
    write(out / "RESIDENT_COMPLETE.json", dict(status="PASS", rows=rows, static_prepare_s=static))


if __name__ == "__main__":
    main()
