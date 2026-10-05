"""Disjoint interleaved workloads avoid recalculating deterministic cache keys."""

from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
import resource
import time

import numpy as np
import psutil

from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBuilder, cached_task
from common import new_run, finish, atomic_json


def main():
    out = new_run(
        "response_table",
        dict(task="worker_benchmark", counts=[4, 6, 7], cases_each=72, seed=20261002),
        ["APECSS v1.8, same focused constitutive manifold"],
    )
    quota = min(psutil.cpu_count(logical=False), len(os.sched_getaffinity(0)))
    q = Path("/sys/fs/cgroup/cpu.max")
    if q.exists():
        a, b = q.read_text().split()
        if a != "max":
            quota = min(quota, int(a) / int(b))
    builder = AdaptiveResponseTableBuilder(
        out,
        "third_party/build/apecss/ultrasound_apecss",
        "configs/stage2/sonovue_parameters.json",
        cache="results/stage2_5/response_table/node_cache",
    )
    rng = np.random.default_rng(20261002)
    points = rng.uniform(builder.bounds[:, 0], builder.bounds[:, 1], (216, 3))
    rows = []
    for offset, workers in enumerate([4, 6, 7]):
        if workers > quota:
            continue
        tasks = [(builder.cache, builder.exe, builder.case(x)) for x in points[offset::3]]
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.perf_counter()
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(cached_task, tasks))
        elapsed = time.perf_counter() - start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
        rows.append(
            dict(
                workers=workers,
                cases=len(tasks),
                elapsed_s=elapsed,
                cases_per_second=len(tasks) / elapsed,
                cpu_seconds=cpu,
                worker_efficiency=cpu / (elapsed * workers),
                new_calls=sum(not r["cache_hit"] for r in results),
                complete=sum(r["status"] == "COMPLETE" for r in results),
                parameters=[r["case"] for r in results],
            )
        )
        print(workers, len(tasks) / elapsed, cpu / (elapsed * workers), flush=True)
    best = max(rows, key=lambda r: r["cases_per_second"])["workers"]
    summary = dict(
        status="BENCHMARK_COMPLETE",
        quota_cores=quota,
        rows=rows,
        selected_workers=best,
        note="Single disjoint random workloads; throughput includes pool startup and raw output I/O. Not identical-case timing; deterministic completed keys never recomputed.",
    )
    atomic_json(Path("results/stage2_5/response_table/ADAPTIVE_TABLE_PERFORMANCE.json"), summary)
    finish(out, summary)


if __name__ == "__main__":
    main()
