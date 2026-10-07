"""Measured layout evidence; task throughput and trajectory speed stay separate."""

from collections import defaultdict
import numpy as np
from scripts.stage5a_r3.common import BASE, read, write, sha


def main():
    source = BASE / "benchmarks/LAYOUT_ISOLATED.json"
    data = read(source)
    groups = defaultdict(list)
    for row in data["rows"]:
        groups[(row["case"], row["window"], row["layout"], row["optimized"])].append(row)
    layouts = []
    for (case, window, layout, optimized), rows in groups.items():
        warm = [r for r in rows if r["repeat"] > 0]
        layouts.append(
            dict(
                case=case,
                window=window,
                layout=layout,
                optimized=optimized,
                repeats=len(rows),
                status="PASS" if all(r["status"] == "PASS" for r in rows) else "FAIL",
                cold_wall_s=rows[0]["wall_s"],
                warm_median_wall_s=float(np.median([r["wall_s"] for r in warm])),
                warm_wall_samples_s=[r["wall_s"] for r in warm],
                process_cpu_over_wall=sum(r["process_cpu_s"] for r in warm) / sum(r["wall_s"] for r in warm),
                persistent_executor=any(r["persistent_executor"] for r in rows),
                actual_changing_states=True,
            )
        )
    gpu = []
    for case in ["B2", "B3"]:
        rows = [r for r in data["gpu"] if r["case"] == case and r["repeat"] > 0]
        metrics = {
            k: float(np.median([r[k] for r in rows]))
            for k in [
                "CPU_s",
                "GPU_end_to_end_s",
                "H2D_including_sync_s",
                "GPU_API_s",
                "remaining_sync_s",
                "CUDA_event_inclusive_s",
                "D2H_s",
            ]
        }
        gpu.append(
            dict(
                case=case,
                repeats=len(rows),
                median=metrics,
                force_only_speedup=metrics["CPU_s"] / metrics["GPU_end_to_end_s"],
                max_force_error_N=max(r["max_force_error_N"] for r in rows),
                assertion="float64 allclose rtol=2e-12, atol=1e-26",
                timing_scope="CUDA event, GPU API and synchronization are overlapping views. GPU end-to-end includes transfers and synchronization.",
            )
        )
    throughput = read(BASE / "benchmarks/INDEPENDENT_THROUGHPUT.json")
    result = dict(
        status="PASS" if all(r["status"] == "PASS" for r in layouts) else "FAIL",
        layouts=layouts,
        pair_force=gpu,
        independent_throughput={
            k: throughput[k] for k in ["serial_wall_s", "parallel_wall_s", "speedup", "workers", "label"]
        },
        selected=dict(
            QP="global",
            component_workers=1,
            BLAS_threads=1,
            pair_force="cupy",
            reason="Lowest measured changing-state substep wall time in all four tested windows. Component threading overhead exceeded useful concurrency at 186/189 particles.",
        ),
        resource_identity=str(BASE / "environment/IDENTITY.json"),
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        raw_layout_sha256=sha(source),
        raw_throughput_sha256=sha(BASE / "benchmarks/INDEPENDENT_THROUGHPUT.json"),
        limits="No claim of 5x trajectory speedup from force-only or independent-task throughput. Per-process CPU/wall and actual worker execution are used; host heatmaps are not evidence of this solver parallelism.",
        accounting="Cold import/setup separate from warm changing-state samples. Formal end-to-end timing includes initialization and all output.",
    )
    write(BASE / "final/RESOURCE_LAYOUT_BENCHMARK.json", result)
    print(dict(layouts=len(layouts), pair=gpu, throughput=result["independent_throughput"]))


if __name__ == "__main__":
    main()
