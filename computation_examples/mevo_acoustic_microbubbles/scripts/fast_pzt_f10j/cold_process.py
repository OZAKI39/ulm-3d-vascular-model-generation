"""One independent process per implementation, warm compiler cache, no repetition."""

import argparse
import subprocess
import sys
from time import perf_counter
from .common import ROOT, BASE, js, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--method", choices=["original", "corrected", "optimized", "cpu"])
    args = ap.parse_args()
    out = BASE / args.label
    if args.method is None:
        assert js(BASE / "optimized_v2/RESULTS.json")["numerical_admitted"]
        out.mkdir(parents=True, exist_ok=False)
        write(out / "IDENTITY.json", identity())
        rows = []
        for method in ["original", "corrected", "optimized", "cpu"]:
            busy = subprocess.check_output(
                ["nvidia-smi", "--query-compute-apps=pid", "--format=csv,noheader"], text=True
            ).strip()
            assert not busy, "Wait for other GPU users; never stop their jobs"
            start = perf_counter()
            with (out / (method + ".log")).open("w") as log:
                result = subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "scripts.fast_pzt_f10j.cold_process",
                        "--label",
                        args.label,
                        "--method",
                        method,
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
            rows.append(
                dict(method=method, process_wall_s=perf_counter() - start, exit_code=result.returncode)
            )
            assert result.returncode == 0
        write(
            out / "RESULTS.json",
            dict(
                rows=rows,
                repeats=1,
                cache="Existing pinned Warp compiler cache; new process, context, buffers and graphs",
                scope="Includes identical shared validation-harness imports, input loading, construction, one solve, independent quality and matched output, JSON evidence save, shutdown. CPU harness imports Warp but creates no CUDA context. Diagnostic startup cost, not steady-state speed.",
            ),
        )
        return
    import numpy as np
    from .performance import common_output
    from mevo_mb.fast_sim.f10g_contact import resistance
    from mevo_mb.fast_sim.f10h_contact import geometric_contact
    from mevo_mb.fast_sim.f10i_gpu import FixedStateB2
    from mevo_mb.fast_sim.f10j_gpu import CorrectedB2, OptimizedB2

    a = js(ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json")["rows"][-1]
    x, free = np.array(a["initial_positions_m"]), np.array(a["free_velocity_m_s"])
    ref = np.array(js(BASE / "corrected_v1/RESULTS.json")["rows"][-1]["independent_reference"]["positions_m"])
    start = perf_counter()
    A = resistance(x, 1.78e-6)["A"]
    assembly = perf_counter() - start
    construction = 0.0
    if args.method == "cpu":
        result = geometric_contact(
            x, free, 1.78e-6, 6 * np.pi * 0.001 * 1.78e-6, 1e-5, 0.04641319518378658, A, atol_m=1e-15
        )
    else:
        import warp as wp

        wp.init()
        t = perf_counter()
        g = {"original": FixedStateB2, "corrected": CorrectedB2, "optimized": OptimizedB2}[args.method](32)
        construction = perf_counter() - t
        result = g.solve(x, free)
    output = common_output(result, A, x, free, ref)
    if args.method != "original":
        assert result["converged"] and output["quality_pass"]
    write(
        out / (args.method + ".json"),
        dict(
            method=args.method,
            cpu_reference_assembly_s=assembly,
            construction_s=construction,
            result=result,
            matched_output=output,
        ),
    )


if __name__ == "__main__":
    main()
