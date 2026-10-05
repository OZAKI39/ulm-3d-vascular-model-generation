"""Render completed P0 profile evidence and retain readable cProfile call tables."""
import io
import json
import pstats

from scripts.stage4i5.common import BASE,new_run,finish


def main():
    out=new_run("profiles_report",dict(method="Exclusive explicit timers plus cProfile; no profiler runtime used as speedup denominator"))
    results=[]
    lines=["# Stage 4I.5 baseline profile","",
        "P0 includes the per-particle exact interior-wall certificate already present in gold. These costs must not be compared to the older unaccelerated Stage 4 profile as if it were current gold.",
        "Initialization, scientific stepping and HDF output are reported separately. No VTK frame output occurs in the original scientific loop. Python overhead includes unwrapped work and instrumentation overhead, not an independently exact physical component.",
        "VTK already builds its original locator once and reuses generic cells. Value and gradient calls each perform a closest-point search.", ""]
    for path in sorted((BASE/"profiles").glob("*/SUMMARY.json")):
        row=json.loads(path.read_text())
        if row.get("status")!="COMPLETE" or row.get("backend")!="P0":
            continue
        results.append(row)
        name=row["segment"]["id"]
        lines += ["## "+name,"",f"Accepted steps: {row['accepted_steps']}; physical interval {row['start_time_s']}–{row['end_time_s']} s. Instrumented wall time: {row['elapsed_s']:.3f} s.",
                  f"Gold endpoint maximum position difference: {row['gold_endpoint_max_position_error_m']:.8g} m. This reconstruction check is separate from candidate-vs-P0 agreement.","",
                  "| Exclusive component | Seconds | Calls |","|---|---:|---:|"]
        for key,value in sorted(row["timers_exclusive_s"].items(),key=lambda kv:-kv[1]):
            lines.append(f"| {key} | {value:.6f} | {row['calls'].get(key,0)} |")
        lines += [f"| Unclassified Python/timer overhead | {row['unclassified_python_and_timer_overhead_s']:.6f} | — |",
                  f"| HDF output (outside loop) | {row['output_HDF_s']:.6f} | 1 |",""]
        stream=io.StringIO()
        pstats.Stats(str(path.parent/"profile.pstats"),stream=stream).strip_dirs().sort_stats("cumulative").print_stats(45)
        (out/(name+"_CPROFILE.txt")).write_text(stream.getvalue())
    (out/"PROFILE_STAGE4I5_BASELINE.md").write_text("\n".join(lines)+"\n")
    finish(out,dict(status="COMPLETE" if len(results)==4 else "PARTIAL",profiles=results))


if __name__=="__main__":
    main()
