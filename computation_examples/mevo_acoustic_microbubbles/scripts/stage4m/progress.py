"""Read-only progress bounds; never change or terminate a numerical run."""

import json
from scripts.stage4m.common import BASE, write


def snapshot():
    rows = []
    a0 = {
        "B2": {"steps": 17385, "wall_s": 2216.0109830279835},
        "B3": {"steps": 18627, "wall_s": 3080.2959012428764},
    }
    for p in sorted((BASE / "dense_probe").glob("*/*/PROGRESS.json")):
        d = json.loads(p.read_text())
        policy = p.parent.parent.name
        reference = a0[d["case"]]
        row = dict(
            case=d["case"],
            policy=policy,
            progress_percent=100 * d["elapsed_probe_s"] / 0.005,
            accepted_steps=d["steps"],
            wall_s=d["wall_s"],
            latest_dt_s=d["latest_dt_s"],
            monotonic_step_reduction_upper_bound=reference["steps"] / d["steps"],
            monotonic_wall_speedup_upper_bound=reference["wall_s"] / d["wall_s"],
        )
        row["fivefold_performance_already_unreachable"] = (
            policy == "A2"
            and row["monotonic_step_reduction_upper_bound"] < 5
            and row["monotonic_wall_speedup_upper_bound"] < 5
        )
        rows.append(row)
    unreachable = any(r["fivefold_performance_already_unreachable"] for r in rows)
    result = dict(
        rows=rows,
        A2_engineering_gate_already_unreachable=unreachable,
        scope="Counters and elapsed time are monotonic. Bounds only rule out later fivefold performance; they do not substitute for completed 5 ms diagnostics or stability checks.",
    )
    write(BASE / "comparisons/MONOTONIC_RUNTIME_BOUND.json", result)
    print(json.dumps(result), flush=True)
    return result


if __name__ == "__main__":
    snapshot()
