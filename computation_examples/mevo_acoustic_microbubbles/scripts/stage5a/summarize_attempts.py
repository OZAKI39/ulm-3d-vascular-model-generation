"""Aggregate accepted/rejected execution ledgers without changing any run."""

from pathlib import Path
import json
import collections
import numpy as np
from scripts.stage4l.dense_probe import write, sha

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"


def dist(values):
    if not values:
        return None
    return dict(
        zip(["min", "p05", "median", "p95", "max"], map(float, np.quantile(values, [0, 0.05, 0.5, 0.95, 1])))
    )


def main():
    (BASE / "final").mkdir(parents=True, exist_ok=True)
    results = []
    for p in sorted((BASE / "runs").glob("*/SUMMARY.json")):
        folder = p.parent
        summary = json.load(open(p))
        accepted = [json.loads(s) for s in (folder / "ACCEPTED.jsonl").open()]
        rejected = [json.loads(s) for s in (folder / "REJECTED.jsonl").open()]
        allrows = accepted + [v["diagnostics"] for v in rejected if "predictor" in v["diagnostics"]]
        modes = {}
        for mode in ["predictor", "corrector"]:
            calls = [r[mode] for r in allrows]
            components = [s for r in calls for s in r["solves"]]
            modes[mode] = dict(
                calls=len(calls),
                time_s=sum(r["qp_seconds"] for r in calls),
                setup_s=sum(r["setup_s"] for r in calls),
                solve_s=sum(r["solve_s"] for r in calls),
                polish_s=sum(r["polish_s"] for r in calls),
                component_qps=len(components),
                primal_fallback_qps=sum(s.get("primal_fallback", False) for s in components),
                active_iterations=dist(
                    [s["active_iterations"] for s in components if "active_iterations" in s]
                ),
                dual_iterations=dist([s["iterations"] for s in components]),
                constraints=dist([r["constraint_count"] for r in calls]),
                component_size_histogram=dict(
                    sorted(
                        collections.Counter(str(n) for r in calls for n in r["component_sizes"]).items(),
                        key=lambda p: int(p[0]),
                    )
                ),
                warm_start_calls=sum(not s["cold_start"] for s in components),
                primal_warm_calls=sum(s["primal_warm_start"] for s in components),
                workspace_reused_calls=sum(s.get("workspace_reused", False) for s in components),
            )
        results.append(
            dict(
                run=str(folder),
                summary_sha256=sha(p),
                case=summary["case"],
                duration_s=summary["duration_s"],
                macro_dt_s=summary["macro_dt_s"],
                status=summary["status"],
                accepted_steps=len(accepted),
                rejected_intervals=len(rejected),
                rejection_reasons=dict(collections.Counter(v["reason"] for v in rejected)),
                accepted_safety=dict(
                    kkt_max=max((r["kkt"] for r in accepted), default=None),
                    linear_residual_max=max((r["linear_residual"] for r in accepted), default=None),
                    projection_calls=sum(r["projection_calls"] for r in accepted),
                    projection_corrected_calls=sum(r["projection_corrected_calls"] for r in accepted),
                    projection_max_m=max((r["projection_max_m"] for r in accepted), default=None),
                    projection_rms_max_m=max((r["projection_rms_m"] for r in accepted), default=None),
                ),
                accepted_error_ratio=dist([r["error_ratio"] for r in accepted]),
                rejected_error_ratio=dist(
                    [r["diagnostics"]["error_ratio"] for r in rejected if "error_ratio" in r["diagnostics"]]
                ),
                solvers=modes,
                failed_solver=summary["solver_failure"],
                coverage="Completed QP calls in accepted and local-error rejected trials; the final failed QP is recorded separately and its incomplete polish time is not included.",
            )
        )
    write(BASE / "final/ATTEMPT_ANALYSIS.json", results)
    log = BASE / "CAMPAIGN.log"
    warnings = [
        dict(line_number=i + 1, text=line.strip())
        for i, line in enumerate(log.read_text().splitlines())
        if "illegal value" in line or "warning" in line.lower()
    ]
    write(
        BASE / "final/NATIVE_WARNINGS.json",
        dict(
            source=str(log),
            source_sha256=sha(log),
            warnings=warnings,
            counts=dict(collections.Counter(row["text"] for row in warnings)),
            interpretation="Preserved native/library messages; exact call site is not established. Independent accepted-state checks do not prove that the native warning is harmless. The scientific numerical gate remains failed.",
        ),
    )


if __name__ == "__main__":
    main()
