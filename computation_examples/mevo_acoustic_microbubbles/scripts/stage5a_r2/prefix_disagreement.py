"""Confirm self convergence on exactly the same clock as a partial legacy comparison."""

import json
import numpy as np
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write, sha


def main():
    rows = []
    for case in ["B2", "B3"]:
        names = [f"{case}_0.005_v4a_partial_{suffix}" for suffix in ["legacy", "refined_legacy"]]
        reports = [BASE / "comparisons" / f"{name}.json" for name in names]
        if not all(p.exists() for p in reports):
            continue
        old_checks = [json.loads(p.read_text()) for p in reports]
        sources = [BASE / "comparisons" / f"{name}_candidate.npz" for name in names]
        nominal, refined = [dict(np.load(p)) for p in sources]
        assert np.array_equal(nominal["requested_time_s"], refined["requested_time_s"])
        assert np.array_equal(nominal["particle_ids"], refined["particle_ids"])
        assert np.array_equal(nominal["positions_m"][0], refined["positions_m"][0])
        rg = abs(nominal["rg_m"][-1] / refined["rg_m"][-1] - 1)
        series = np.linalg.norm(nominal["rg_m"] - refined["rg_m"]) / np.linalg.norm(refined["rg_m"])
        quart = abs(nominal["pair_quantiles_m"][-1, 1:4] / refined["pair_quantiles_m"][-1, 1:4] - 1)
        center = float(np.linalg.norm(nominal["centroid_m"][-1] - refined["centroid_m"][-1]))
        travel = float(np.linalg.norm(np.diff(refined["centroid_m"], axis=0), axis=1).sum())
        budget = max(0.1 * float(refined["radius_m"]), 0.02 * travel)
        # Each old comparison's bound also includes its reference's error;
        # summing them is conservative for these two candidate samplings.
        uncertainty = sum(r["endpoint_sampling_uncertainty_m"] for r in old_checks)
        bound_max = sum(r["candidate_sampling_displacement_bound_max_m"] for r in old_checks)
        series_uncertainty = np.sqrt(len(refined["rg_m"])) * bound_max / np.linalg.norm(refined["rg_m"])
        robust = bool(
            rg + uncertainty / refined["rg_m"][-1] <= 0.02
            and series + series_uncertainty <= 0.02
            and np.all(quart + 2 * uncertainty / refined["pair_quantiles_m"][-1, 1:4] <= 0.02)
            and center + uncertainty <= budget
        )
        disagreement = robust and all(r["status"] == "FAIL" for r in old_checks)
        rows.append(
            dict(
                case=case,
                compared_duration_s=old_checks[0]["compared_duration_s"],
                interpretation="BASELINE_DISAGREEMENT"
                if disagreement
                else "NO_CONFIRMED_BASELINE_DISAGREEMENT",
                same_clock_self_pass_with_sampling_bounds=robust,
                self_Rg_endpoint_relative=float(rg),
                self_Rg_series_relative_L2=float(series),
                self_pair_quartiles_relative=quart.tolist(),
                self_centroid_difference_m=center,
                self_centroid_budget_m=budget,
                endpoint_sampling_uncertainty_m=uncertainty,
                Rg_series_sampling_relative_L2_uncertainty=float(series_uncertainty),
                legacy_statuses=[r["status"] for r in old_checks],
                legacy_centroid_differences_m=[r["centroid_difference_m"] for r in old_checks],
                legacy_centroid_budgets_m=[r["centroid_budget_m"] for r in old_checks],
                source_sha256={str(p.relative_to(BASE)): sha(p) for p in reports + sources},
                eligible_for_complete_5ms_gate=False,
            )
        )
    write(
        BASE / "final/PARTIAL_BASELINE_DISAGREEMENT.json",
        dict(
            rows=rows,
            script_sha256=sha(__file__),
            environment_identity_sha256=sha(BASE / "environment/RENDER_ENVIRONMENT.json"),
            execution_location="Local float64 postprocessing of archived actual samples",
            scope="Additional same-clock analysis of actual partial-prefix samples, using original metrics and conservative sampling bounds. No coordinate interpolation, parameter change or full-duration gate upgrade.",
        ),
    )
    print([(r["case"], r["interpretation"], r["self_centroid_difference_m"]) for r in rows])


if __name__ == "__main__":
    main()
