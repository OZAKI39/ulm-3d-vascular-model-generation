from __future__ import annotations

from utils.cfd_flow.tau1_grid_convergence import (
    PRIMARY_FLUX_DEFINITION,
    PRIMARY_METRICS,
    build_primary_analyses,
    evaluate_repaired_grid_gate,
)


def test_all_seven_primary_metrics_must_be_asymptotic_and_within_five_percent() -> None:
    records = {}
    for grid, value in (("coarse", 1.25), ("base", 1.1), ("fine", 1.05)):
        records[grid] = {"flux_definition": PRIMARY_FLUX_DEFINITION, **{metric: value for metric in PRIMARY_METRICS}}
    analyses = build_primary_analyses(records)
    result = evaluate_repaired_grid_gate(analyses)
    assert result["status"] == "PASS"
    assert result["primary_trends_asymptotic_monotonic"]
    assert result["base_fine_primary_within_5_percent"]
