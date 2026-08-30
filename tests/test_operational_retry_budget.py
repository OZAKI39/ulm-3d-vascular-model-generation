from __future__ import annotations

import pytest

from utils.cfd_flow.dimensionless_geometry_kernel import OperationalRetryBudget


def test_operational_retry_budget_enforces_stage_and_total_limits() -> None:
    budget = OperationalRetryBudget(stage_limit=2, total_limit=3)
    assert budget.consume("tiny") == (1, 1)
    assert budget.consume("tiny") == (2, 2)
    with pytest.raises(
        RuntimeError, match="CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
    ):
        budget.consume("tiny")
