from __future__ import annotations

import numpy as np

from utils.cfd_flow.qvalue_contract_forensics import _summary


def test_fallback_only_distribution_fails_continuous_contract() -> None:
    values = np.asarray([0.5] * 98 + [1.0] * 2, dtype=np.float64)
    summary = _summary(values)

    assert summary["valid_q_fraction"] == 1.0
    assert summary["near_0.5_fraction"] == 0.98
    assert summary["near_1.0_fraction"] == 0.02
    assert summary["non_fallback_continuous_fraction"] == 0.0
    assert summary["unique_q_count_at_1e-10"] == 2
