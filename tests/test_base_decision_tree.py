from __future__ import annotations

from utils.cfd_flow.tau1_base import base_decision


def test_referee_failure_forbids_long_base() -> None:
    assert base_decision(referee_pass=False, base_pass=None) == (
        "CFD_FLOW_CONTINUOUS_Q_REFEREE_FAILED"
    )
    assert base_decision(referee_pass=True, base_pass=None) == "RUN_FRESH_BASE_TAU1"
    assert base_decision(referee_pass=True, base_pass=True).endswith("STEADY_PASS")
    assert base_decision(referee_pass=True, base_pass=False).endswith("STEADY_FAILED")
