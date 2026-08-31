from __future__ import annotations

from utils.cfd_flow.tau1_base import (
    OLD_DT_S,
    Tau1BaseRuntimeContract,
    rescale_physical_window,
)


def test_old_to_new_window_preserves_physical_duration() -> None:
    contract = Tau1BaseRuntimeContract()
    for old_iterations in (100, 5_000, 10_000, 20_000):
        new_iterations = rescale_physical_window(old_iterations)
        error_s = abs(new_iterations * contract.dt_s - old_iterations * OLD_DT_S)
        assert error_s <= 0.5 * contract.dt_s
