from __future__ import annotations

from utils.cfd_flow.tau1_grid_convergence import GRID_SPECS


def test_windows_are_rescaled_by_physical_timestep() -> None:
    assert GRID_SPECS["coarse"].short_window_iterations == 70859
    assert GRID_SPECS["coarse"].long_window_iterations == 141717
    assert GRID_SPECS["base"].short_window_iterations == 119751
    assert GRID_SPECS["base"].long_window_iterations == 239502
    assert GRID_SPECS["fine"].short_window_iterations == 202379
    assert GRID_SPECS["fine"].long_window_iterations == 404758
