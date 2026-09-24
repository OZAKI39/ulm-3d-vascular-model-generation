from __future__ import annotations

import pytest

from utils.cfd_flow.tau1_grid_convergence import assert_launch_allowed


def test_accepted_base_rejects_new_seeder_and_long_musubi() -> None:
    with pytest.raises(PermissionError):
        assert_launch_allowed("base", "seeder")
    with pytest.raises(PermissionError):
        assert_launch_allowed("base", "long_musubi")
    assert_launch_allowed("base", "referee_one_step")
    assert_launch_allowed("coarse", "long_musubi")
