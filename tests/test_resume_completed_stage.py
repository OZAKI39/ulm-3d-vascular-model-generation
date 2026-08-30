from __future__ import annotations

from utils.cfd_flow.dimensionless_geometry_kernel import completed_stage_reusable


def test_resume_reuses_only_pass_stage_with_identical_inputs() -> None:
    checkpoint = {"status": "PASS", "input_hashes": {"mesh": "abc"}}

    assert completed_stage_reusable(checkpoint, {"mesh": "abc"}) is True
    assert completed_stage_reusable(checkpoint, {"mesh": "changed"}) is False
    assert (
        completed_stage_reusable({**checkpoint, "status": "FAIL"}, {"mesh": "abc"})
        is False
    )
