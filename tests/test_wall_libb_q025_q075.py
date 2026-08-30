from __future__ import annotations

import pytest

from utils.cfd_flow.musubi_wall_force_diagnostics import (
    bouzidi_coefficients,
    wall_libb_post_pdf,
)


@pytest.mark.parametrize(
    ("q_value", "expected_coefficients", "expected_pdf"),
    [
        (0.25, {"c_in": 0.0, "c_out": 0.5, "c_neighbor": 0.5}, 5.5),
        (
            0.75,
            {"c_in": 1.0 / 3.0, "c_out": 2.0 / 3.0, "c_neighbor": 0.0},
            7.0 / 3.0,
        ),
    ],
)
def test_wall_libb_q025_q075(
    q_value: float, expected_coefficients: dict[str, float], expected_pdf: float
) -> None:
    assert bouzidi_coefficients(q_value) == pytest.approx(expected_coefficients)
    assert wall_libb_post_pdf(
        q_value=q_value, f_in=1.0, f_out=3.0, f_neighbor=8.0
    ) == pytest.approx(expected_pdf)
