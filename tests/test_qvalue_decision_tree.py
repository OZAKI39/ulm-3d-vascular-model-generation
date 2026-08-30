from __future__ import annotations

import pytest

from utils.cfd_flow.qvalue_contract_forensics import classify_qvalue_contract


@pytest.mark.parametrize(
    ("reader_pass", "pipe_abnormal", "vascular_abnormal", "case"),
    [
        (False, False, False, "CASE_A"),
        (True, True, False, "CASE_B"),
        (True, True, True, "CASE_C"),
        (True, False, False, "CASE_D"),
    ],
)
def test_qvalue_decision_tree_is_exclusive(
    reader_pass: bool,
    pipe_abnormal: bool,
    vascular_abnormal: bool,
    case: str,
) -> None:
    _, _, actual_case = classify_qvalue_contract(
        reader_pass=reader_pass,
        pipe_abnormal=pipe_abnormal,
        vascular_abnormal=vascular_abnormal,
    )

    assert actual_case == case
