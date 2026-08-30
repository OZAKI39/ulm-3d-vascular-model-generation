from __future__ import annotations

from utils.cfd_flow.qvalue_contract_forensics import (
    build_qvalue_file_format_contract,
)


def test_source_proves_qvalue_file_contract() -> None:
    contract = build_qvalue_file_format_contract()

    assert contract["status"] == "PASS"
    assert contract["contract"]["qval_lsb_numpy_dtype"] == "<f8"
    assert contract["contract"]["n_sides"].endswith("qQQQ=26")
    assert contract["d3q19_direction_mapping_pass"] is True
