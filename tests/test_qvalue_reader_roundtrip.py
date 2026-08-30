from __future__ import annotations

from utils.cfd_flow.qvalue_contract_forensics import synthetic_reader_roundtrip


def test_qvalue_reader_known_pattern_roundtrip() -> None:
    result = synthetic_reader_roundtrip()

    assert result["status"] == "PASS"
    assert result["row_mapping_pass"] is True
    assert result["direction_mapping_pass"] is True
    assert result["dtype_pass"] is True
    assert result["endian_pass"] is True
