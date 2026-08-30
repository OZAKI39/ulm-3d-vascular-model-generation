from __future__ import annotations

import numpy as np

from utils.cfd_flow.qvalue_contract_forensics import (
    MUSUBI_ROOT_WSL,
    _unc,
)
from utils.cfd_flow.port_flux_audit import parse_treelm_side_contract
from utils.cfd_flow.restart_decode import D3Q19_DIRECTIONS


def test_treelm_to_d3q19_direction_mapping() -> None:
    source = (
        _unc(MUSUBI_ROOT_WSL)
        / "tem/source/tem_param_module.f90"
    ).read_text(encoding="utf-8")
    names, offsets = parse_treelm_side_contract(source)

    assert names[:6] == ("W", "S", "B", "E", "N", "T")
    assert np.array_equal(offsets[:18], D3Q19_DIRECTIONS[:18])
