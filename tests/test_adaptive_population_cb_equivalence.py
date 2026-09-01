from __future__ import annotations

from pathlib import Path

import numpy as np

from utils.cfd_flow.adaptive_population_cb_equivalence import (
    equivalence_verdict,
    compare_pdf_payloads,
)


def _write_pdf(path: Path, values: np.ndarray) -> None:
    values.astype("<f8").tofile(path)


def test_pdf_payload_comparator_reports_bitwise_identity(tmp_path: Path) -> None:
    values = np.full((2, 19), 1.0 / 19.0)
    old = tmp_path / "old.lsb"
    new = tmp_path / "new.lsb"
    _write_pdf(old, values)
    _write_pdf(new, values)
    result = compare_pdf_payloads(
        old,
        new,
        cells=2,
        dx_m=2.0e-7,
        dt_s=2.0e-9,
        pressure_reference_pa=3.0e6,
    )
    assert result["bitwise_equivalent"] is True
    assert result["max_abs_pdf_diff"] == 0.0
    assert result["max_abs_rho_diff"] == 0.0
    assert result["max_abs_u_lat_diff"] == 0.0


def test_equivalence_verdict_uses_machine_precision_gates() -> None:
    metrics = {
        "bitwise_equivalent": False,
        "max_abs_pdf_diff": 5.0e-15,
        "max_abs_rho_diff": 5.0e-14,
        "max_abs_u_lat_diff": 5.0e-14,
    }
    controller = {
        "target_relative_difference": 5.0e-13,
        "controlled_relative_difference": 5.0e-13,
    }
    assert equivalence_verdict(metrics, controller) == "PASS_MACHINE_PRECISION_EQUIVALENT"
    metrics["max_abs_pdf_diff"] = 2.0e-14
    assert equivalence_verdict(metrics, controller) == "FAIL"


def test_bitwise_equivalence_is_preferred() -> None:
    metrics = {
        "bitwise_equivalent": True,
        "max_abs_pdf_diff": 0.0,
        "max_abs_rho_diff": 0.0,
        "max_abs_u_lat_diff": 0.0,
    }
    controller = {"target_relative_difference": 0.0, "controlled_relative_difference": 0.0}
    assert equivalence_verdict(metrics, controller) == "PASS_BITWISE_EQUIVALENT"
