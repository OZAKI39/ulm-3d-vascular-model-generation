from __future__ import annotations

import numpy as np

from utils.cfd_flow.dimensionless_geometry_kernel import compare_qvalue_support


def test_numeric_error_excludes_intersection_classification_mismatches() -> None:
    result = compare_qvalue_support(
        np.asarray((0.2, -1.0, 0.7, -1.0)),
        np.asarray((0.21, 0.4, -1.0, -1.0)),
    )

    assert result["TRUE_INTERSECTION_BOTH"] == 1
    assert result["SEEDER_ONLY_INTERSECTION"] == 1
    assert result["ORACLE_ONLY_INTERSECTION"] == 1
    assert result["NO_INTERSECTION_BOTH"] == 1
    assert np.isclose(result["numeric_q_error"]["rms"], 0.01)
