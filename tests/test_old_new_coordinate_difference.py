from __future__ import annotations

import numpy as np

from utils.cfd_flow.repaired_topology_forensics import coordinate_set_difference


def test_coordinate_difference_reports_common_old_only_and_new_only() -> None:
    result = coordinate_set_difference(
        np.asarray((1, 2, 3, 5), dtype=np.int64),
        np.asarray((2, 3, 4, 6), dtype=np.int64),
    )
    assert result["common"].tolist() == [2, 3]
    assert result["old_only"].tolist() == [1, 5]
    assert result["new_only"].tolist() == [4, 6]
