from __future__ import annotations

import numpy as np

from utils.cfd_flow.repaired_topology_forensics import port_component_membership


def test_port_membership_counts_each_component_without_majority_assignment() -> None:
    labels = np.asarray((0, 0, 1, 1, 1), dtype=np.int64)
    boundaries = {
        "inlet": np.asarray((0, 2)),
        "outlet_01": np.asarray((1,)),
        "outlet_02": np.asarray((3, 4)),
        "outlet_03": np.asarray((0, 1)),
    }
    result = port_component_membership(labels, boundaries, component_count=2)
    assert result["inlet"]["component_boundary_cell_counts"] == {"0": 1, "1": 1}
    assert result["outlet_02"]["component_ids"] == [1]
    assert result["all_ports_share_a_component"] is False
