import numpy as np


def test_real_nodes_restore_point_data(node_rows):
    assert len(node_rows) >= 1000
    assert {"fixed", "random", "vessel_interior", "INLET", "WALL", "OUTLET_01", "OUTLET_02", "OUTLET_03"} <= {r["group"] for r in node_rows}
    for row in node_rows:
        assert row["inside_lumen"], row
        assert row["velocity_error_m_s"] <= row["velocity_bound_m_s"], row
        assert row["pressure_error_pa"] <= row["pressure_bound_pa"], row
    print("max node velocity / pressure error", max(r["velocity_error_m_s"] for r in node_rows), max(r["pressure_error_pa"] for r in node_rows))


def test_node_tetra_id_is_smallest_incident_cell(real_field, node_rows):
    for row in node_rows[::13]:
        incident = np.flatnonzero(np.any(real_field.tetra == row["node_id"], axis=1))
        assert row["tetra_id"] == int(incident.min())
