import numpy as np


def test_closed_domain_all_five_boundary_roles(classifications):
    rows = [r for r in classifications if r["kind"] == "boundary"]
    assert {r["role"] for r in rows} == {"WALL", "INLET", "OUTLET_01", "OUTLET_02", "OUTLET_03"}
    assert all(r["inside_lumen"] and r["tetra_id"] >= 0 for r in rows)


def test_shared_face_and_edge_choose_minimal_id_repeatably(real_field, faces):
    for face in faces:
        points = [face["center"], real_field.points_m[face["nodes"][:2]].mean(axis=0)]
        for point in points:
            # Exhaustive independent candidate enumeration checks locator completeness.
            weights = real_field.geometry.weights(point, np.arange(len(real_field.tetra)))
            ids = np.flatnonzero(real_field.geometry.contains(weights, np.arange(len(real_field.tetra))))
            assert len(ids) >= 2
            assert [real_field.sample(point).tetra_id for _ in range(4)] == [int(ids.min())] * 4
