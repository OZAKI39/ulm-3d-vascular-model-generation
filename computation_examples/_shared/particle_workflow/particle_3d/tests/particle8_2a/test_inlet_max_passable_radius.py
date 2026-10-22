import numpy as np
from particle_3d.particle82a_geometry import maximum_handoff_radius, lower_gap, perimeter_edges


def test_inverse_handoff_constraint_and_maximality():
    for clearance in [0., 1e-9, 1e-6, 2.1e-6, 5e-6]:
        a = float(maximum_handoff_radius(clearance))
        if a:
            assert abs(a+lower_gap(a)-clearance) < 1e-20
            assert a+1e-12+lower_gap(a+1e-12) > clearance


def test_internal_triangle_edge_is_not_aperture_edge():
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float)
    edges, ids = perimeter_edges(points, np.array([[0, 1, 2], [0, 2, 3]]))
    assert len(edges) == 4
    assert [0, 2] not in ids.tolist()
