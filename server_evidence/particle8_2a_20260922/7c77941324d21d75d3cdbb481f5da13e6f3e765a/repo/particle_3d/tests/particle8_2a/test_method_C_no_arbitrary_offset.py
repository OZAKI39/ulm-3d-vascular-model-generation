import numpy as np
from particle_3d.particle82a_geometry import first_admissible


def test_search_follows_geometry_not_radius_multiple():
    path = np.array([[0., 0, 0, 0], [1., 0, 0, 10.]])
    for geometric_transition in [1.37, 4.29, 8.13]:
        r = first_admissible(path, lambda p: p[2]-geometric_transition, horizon=10, tolerance=1e-8)
        assert abs(r['s_birth_m']-geometric_transition) <= 1e-8


def test_never_admissible_and_unresolved_are_not_fabricated_births():
    path = np.array([[0., 0, 0, 0], [1., 0, 0, 10.]])
    assert first_admissible(path, lambda p: -20., horizon=10)['status'] == 'NO_ADMISSIBLE_INWARD_LOCATION'
    assert first_admissible(path, lambda p: -np.inf, horizon=10)['status'] == 'NUMERICAL_GEOMETRY_UNRESOLVED'
