import numpy as np
from particle_3d.particle82a_geometry import first_admissible


def test_straight_cylinder_first_full_center_is_one_radius_downstream():
    a = 1.2345e-6
    path = np.array([[0., 0, 0, 0], [.01, 0, 0, 10e-6]])
    result = first_admissible(path, lambda p: min(p[2]-a, 5e-6-a-2e-9), horizon=8e-6, tolerance=1e-11)
    assert result['status'] == 'ACCEPTED'
    assert a <= result['s_birth_m'] <= a+1e-11
    assert result['bracket_m'][1]-result['bracket_m'][0] <= 1e-11


def test_left_first_search_does_not_skip_earlier_nonmonotone_window():
    path = np.array([[0., 0, 0, 0], [1., 0, 0, 10.]])
    # Two windows; their boundaries obey the signed distance Lipschitz bound.
    margin = lambda p: max(min(p[2]-2, 3-p[2]), min(p[2]-7, 9-p[2]))
    result = first_admissible(path, margin, horizon=10, tolerance=1e-8)
    assert result['status'] == 'ACCEPTED'
    assert 2 <= result['s_birth_m'] <= 2+1e-8
