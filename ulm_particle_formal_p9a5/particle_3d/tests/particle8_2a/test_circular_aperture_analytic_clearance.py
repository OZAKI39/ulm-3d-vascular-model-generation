import numpy as np
import pytest
from particle_3d.particle82a_geometry import cylinder_triangles, segment_distance
from particle_3d.wall_geometry import WallGeometry
from particle_3d.injection_admission import FiniteSizeAdmission


@pytest.mark.parametrize('sides', [32, 64, 128])
def test_inscribed_polygon_has_known_analytic_clearance(sides):
    R = 5e-6
    wall, cap = cylinder_triangles(radius=R, sides=sides)
    edges = cap[:, [1, 2]]
    exact = R*np.cos(np.pi/sides)
    assert abs(segment_distance(np.zeros(3), edges)-exact) < 1e-18
    assert abs(WallGeometry(wall).nearest_center_triangle(np.zeros(3))[1]-exact) < 1e-18


def test_near_edge_sphere_rejected():
    wall, _ = cylinder_triangles(radius=5e-6)
    checker = FiniteSizeAdmission(None, None, wall=WallGeometry(wall))
    event = dict(species='MB', radius_m=1e-6, particle_id=1, q=[1, 0, 0, 0])
    assert checker.check(event, np.array([4.5e-6, 0, 0]), {})[1] == 'WALL_REJECTED'


@pytest.mark.parametrize('angle', [0., np.pi/128, np.pi/64])
def test_face_edge_vertex_do_not_create_false_rejections(angle):
    wall, _ = cylinder_triangles(radius=5e-6, sides=128)
    checker = FiniteSizeAdmission(None, None, wall=WallGeometry(wall))
    event = dict(species='MB', radius_m=.5e-6, particle_id=1, q=[1, 0, 0, 0])
    p = np.array([3e-6*np.cos(angle), 3e-6*np.sin(angle), 0.])
    assert checker.check(event, p, {})[1] == 'ACCEPTED'
