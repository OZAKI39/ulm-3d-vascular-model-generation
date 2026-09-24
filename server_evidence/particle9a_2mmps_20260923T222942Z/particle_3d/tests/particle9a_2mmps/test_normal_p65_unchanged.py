import numpy as np
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle9a_motion import augment_planar_system
from particle_3d.wall_gap import wall_gap
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle3_cases import plane_triangle
from particle_3d.particle_shapes import Sphere
from particle_3d.resistance_solver import solve_resistance

def test_normal_p65_unchanged():
    wall=WallGeometry([plane_triangle()]);mu=.00345312
    for h in [2e-9,3e-9,1e-8,1e-7,1e-6]:
        shapes={1:Sphere([0,0,1e-6+h],1e-6)};free={1:np.array([0,0,-1e-4,0,0,.2])}
        p65=assemble_v1(shapes,free,mu,wall)
        p9=augment_planar_system(p65,shapes,mu,wall,lambda x:np.zeros((3,3)),wall_gap)
        u65=solve_resistance(p65).velocity;u9=solve_resistance(p9).velocity
        np.testing.assert_allclose(u9,u65,rtol=4e-15,atol=1e-20)
        np.testing.assert_array_equal(p9.planar_matrix.toarray()[:,[2,5]],np.zeros((6,2)))
