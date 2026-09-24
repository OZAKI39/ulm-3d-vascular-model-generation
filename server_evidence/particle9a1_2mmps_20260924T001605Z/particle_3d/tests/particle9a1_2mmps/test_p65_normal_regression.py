import numpy as np
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle9a_motion import augment_planar_system
from particle_3d.resistance_solver import solve_resistance
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle3_cases import plane_triangle
from particle_3d.wall_gap import wall_gap

def test_unchanged_normal_resistance():
 w=WallGeometry([plane_triangle()])
 for h in [2e-9,3e-9,1e-8,1e-7,1e-6]:
  s={1:Sphere([0,0,1e-6+h],1e-6)};u={1:np.array([0,0,-1e-4,0,0,.2])}
  base=assemble_v1(s,u,.00345312,w);new=augment_planar_system(base,s,.00345312,w,lambda x:np.zeros((3,3)),wall_gap)
  np.testing.assert_allclose(solve_resistance(base).velocity,solve_resistance(new).velocity,rtol=4e-15,atol=1e-20)
