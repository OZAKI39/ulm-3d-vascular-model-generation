import numpy as np
from particle_3d.flowfield_conservation_diagnosis import clip_polygon,polygon_flux

def test_wall_flux_zero_for_zero_wall_nodes(cube):
    tri=cube['points'][cube['boundaries']['WALL'][0][0]]
    x,u=clip_polygon(tri,np.zeros((3,3)),cube['center'],cube['normal'])
    assert polygon_flux(x,u)[0]==0.
    # A nonzero normal field must be detected: zero result cannot be hard-coded.
    n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n)
    q,area=polygon_flux(tri,np.tile(n,(3,1)))
    assert np.isclose(q,area,rtol=1e-13,atol=0)
