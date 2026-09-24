import numpy as np
from conftest import balance,grid
from particle_3d.interior_section import cut_tetrahedra,p1_flux

def test_section_two_independent_integrators_agree(cube):
    u=.001*(cube['unit']@np.array([[2,3,4],[-1,2,1],[0,1,-3]]))+[.001,.002,.003]
    r=balance(cube,u)
    q=p1_flux(cut_tetrahedra(grid(cube,u),cube['center'],cube['normal']),cube['normal'])['signed_Q_m3_s']
    assert np.isclose(q,r['Q_section_numpy_m3_s'],rtol=2e-13,atol=0)
