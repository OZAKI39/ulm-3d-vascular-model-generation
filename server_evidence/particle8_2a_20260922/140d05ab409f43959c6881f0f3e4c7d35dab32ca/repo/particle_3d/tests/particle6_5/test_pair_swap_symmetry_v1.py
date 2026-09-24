import numpy as np,pytest
from particle_3d.particle_shapes import Sphere
from particle_3d.particle65_motion import assemble_v1
from particle_3d.resistance_solver import solve_resistance

@pytest.mark.parametrize('ratio',[1,2,4])
def test_swap(policy,ratio):
 a=Sphere([0,0,0],1e-6);b=Sphere([1e-6*(1+ratio)+3e-9,0,0],ratio*1e-6)
 u=np.array([1e-5,2e-5,0,1,2,3]);v=-u
 s=assemble_v1({1:a,2:b},{1:u,2:v},.00345312,policy=policy);t=assemble_v1({1:b,2:a},{1:v,2:u},.00345312,policy=policy)
 perm=list(range(6,12))+list(range(6));assert np.array_equal(s.matrix.toarray(),t.matrix.toarray()[np.ix_(perm,perm)])
 assert np.allclose(solve_resistance(s).velocity,solve_resistance(t).velocity[perm],rtol=1e-12,atol=1e-20)
