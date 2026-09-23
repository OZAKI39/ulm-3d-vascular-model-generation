import numpy as np
from audit_math import wall_coefficient
from particle_3d.particle_shapes import Sphere
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.hydrodynamic_resistance import PhysicalNearField

def test_all_regularization_regions_against_formal_policy():
    a=1e-6;mu=.00345312;sphere=Sphere(np.zeros(3),a)
    for h in [2e-9,3e-9,1e-8,3e-8,4.9999e-8,5e-8,1e-6]:
        spec=PhysicalNearField(1,None,h,np.array([0,1.,0]),1e-19)
        z,n,record=NearFieldRegularizationV1().evaluate(spec,{1:sphere},mu)
        np.testing.assert_allclose(wall_coefficient(a,h,mu),z,rtol=1e-12,atol=0)
        v=np.array([.1,0,.3]);assert z*abs(v@n)==0
