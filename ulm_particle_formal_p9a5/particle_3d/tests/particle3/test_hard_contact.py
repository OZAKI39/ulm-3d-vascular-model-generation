import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere,Ellipsoid
from particle_3d.rbc_orientation import rotation_matrix,quaternion_from_short_axis
from particle_3d.wall_gap import wall_gap,touching_contacts
from particle_3d.wall_contact import contact_velocity,constrain_contacts


@pytest.mark.parametrize('v,expected',[([2,3,-4],[2,3,0]),([2,3,0],[2,3,0]),([2,3,4],[2,3,4])])
def test_sphere_incoming_parallel_separating_preserves_tangent(plane_wall,v,expected):
    sphere=Sphere([0,0,1e-6],1e-6);gap=wall_gap(sphere,plane_wall);v=np.asarray(v)*1e-4
    omega=np.array([3.,-4.,5.])
    result,w=contact_velocity(v,omega,sphere.center_m,gap)
    np.testing.assert_allclose(result,np.array(expected)*1e-4,atol=1e-18,rtol=0)
    np.testing.assert_array_equal(w,omega)
    general=constrain_contacts(sphere,v,omega,touching_contacts(sphere,plane_wall))[0]
    np.testing.assert_allclose(general,result,atol=1e-18,rtol=0)


def test_rotating_ellipsoid_contact_point_velocity(plane_wall):
    rot=rotation_matrix(quaternion_from_short_axis([1,0,1]));axes=np.array([2,2,.5])*1e-6
    height=np.linalg.norm(axes*(rot.T@np.array([0,0,1])))
    shape=Ellipsoid([0,0,height],axes,rot);gap=wall_gap(shape,plane_wall)
    v=np.array([1e-4,2e-4,-1e-4]);w=np.array([0.,20.,0.])
    corrected,omega=contact_velocity(v,w,shape.center_m,gap)
    normal_speed=gap.normal_inward@(corrected+np.cross(omega,gap.particle_point_m-shape.center_m))
    assert abs(normal_speed)<1e-18
    np.testing.assert_array_equal(omega,w)
    np.testing.assert_allclose(corrected[:2],v[:2],atol=1e-18,rtol=0)


def test_no_contact_has_exactly_no_correction(plane_wall):
    s=Sphere([0,0,2e-6],1e-6);g=wall_gap(s,plane_wall)
    v=np.array([1.,2.,-3.]);omega=np.array([3.,2.,1.])
    got,w=contact_velocity(v,omega,s.center_m,g)
    np.testing.assert_array_equal(got,v);np.testing.assert_array_equal(w,omega)
