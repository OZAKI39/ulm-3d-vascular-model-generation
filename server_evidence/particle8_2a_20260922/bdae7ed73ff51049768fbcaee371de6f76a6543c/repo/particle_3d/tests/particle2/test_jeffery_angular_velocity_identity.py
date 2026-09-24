import numpy as np
from particle_3d.rbc_orientation import angular_velocity,advance_orientation,quaternion_from_short_axis,short_axis


def test_cross_omega_p_matches_independent_jeffery_rhs(geometries):
    rng=np.random.default_rng(222)
    for geometry in geometries:
        for _ in range(30):
            g=rng.normal(size=(3,3))*20
            p=rng.normal(size=3);p/=np.linalg.norm(p)
            e=(g+g.T)/2;w=(g-g.T)/2
            expected=w@p+geometry.jeffery_lambda*(e@p-np.dot(p,e@p)*p)
            omega=angular_velocity(p,geometry.jeffery_lambda,g)
            np.testing.assert_allclose(np.cross(omega,p),expected,atol=1024*np.finfo(float).eps*np.linalg.norm(g),rtol=0)


def test_quaternion_axis_derivative_matches_jeffery_not_only_formula(geometries):
    g=np.array([[1.,12.,-3.],[2.,-4.,5.],[-1.,6.,3.]])
    e=(g+g.T)/2;w=(g-g.T)/2
    q=quaternion_from_short_axis([1,2,3]);p=short_axis(q)
    h=1e-6/np.linalg.norm(g)
    for geometry in geometries:
        omega=angular_velocity(p,geometry.jeffery_lambda,g)
        forward=short_axis(advance_orientation(q,omega,h))
        backward=short_axis(advance_orientation(q,-omega,h))
        measured=(forward-backward)/(2*h)
        expected=w@p+geometry.jeffery_lambda*(e@p-np.dot(p,e@p)*p)
        bound=64*np.finfo(float).eps/h+8*h*h*np.linalg.norm(g)**3
        np.testing.assert_allclose(measured,expected,rtol=0,atol=bound)
