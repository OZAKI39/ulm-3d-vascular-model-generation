import numpy as np
from fem3d.cap_remesh import plane_basis,project,lift


def test_projection_lift_preserves_tangential_coordinates_and_known_normal_offset():
    basis=plane_basis([.3,.7,.2]); origin=np.array([1e-4,2e-4,3e-4])
    xy=np.array([[0.,0.],[2e-6,-1e-6],[-3e-6,1e-6]])
    points=lift(xy,origin,basis)
    projected,normal_offset=project(points,origin,basis)
    np.testing.assert_allclose(projected,xy,rtol=0,atol=6e-20)
    assert np.max(abs(normal_offset))<6e-20
    noisy=points+np.array([0,1e-11,-2e-11])[:,None]*basis[2]
    _,offset=project(noisy,origin,basis)
    np.testing.assert_allclose(offset,[0,1e-11,-2e-11],atol=6e-20)
    # Lifting an original noisy rim would move it: production reuses original bytes.
    assert not np.array_equal(lift(project(noisy,origin,basis)[0],origin,basis),noisy)
