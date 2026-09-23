import numpy as np
import pytest
from fem3d.cap_remesh import plane_basis


@pytest.mark.parametrize('normal',[[0,0,1],[1,0,0],[.1,.2,.9],[-.9,.02,.13]])
def test_source_normal_basis_is_orthonormal_right_handed_deterministic(normal):
    basis=np.array(plane_basis(normal))
    np.testing.assert_allclose(basis@basis.T,np.eye(3),atol=5e-16)
    np.testing.assert_allclose(np.cross(basis[0],basis[1]),basis[2],atol=5e-16)
    np.testing.assert_array_equal(basis,plane_basis(normal))
