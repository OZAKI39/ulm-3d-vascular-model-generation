import numpy as np
from particle_3d.microbubble import MicrobubbleState,equilibrium_state
from particle_3d.particle1_cases import AffineValidationField


def test_missing_half_factor_cannot_pass():
    # Nonzero, distinct components make Omega=curl(u) fail by O(1), not roundoff.
    gradient=np.array([[0.,-3.,-2.],[3.,0.,-1.25],[2.,1.25,0.]])
    sample=AffineValidationField(np.zeros(3),gradient).sample([1e-5,2e-5,3e-5])
    state=equilibrium_state(MicrobubbleState(0,[1e-5,2e-5,3e-5],1e-6,np.zeros(3),np.zeros(3)),sample)
    np.testing.assert_array_equal(state.angular_velocity_s_inv,[1.25,-2.,3.])
    assert not np.allclose(state.angular_velocity_s_inv,sample.vorticity_s_inv,rtol=256*np.finfo(float).eps,atol=256*np.finfo(float).eps)
