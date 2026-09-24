import numpy as np
import pytest
from particle_3d.microbubble import MicrobubbleState,equilibrium_state,stokes_drag_force_n
from particle_3d.particle1_cases import AffineValidationField


@pytest.mark.parametrize('radius',[.25e-6,1e-6,3e-6])
def test_force_free_translation_is_radius_independent(radius):
    velocity=np.array([2e-4,-1e-4,.5e-4]); field=AffineValidationField(velocity,np.zeros((3,3)))
    state=MicrobubbleState(0,[0.,0.,0.],radius,[999.,0.,0.],[0.,0.,0.])
    mobile=equilibrium_state(state,field.sample(state.position_m))
    np.testing.assert_array_equal(mobile.velocity_m_s,velocity)
    np.testing.assert_array_equal(stokes_drag_force_n(mobile.velocity_m_s,velocity,.00345312,radius),np.zeros(3))


def test_drag_opposes_slip_and_is_linear():
    slip=np.array([1e-4,-2e-4,3e-4])
    force=stokes_drag_force_n(slip,[0.,0.,0.],.00345312,1e-6)
    assert force@slip<0
    np.testing.assert_allclose(stokes_drag_force_n(2*slip,[0.,0.,0.],.00345312,1e-6),2*force,rtol=8*np.finfo(float).eps,atol=0)
