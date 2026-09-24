from dataclasses import FrozenInstanceError
import numpy as np
import pytest
from particle_3d.microbubble import MicrobubbleState


def test_si_state_copies_input_and_is_immutable():
    position=np.array([1.,2.,3.]); velocity=np.array([4.,5.,6.]); omega=np.array([7.,8.,9.])
    state=MicrobubbleState(3,position,1e-6,velocity,omega)
    position[:]=0; velocity[:]=0; omega[:]=0
    np.testing.assert_array_equal(state.position_m,[1.,2.,3.])
    np.testing.assert_array_equal(state.velocity_m_s,[4.,5.,6.])
    np.testing.assert_array_equal(state.angular_velocity_s_inv,[7.,8.,9.])
    for array in [state.position_m,state.velocity_m_s,state.angular_velocity_s_inv]:
        assert array.dtype == np.float64 and array.shape == (3,)
        with pytest.raises(ValueError): array.setflags(write=True)
    assert isinstance(state.radius_m,np.float64)
    with pytest.raises(FrozenInstanceError): state.radius_m=2e-6


@pytest.mark.parametrize('radius',[0.,-1.,np.nan,np.inf,[1.],True])
def test_invalid_radius_rejected(radius):
    with pytest.raises(ValueError): MicrobubbleState(0,[0.,0.,0.],radius,[0.,0.,0.],[0.,0.,0.])


@pytest.mark.parametrize('field',['position_m','velocity_m_s','angular_velocity_s_inv'])
@pytest.mark.parametrize('bad',[[0.,np.nan,0.],[np.inf,0.,0.],[1.,2.]])
def test_invalid_vectors_rejected(field,bad):
    kwargs=dict(particle_id=0,position_m=[0.,0.,0.],radius_m=1e-6,velocity_m_s=[0.,0.,0.],angular_velocity_s_inv=[0.,0.,0.])
    kwargs[field]=bad
    with pytest.raises(ValueError): MicrobubbleState(**kwargs)
