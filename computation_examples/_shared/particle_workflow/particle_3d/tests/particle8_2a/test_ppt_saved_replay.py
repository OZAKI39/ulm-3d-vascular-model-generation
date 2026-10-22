import numpy as np
import pytest
from particle_3d.particle82a_ppt_visuals import sample_at_age


def test_ppt_replay_interpolates_only_saved_positions_and_held_velocity():
    a=np.zeros((3,7));a[:,0]=[0,.1,.4];a[:,1]=[0,1,4];a[:,4]=[0,10,10]
    p,v,k=sample_at_age(a,.25)
    np.testing.assert_allclose(p,[2.5,0,0]);np.testing.assert_array_equal(v,[10,0,0]);assert k==1


@pytest.mark.parametrize('age',[-.01,.41])
def test_ppt_replay_does_not_invent_after_stopped_or_before_birth(age):
    a=np.zeros((2,7));a[:,0]=[0,.4]
    with pytest.raises(ValueError,match='extrapolation'):sample_at_age(a,age)


def test_ppt_endpoint_is_exact_saved_sample():
    a=np.zeros((2,7));a[:,0]=[0,.4];a[-1,1:4]=[1,2,3]
    p,_,_=sample_at_age(a,.4);np.testing.assert_array_equal(p,a[-1,1:4])
