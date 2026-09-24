import numpy as np
from particle_3d.rbc_orientation import quaternion_from_short_axis,rotation_matrix,advance_orientation


def test_body_to_world_short_axis_and_positive_z_rotation_direction():
    q=quaternion_from_short_axis([1,0,0])
    np.testing.assert_allclose(rotation_matrix(q)@np.array([0,0,1]),[1,0,0],atol=8*np.finfo(float).eps)
    moved=advance_orientation(q,[0,0,np.pi/2],1.)  # SYNTHETIC_ONLY validation dt
    np.testing.assert_allclose(rotation_matrix(moved)@np.array([0,0,1]),[0,1,0],atol=16*np.finfo(float).eps)


def test_world_increment_left_multiplies_noncommuting_rotation():
    q=quaternion_from_short_axis([1,2,3])
    angle=.3
    rz=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    actual=rotation_matrix(advance_orientation(q,[0,0,3.],.1))
    np.testing.assert_allclose(actual,rz@rotation_matrix(q),atol=32*np.finfo(float).eps,rtol=0)
    assert not np.allclose(actual,rotation_matrix(q)@rz)
