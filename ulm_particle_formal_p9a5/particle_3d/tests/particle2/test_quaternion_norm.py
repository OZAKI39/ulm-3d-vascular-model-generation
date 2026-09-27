import numpy as np
from particle_3d.rbc_orientation import advance_orientation,quaternion_from_short_axis,short_axis
from particle_3d.particle2_cases import NORM_ATOL


def test_repeated_noncommuting_updates_preserve_q_and_p_norm():
    q=quaternion_from_short_axis([1,2,3])
    for i in range(1500):
        q=advance_orientation(q,np.array([2.,-3.,4.]) if i%2 else np.array([-.2,5.,1.]),.002)
        assert abs(np.linalg.norm(q)-1)<=NORM_ATOL
        assert abs(np.linalg.norm(short_axis(q))-1)<=NORM_ATOL


def test_zero_omega_preserves_rotation_and_normalizes_input():
    q=advance_orientation([2.,0.,0.,0.],[0.,0.,0.],.001)
    np.testing.assert_array_equal(q,[1.,0.,0.,0.])
