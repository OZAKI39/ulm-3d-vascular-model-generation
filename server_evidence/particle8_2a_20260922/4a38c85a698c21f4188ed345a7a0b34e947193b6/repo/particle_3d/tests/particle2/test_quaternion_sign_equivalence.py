import numpy as np
from particle_3d.rbc_orientation import rotation_matrix,advance_orientation,quaternion_from_short_axis,shape_axis_distance


def test_q_sign_and_axis_sign_are_handled_as_different_equivalences():
    q=quaternion_from_short_axis([1,2,3])
    np.testing.assert_array_equal(rotation_matrix(q),rotation_matrix(-q))
    a=advance_orientation(q,[1,2,3],.002);b=advance_orientation(-q,[1,2,3],.002)
    np.testing.assert_array_equal(rotation_matrix(a),rotation_matrix(b))
    p=np.array([1.,2.,3.]);p/=np.linalg.norm(p)
    assert shape_axis_distance(p,-p)<=4*np.sqrt(np.finfo(float).eps)
    # p and -p are a shape equivalence, not the same directed vector.
    assert np.linalg.norm(p-(-p))>1.9
