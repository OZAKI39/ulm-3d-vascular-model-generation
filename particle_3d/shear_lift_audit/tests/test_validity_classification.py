import numpy as np
from audit_math import validity,dimensionless

def test_wall_length_rejects_even_large_gap_over_radius():
    z=dimensionless(1e-6,1e-8,1.,100e-6,3e-6)
    assert z['h_over_a']>99
    assert not validity(z['Re_p'],z['Re_G'],z['wall_margin'])
    assert validity(1e-8,1e-6,100.)
    assert not validity(1.,1e-6,100.)
    assert not validity(0.,0.,100.)
