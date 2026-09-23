import numpy as np
from audit_math import force_ratio

def test_denominator_never_floored_to_fabricate_ratio():
    r,s=force_ratio([0,1e-32,1e-12,1e-14],[0,1e-30,0,1e-12],1e-25)
    assert np.isnan(r[:3]).all()
    assert np.array_equal(s,[1,1,2,0])
    assert r[3]==.01

def test_lubrication_operator_amplified_roundoff():
    # A resolved candidate must not be called zero merely because lubrication
    # amplifies normal-velocity uncertainty; its ratio remains undefined.
    r,s=force_ratio([1e-16],[1e-20],[1e-18],[1e-25])
    assert np.isnan(r[0]) and s[0]==2
