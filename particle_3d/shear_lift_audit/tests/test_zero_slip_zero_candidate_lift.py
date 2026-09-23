import numpy as np
from audit_math import candidate

def test_zero_slip_or_zero_shear_is_exact_zero():
    assert candidate(1e-6,0,100,.003,1000)==0
    assert candidate(1e-6,1,0,.003,1000)==0
