import numpy as np
def test_no_nan_inf_new_smoke(newpaths):
    for e,a in newpaths:assert np.isfinite(a).all()
