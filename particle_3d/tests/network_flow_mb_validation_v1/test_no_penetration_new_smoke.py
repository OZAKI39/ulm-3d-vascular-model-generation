import numpy as np
def test_no_penetration_new_smoke(newpaths,newenv):
    assert len(newpaths)==30
    for e,a in newpaths:assert np.min(a[:,14])>=-newenv.wall.roundoff_m
