import numpy as np
from method_c_test_helpers import source,FixedDistribution

def test_size_bits_frozen_across_forced_rejections():
    dist=FixedDistribution(1e-6);s=source(dist);checker=s.checker
    class RejectFirst:
        def __init__(self):self.calls=0
        def check(self,e,p,active):
            self.calls+=1
            if self.calls<8:return None,'FORCED_TEST_REJECTION',{}
            return checker.check(e,p,active)
    s.checker=RejectFirst();e=s.event(10)
    assert len(e['position_attempts'])>=8 and dist.calls==1
    assert len({a['diameter_float64_hex'] for a in e['position_attempts']})==1
    assert bytes.fromhex(e['position_attempts'][0]['diameter_float64_hex'])==np.float64(e['diameter_m']).tobytes()
