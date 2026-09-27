import numpy as np

def test_conditional_inverse_and_partial_bin(distribution):
    class R:
        def __init__(self,u):self.u=u
        def random(self):return self.u
    u=np.r_[np.linspace(0,.99999,5000),np.nextafter(1.,0.)]
    d=np.array([distribution.sample(R(v))[0] for v in u])
    assert np.count_nonzero(d==4e-6)<=1  # at most one endpoint rounding, never tail mass
    assert np.allclose(distribution.cdf(d),u,rtol=0,atol=2e-15)
    bins=distribution.contract['bins'];last=bins[-1]
    assert last['low_um']==3.95 and last['high_um']==4.
    assert abs(last['retained_probability']*distribution.mass/last['original_probability']-.5)<1e-13
    assert distribution.contract['clipping'] is False
