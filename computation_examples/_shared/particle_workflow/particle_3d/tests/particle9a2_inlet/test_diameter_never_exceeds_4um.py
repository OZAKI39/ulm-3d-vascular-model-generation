import numpy as np

def test_source_upper_bound(distribution):
    rng=np.random.default_rng(891)
    a=np.array([distribution.sample(rng)[0] for _ in range(12000)])
    assert np.all((a>0)&(a<=4e-6))
