import numpy as np

def test_conditional_cdf_parity_and_no_clipping(distribution):
 rng=np.random.default_rng(1994);a=np.array([distribution.sample(rng)[0] for _ in range(20000)])
 assert np.all((a>0)&(a<=4e-6)) and not np.any(a==4e-6)
 for d in [.6,1.,1.5,2.,3.,4.]:
  assert distribution.cdf(d*1e-6)==min(1.,distribution.original.cdf(d)/distribution.mass)
  assert abs(np.mean(a<=d*1e-6)-distribution.cdf(d*1e-6))<.018
