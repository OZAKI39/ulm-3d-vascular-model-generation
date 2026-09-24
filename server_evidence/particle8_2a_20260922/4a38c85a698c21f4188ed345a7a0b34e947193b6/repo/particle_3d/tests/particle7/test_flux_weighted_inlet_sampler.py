import numpy as np
import pytest
from particle_3d.inlet_flux import InletFluxSampler
@pytest.mark.parametrize('q,mean',[([1,1,1],[1/3,1/3,0]),([0,1,0],[.5,.25,0])])
def test_analytic(q,mean):
 s=InletFluxSampler([[[0,0,0],[1,0,0],[0,1,0]]],[q]); p,_=s.sample(np.random.default_rng(8),100000)
 assert np.max(abs(p.mean(0)-mean))<.004
 assert np.max(abs(s.expectation()-mean))<1e-14
 assert (p[:,:2]>=0).all() and (p[:,:2].sum(1)<=1+1e-14).all()
