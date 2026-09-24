import numpy as np
import pytest
from particle_3d.inlet_flux import InletFluxSampler
T=np.array([[[0,0,0],[1,0,0],[0,1,0]]],float)
@pytest.mark.parametrize('q,expected',[([1,1,1],.5),([-1,1,1],5/24),([-1,-1,1],1/24),([0,0,1],1/6)])
def test_clipped(q,expected):
 assert InletFluxSampler(T,[q]).Q_m3_s==pytest.approx(expected,abs=2e-16)
def test_no_flux():
 with pytest.raises(ValueError): InletFluxSampler(T,[[-1,0,0]])
def test_frozen_consistent(real):
 audit=real[4]
 for role,d in audit.items():
  if role=='mass_balance': continue
  assert d['positive_Q_m3_s']>=d['signed_Q_m3_s']-1e-28
 assert abs(audit['mass_balance']['relative_signed_residual'])<1e-5
