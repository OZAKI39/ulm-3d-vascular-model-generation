import numpy as np
import pytest
from particle_3d.microbubble import MicrobubbleState
from particle_3d.integrator import advance_single_microbubble


@pytest.mark.parametrize('dt',[0.,-1.,np.nan,np.inf,-np.inf,True,[.001]])
def test_invalid_dt_rejected_before_sampling(dt):
    class DoNotSample:
        def sample(self,p): raise AssertionError('bad dt must fail before field access')
    state=MicrobubbleState(0,[0.,0.,0.],1e-6,[0.,0.,0.],[0.,0.,0.])
    with pytest.raises(ValueError,match='dt_s'): advance_single_microbubble(state,DoNotSample(),dt)
