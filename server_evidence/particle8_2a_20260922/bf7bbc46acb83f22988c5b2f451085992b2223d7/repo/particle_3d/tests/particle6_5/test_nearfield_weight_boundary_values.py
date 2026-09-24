import pytest

@pytest.mark.parametrize('chi,w',[(0,1),(.005,1),(.01,1),(.02,.84375),(.03,.5),(.04,.15625),(.05,0),(.06,0)])
def test_values(policy,chi,w):assert policy.activation_weight(chi)==pytest.approx(w,abs=3e-16)
