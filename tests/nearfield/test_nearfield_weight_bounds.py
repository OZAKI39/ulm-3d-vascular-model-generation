import numpy as np,pytest

def test_bounds(policy):
 assert all(0<=policy.activation_weight(x)<=1 for x in np.linspace(-.01,.1,1001))
@pytest.mark.parametrize('x',[float('nan'),float('inf')])
def test_nonfinite(policy,x):
 with pytest.raises(ValueError):policy.activation_weight(x)
