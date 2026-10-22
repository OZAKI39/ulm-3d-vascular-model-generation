import numpy as np

def test_monotonic(policy):
 assert np.all(np.diff([policy.activation_weight(x) for x in np.linspace(0,.06,6001)])<=0)
