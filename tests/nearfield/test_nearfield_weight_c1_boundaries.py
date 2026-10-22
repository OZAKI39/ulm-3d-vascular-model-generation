import numpy as np

def test_c1_numerical_endpoints(policy):
 for c in [.01,.05]:
  errors=[]
  for eps in [1e-5,1e-6,1e-7,1e-8]:
   derivatives=[(policy.activation_weight(c+eps)-policy.activation_weight(c))/eps,(policy.activation_weight(c)-policy.activation_weight(c-eps))/eps]
   errors.append(max(abs(x) for x in derivatives))
  assert errors[-1]<2e-5 and np.all(np.diff(errors)<0)
  assert policy.activation_derivative(c)==0
