import math,numpy as np,pytest
from particle_3d.continuous_infusion import poisson_interval,counter_rng,P9A4_ARRIVAL

def test_inverse_transform_every_value():
 for i in range(1,1001):
  u=float(counter_rng(42,i,P9A4_ARRIVAL).random())
  assert poisson_interval(u,.13)==-math.log1p(-u)/.13

@pytest.mark.parametrize('u,rate',[(-.1,1),(1.,1),(float('nan'),1),(.2,-1),(.2,float('inf'))])
def test_invalid_inputs(u,rate):
 with pytest.raises(ValueError):poisson_interval(u,rate)

def test_zero_and_near_zero():
 with pytest.raises(StopIteration):poisson_interval(.5,0.)
 assert math.isfinite(poisson_interval(.5,1e-300))
 with pytest.raises(FloatingPointError):poisson_interval(.5,1e-320)
 assert poisson_interval(0.,1.)==0.

def test_fixed_seed_moments():
 a=np.array([poisson_interval(counter_rng(2026092594,i,940).random(),3.) for i in range(1,50001)])
 # Predeclared 5% relative bounds, many standard errors at fixed N/seed.
 assert abs(a.mean()*3-1)<.05
 assert abs(a.var()*9-1)<.05
