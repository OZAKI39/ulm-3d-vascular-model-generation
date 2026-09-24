import numpy as np
def test_all_times(scheduler):
 for t in np.linspace(0,.35,31):
  list(scheduler.through(float(t)))
  error=scheduler.mb_clock.cumulative(t)-scheduler.mb_count
  assert 0<=error<1
