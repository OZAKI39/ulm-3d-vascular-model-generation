import numpy as np,math
def test_residual(scheduler):
 events=[]
 for t in np.linspace(.0001,.12,81):
  events.extend(scheduler.through(float(t)))
  volume=math.fsum(e['volume_m3'] for e in events if e['species']=='RBC')
  residual=scheduler.rbc_clock.cumulative(t)-volume
  assert -1e-28<=residual<scheduler.next_rbc['volume_m3']
