import numpy as np,pytest
def test_actual_volumes(scheduler):
 rows=list(scheduler.through(.01)); t=0.
 for e in rows:
  if e['species']=='RBC':
   assert (e['scheduled_time_s']-t)*.45e-12==pytest.approx(e['volume_m3'],rel=2e-13)
   t=e['scheduled_time_s']
 assert np.std([e['volume_m3'] for e in rows])>1e-18
