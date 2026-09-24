import numpy as np

def test_every_position_is_time_integrated(wall_run,pair_run,real_saved):
 for d in [wall_run,pair_run,real_saved]:
  prior={p['particle_id']:np.asarray(p['center_m']) for p in d['positions'] if p['time_s']==0};t=0.
  for s in d['states']:
   assert not s['projection']['position_projection'];dt=s['time_s']-t
   for p in s['particles']:
    i=p['particle_id'];expected=prior[i]+dt*np.asarray(p['velocity_m_s'])
    assert np.array_equal(expected,p['center_m']);prior[i]=np.asarray(p['center_m'])
   t=s['time_s']
