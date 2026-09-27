from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_extended_trajectory_preserves_old_prefix(monkeypatch,tmp_path,core,real_env):
 import numpy as np
 import particle_3d.formal_dynamics_p9a5 as d
 # Short physical-age schedule exercises exactly the production continuation
 # path on the real NEW field without a 12-second permanent-test cost.
 schedule=(.002,.004,.008);original=HorizonSchedule
 monkeypatch.setattr(d,'HORIZONS',schedule)
 monkeypatch.setattr(d,'HorizonSchedule',lambda **kw:original(horizons=schedule,**kw))
 monkeypatch.setattr(d,'ENV',real_env);monkeypatch.setattr(d,'IDENTITY',{'role':'SHORT_PREFIX_PERMANENT_TEST'})
 row=d.job((core['events'][0],str(tmp_path/'track')))
 a=np.load(tmp_path/'track/trajectory.npz')['samples']
 for h in [.002,.004]:
  prefix=np.load(tmp_path/'track/horizons'/f'age_{h:g}s.npz')['samples']
  assert np.array_equal(prefix,a[:len(prefix)])
 assert row['horizon_extensions']==[[.002,.004],[.004,.008]]
 assert row['status']=='LONG_RESIDENCE_CENSORED'
