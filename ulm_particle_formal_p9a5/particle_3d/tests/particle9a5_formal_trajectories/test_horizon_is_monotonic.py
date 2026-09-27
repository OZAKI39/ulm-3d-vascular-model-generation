from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_horizon_is_monotonic():
 visited=[];s=HorizonSchedule(checkpoint=lambda a,b:visited.append((a,b)))
 steps=list(s.steps());assert steps==list(range(1,round(12./DT)+1))
 assert visited==[(3.,6.),(6.,12.)] and s.used==12.
 for schedule in [(12.,6.),(3.,1.5),(3.,3.)]:
  with pytest.raises(ValueError,match='increase'):HorizonSchedule(schedule)
 # Early natural termination cannot claim later horizons were used.
 early=HorizonSchedule();iterator=early.steps();next(iterator)
 assert early.used==3. and not early.extensions
