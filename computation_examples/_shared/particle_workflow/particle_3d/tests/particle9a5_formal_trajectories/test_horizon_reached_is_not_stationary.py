from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_horizon_reached_is_not_stationary():
 m=dict(completed=False,end_reason='PHYSICAL_RESIDENCE_HORIZON_REACHED',last_elapsed_time_s=3.)
 assert terminal_status(m,True,12.)=='HORIZON_REACHED'
 m['last_elapsed_time_s']=12.
 assert terminal_status(m,True,12.)=='LONG_RESIDENCE_CENSORED'
