from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_new_formal_horizon_not_shorter_than_historical_formal(root):
 import csv
 rows=list(csv.DictReader((root/'particle_3d/outputs/particle8_2a_ppt/data/trajectory_catalog.csv').open()))
 old={float(r['last_age_s']) for r in rows if r['end_reason']=='PHYSICAL_RESIDENCE_HORIZON_REACHED'}
 assert old=={1.5} and HORIZONS[0]>=max(old) and DT==.001
 recorded=json.loads((root/REL/'data/historical_horizon_audit.json').read_text())
 assert recorded['previous_formal_horizon_s']==max(old)
