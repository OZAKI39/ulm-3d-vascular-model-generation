from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_core500_does_not_use_future_outlet(population,core):
 p=deepcopy(population)
 for e in p['events']:e['future_outlet']='O1' if e['particle_id']>500 else 'O3'
 selected=cohort(p,500,core['master_seed'])
 assert [e['particle_id'] for e in selected['events']]==list(range(1,501))
 assert 'point' not in str(__import__('inspect').signature(cohort))
 assert all('outlet' not in k.lower() and 'basin' not in k.lower() for e in core['events'] for k in e)
