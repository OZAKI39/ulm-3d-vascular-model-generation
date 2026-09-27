from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

@pytest.mark.parametrize('n',[500,750,2000,2500,5000])
def test_final_formal_cohort_is_prefix_of_p9a4_births(population,core,n):
 actual=cohort(population,n,core['master_seed'])
 assert actual['events']==population['events'][:n]
 assert actual['events'][:500]==core['events']
