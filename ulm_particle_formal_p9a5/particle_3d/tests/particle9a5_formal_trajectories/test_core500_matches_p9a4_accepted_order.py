from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_core500_matches_p9a4_accepted_order(core,population):
 assert core['events']==population['events'][:500]
 assert [e['particle_id'] for e in core['events']]==list(range(1,501))
 assert core['event_sha256']==[content_sha(e) for e in population['events'][:500]]
 assert all(e['source_event_id']>=e['particle_id'] for e in core['events'])
