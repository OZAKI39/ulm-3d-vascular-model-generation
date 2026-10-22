from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_point_tracer_never_changes_formal_membership(population,core,tmp_path,root):
 frozen=tmp_path/'cohort.json';before=write_new(frozen,core)
 # Even an oracle predicting a rare outlet only beyond CORE500 cannot alter it.
 point={e['particle_id']:('O1' if e['particle_id']>500 else 'O2') for e in population['events']}
 assert all(point[e['particle_id']]=='O2' for e in core['events'])
 assert cohort(population,500,core['master_seed'])==core and digest(frozen)==before
 saved=root/REL/'data/point_metrics.json'
 if saved.exists():
  final=root/REL/'data/FINAL_FORMAL_COHORT.json';events=json.loads(final.read_text())['events'];points=json.loads(saved.read_text())
  assert len(points)==len(events)
  for e,p in zip(events,points):
   assert p['cohort_sha256']==digest(final) and p['cohort_preexists_point_execution']
   assert p['particle_id']==e['particle_id'] and p['initial_center_m']==e['birth_center_m']
