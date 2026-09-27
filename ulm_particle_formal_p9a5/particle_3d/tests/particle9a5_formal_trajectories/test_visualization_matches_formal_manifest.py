from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_visualization_matches_formal_manifest(core,root):
 rows=[{'particle_id':e['particle_id'],'status':'COMPLETED','outlet':'O3'} for e in core['events']]
 # Display exporter uses the same strict merge: all members must be present.
 assert len(merge_rows(rows,[e['particle_id'] for e in core['events']]))==500
 with pytest.raises(ValueError):merge_rows(rows[:48],[e['particle_id'] for e in core['events']])
 manifest=root/REL/'data/visualization_manifest.json'
 if manifest.exists():
  display=json.loads(manifest.read_text());formal=root/REL/'data/FINAL_FORMAL_COHORT.json'
  frozen=json.loads(formal.read_text())
  assert display['formal_cohort_sha256']==digest(formal)
  assert display['all_drawn_particle_ids']==[e['particle_id'] for e in frozen['events']]
  assert display['formal_N']==frozen['count'] and display['figure01_all_members']
  for name,sha in display['output_files'].items():assert digest(root/REL/'figures'/name)==sha
