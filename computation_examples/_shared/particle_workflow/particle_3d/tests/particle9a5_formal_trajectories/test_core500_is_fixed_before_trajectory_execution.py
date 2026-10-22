from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_core500_is_fixed_before_trajectory_execution(root,core,tmp_path):
 p=root/REL/'data/CORE500_COHORT.json';expected=json.loads((root/'particle_3d/contracts/P9A5_FORMAL_PRODUCTION_V1.json').read_text())['core_cohort_sha256']
 assert digest(p)==expected and len(core['events'])==500
 local=tmp_path/'CORE500_COHORT.json';write_new(local,core)
 changed=deepcopy(core);changed['events'][0]['diameter_um']+=.1
 with pytest.raises(ValueError,match='Refuse replacement'):write_new(local,changed)
 assert digest(local)==expected
