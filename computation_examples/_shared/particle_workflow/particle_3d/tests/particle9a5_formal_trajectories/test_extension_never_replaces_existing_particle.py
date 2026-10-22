from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_extension_never_replaces_existing_particle(population,core):
 new=cohort(population,750,core['master_seed']);new=deepcopy(new)
 new['events'][10]['q'][0]+=.001
 with pytest.raises(ValueError,match='replace'):require_prefix(core,new)
 with pytest.raises(ValueError):require_prefix(cohort(population,750,core['master_seed']),core)
