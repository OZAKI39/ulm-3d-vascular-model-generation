from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_restart_uses_only_hash_matching_output(valid_completion):
 folder,identity,event=valid_completion
 assert completion_matches(folder,identity,event)
 assert not completion_matches(folder,{'source':'wrong'},event)
 changed=deepcopy(event);changed['diameter_um']+=.001
 assert not completion_matches(folder,identity,changed)
 (folder/'trajectory.npz').write_bytes(b'corrupt')
 assert not completion_matches(folder,identity,event)
