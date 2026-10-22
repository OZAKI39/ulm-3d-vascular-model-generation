from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_partial_output_not_counted_as_completed(tmp_path,core):
 for name in ['trajectory.npz','trajectory.json']:(tmp_path/name).write_bytes(b'partial')
 assert not completion_matches(tmp_path,{'role':'test'},core['events'][0])
 (tmp_path/'COMPLETE.json').write_text('{')
 assert not completion_matches(tmp_path,{'role':'test'},core['events'][0])
