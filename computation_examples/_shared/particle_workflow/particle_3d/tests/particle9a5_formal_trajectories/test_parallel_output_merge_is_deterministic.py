from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_parallel_output_merge_is_deterministic():
 rows=[dict(particle_id=i,status='COMPLETED',outlet='O3') for i in range(1,25)]
 expected=content_sha(merge_rows(rows,range(1,25)))
 for seed in range(6):
  __import__('random').Random(seed).shuffle(rows)
  assert content_sha(merge_rows(rows,range(1,25)))==expected
 with pytest.raises(ValueError,match='Missing'):merge_rows(rows[:-1],range(1,25))
 with pytest.raises(ValueError):merge_rows(rows+[rows[0]],range(1,25))
