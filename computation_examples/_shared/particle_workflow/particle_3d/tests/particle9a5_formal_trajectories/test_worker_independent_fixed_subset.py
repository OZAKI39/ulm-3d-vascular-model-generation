from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_worker_independent_fixed_subset(parallel_subset,root):
 a,b=parallel_subset
 assert a['hashes']==b['hashes'] and len(a['hashes'])==2
 assert all(r['status']=='COMPLETED' for r in a['rows']+b['rows'])
 p=root/REL/'data/worker_scaling_results.json'
 if p.exists():
  evidence=json.loads(p.read_text());configs=evidence['configurations']
  assert [r['workers'] for r in configs]==[1,2,4,6,8,12,16]
  assert all(r['count']==24 and r['solver_failures']==0 for r in configs)
  assert len({r['scientific_results_sha256'] for r in configs})==1
  assert evidence['exact_same_dt_reference_24_trajectory_parity']
