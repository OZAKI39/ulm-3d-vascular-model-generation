from pathlib import Path
import sys,json,multiprocessing
import pytest
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *

@pytest.fixture(scope='session')
def root():return ROOT

@pytest.fixture(scope='session')
def population():return load_births(ROOT)

@pytest.fixture(scope='session')
def core():return json.loads((ROOT/REL/'data/CORE500_COHORT.json').read_text())

@pytest.fixture(scope='session')
def real_env():
 from particle_3d.population_inlet_p9a4 import load_new_environment
 return load_new_environment(ROOT)

@pytest.fixture(scope='session')
def parallel_subset(tmp_path_factory,core,real_env):
 """Two actual, frozen NEW-field trajectories; independent fresh processes."""
 from concurrent.futures import ProcessPoolExecutor
 import particle_3d.formal_dynamics_p9a5 as d
 folder=tmp_path_factory.mktemp('p9a5_worker_parity')
 d.ENV=real_env;d.IDENTITY={'role':'PERMANENT_TEST_ONLY_NOT_FORMAL_PRODUCTION'}
 results=[]
 for workers in [1,2]:
  target=folder/str(workers)
  with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork')) as pool:
   rows=list(pool.map(d.job,[(e,str(target/str(e['particle_id']))) for e in core['events'][:2]]))
  receipts=[json.loads((target/str(e['particle_id'])/'receipt.json').read_text()) for e in core['events'][:2]]
  results.append(dict(rows=rows,hashes=[r['scientific_result_sha256'] for r in receipts],target=target))
 return results

@pytest.fixture
def valid_completion(tmp_path,core):
 identity={'source':'test'};event=core['events'][0];files={}
 for name in ['trajectory.npz','trajectory.json','audit.jsonl.gz','metrics.json','support.json']:
  p=tmp_path/name;p.write_bytes(b'test science bytes');files[name]=digest(p)
 write_new(tmp_path/'COMPLETE.json',dict(identity=identity,event_sha256=content_sha(event),files=files))
 return tmp_path,identity,event
