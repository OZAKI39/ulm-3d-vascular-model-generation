from stage016_helpers import *
import pytest
from fem3d.audit import sha256

def test_pre_run_frozen_sizes_formula_and_hashes():
 lock=read(REPORT/'freeze_lock.json')
 assert POLICY['candidates']=={'sparse_A':3.5e-7,'sparse_B':4.5e-7,'sparse_C':5.5e-7}
 assert POLICY['sizing']['transition_R_eq']==.5
 assert sha256(REPORT/'acceptance_policy.json')==lock['policy_sha256']
 assert sha256(REPORT/'planar_port_contract_v2.json')==lock['contract_sha256']
 for name in CANDIDATES:
  m=read(OUT/name/'metadata/cap_remesh.json')
  assert lock['timestamp']<m['timestamp'] and m['config_sha256']==lock['policy_sha256']
  assert m['central_target_m']==POLICY['candidates'][name]
  assert m['sizing']==POLICY['sizing'] and m['gmsh_options']['Mesh.Algorithm']==6
  assert m['gpu_used'] is False

def test_no_candidate_search_expansion():
 assert not any(p.name.startswith('sparse_') and p.name not in CANDIDATES for p in OUT.iterdir())
 assert not read(REPORT/'candidate_selection.json')['parameter_search_expanded']
