from stage017_helpers import *
import pytest
from fem3d.audit import sha256

def test_pre_search_immutable_policy_and_direct_v2_reference():
 lock=read(REPORT/'freeze_lock.json')
 assert sha256(REPORT/'acceptance_policy.json')==lock['policy_sha256']
 assert sha256(ROOT/'reports/stage01_6/planar_port_contract_v2.json')==lock['contract_sha256']
 assert POLICY['search']['search_factor']==1.25 and POLICY['search']['grading_slope']==1.
 assert POLICY['limits']['maximum_surface_trials_per_port']==8
 assert POLICY['limits']['maximum_volume_iterations']==POLICY['limits']['maximum_total_volume_meshes']==5
 assert POLICY['cost']=={'C_P2_max':1.35,'C_tetra_max':1.35,'cap_triangle_count_is_hard_gate':False}
 for path in (OUT/'surface_trials').glob('*/*/trial.json'):
  trial=read(path);assert lock['timestamp']<trial['timestamp'] and trial['policy_sha256']==lock['policy_sha256']

def test_production_optimizer_has_no_model_specific_cap_lengths():
 import re,ast
 for p in (ROOT/'src/fem3d').glob('adaptive_*.py'):
  text=p.read_text()
  values=[node.value for node in ast.walk(ast.parse(text)) if isinstance(node,ast.Constant) and isinstance(node.value,float)]
  assert not any(v in (3.5e-7,4.5e-7,5.5e-7,6.5e-7) for v in values)
  assert not re.search(r'0\.(35|45|55|65)\s*e-6',text)
  assert not re.search(r'(3\.5|4\.5|5\.5|6\.5)\s*e-7',text)
  assert 'radius_um' not in text and '800' not in text and '147569' not in text
