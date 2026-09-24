from stage016_helpers import *
import pytest
from fem3d.audit import sha256
from PIL import Image

def test_all_ten_figures_have_actual_evidence_and_explicit_missing_panels():
 m=read(REPORT/'visualization_manifest.json');assert len(m['images'])==10 and m['selected_candidate'] is None
 expected={'port_contract_explanation','cap_triangulation_before_after','cap_density_quality_tradeoff','tetra_quality_before_after','low_quality_count_by_boundary','worst_elements_before_after','rim_overlay','boundary_tags_selected','tetrahedral_cutaway_selected','mesh_cost_comparison'}
 assert {Path(p).stem for p in m['images']}==expected
 for name,data in m['images'].items():
  assert sha256(REPORT/name)==data['sha256']
  with Image.open(REPORT/name) as im:assert im.width>=1000 and im.height>=700
 assert 'unavailable' in m['limitation'] and '10000x' in m['images']['port_contract_explanation.png']['role']
