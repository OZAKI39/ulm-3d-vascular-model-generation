from stage017_helpers import *
import pytest
from fem3d.audit import sha256
from PIL import Image

def test_twelve_actual_figures_and_hashes():
 manifest=read(REPORT/'visualization_manifest.json');assert len(manifest['images'])==12
 assert manifest['selected_iteration']==result()['selected_iteration']
 expected={'adaptive_surface_search','port_mesh_before_after','surface_quality_before_after','adaptive_volume_trace','tetra_quality_before_after','low_quality_count_by_boundary','quality_cost_pareto','mesh_cost_comparison','worst_elements_before_after','rim_overlay','boundary_tags_selected','tetrahedral_cutaway_selected'}
 assert {Path(p).stem for p in manifest['images']}==expected
 for name,data in manifest['images'].items():
  assert sha256(REPORT/name)==data['sha256']
  with Image.open(REPORT/name) as im:assert im.width>=1000 and im.height>=700
