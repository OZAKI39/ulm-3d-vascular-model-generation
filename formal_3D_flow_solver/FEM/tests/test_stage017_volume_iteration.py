from stage017_helpers import *
import pytest
from fem3d.adaptive_port import volume_acceptance

def test_actual_volume_generation_and_DOLFINx_conversion():
 for r in result()['iterations']:
  p=OUT/r['iteration'];m=read(p/'metadata/meshing.json')
  assert m['effective_gmsh_options']==BASE['effective_gmsh_options'] and not m['gpu_used']
  assert read(p/'qc/dolfinx_loadable.json')['status']=='PASS'
  assert (p/'mesh/fluid.msh').exists() and r['surface_status']=='PASS'

def test_selected_validity_quality_and_relative_cost_recomputed():
 measured=selected_recomputed();a=volume_acceptance(measured,BASE,POLICY,True)
 assert a['status']=='PASS' and a['quality_pass'] and a['cost_pass']
 assert not any(measured['validity'].values())
 assert measured['topology']['connected_fluid_components']==1
