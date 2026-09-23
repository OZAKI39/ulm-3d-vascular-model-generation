from stage016_helpers import *
import pytest
def test_no_rejected_surface_enters_volume_generation():
 for name in CANDIDATES:
  s=read(OUT/name/'qc/surface_invariants.json')
  assert not s['volume_meshing_permitted']
  assert not (OUT/name/'mesh/volume_mesh.npz').exists() and not (OUT/name/'mesh/fluid.msh').exists()
 assert not (OUT/'selected').exists()

def test_geometry_audit_reproduces_frozen_baseline():
 r=baseline_volume_audit()
 assert r['boundary_fidelity']['maximum_boundary_displacement_m']==0
 assert r['validity']=={'zero_volume':0,'negative_volume':0,'nonfinite_volume':0}
 assert r['topology']['connected_fluid_components']==1
 assert all(p['status']=='PASS' for p in r['ports'].values())

def test_volume_controls_exactly_match_stage1_medium():
 original=read(ROOT/'outputs/stage01/medium/metadata/meshing.json')['effective_gmsh_options']
 assert POLICY['volume_meshing']['effective_gmsh_options']==original
 assert original['Mesh.Algorithm3D']==1 and original['Mesh.MeshSizeMax']==5e-7 and original['Mesh.OptimizeNetgen']==0
