from pathlib import Path
import importlib.util,json
import numpy as np
import pytest
from particle_3d.wall_geometry import WallGeometry
from particle_3d.open_boundary_rim import build_rim_topology
from particle_3d.particle_shapes import Sphere
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle9a_motion import augment_planar_system
from particle_3d.wall_gap import wall_gap
ROOT=Path(__file__).resolve().parents[2]
@pytest.fixture
def old_planar():
 p=ROOT/'reports/particle9a1_2mmps/reference/old_planar_wall_hydrodynamics.py'
 spec=importlib.util.spec_from_file_location('old_planar',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
@pytest.fixture
def id7():return json.loads((ROOT/'reports/particle9a1_2mmps/reference/id7_duplicate_contact.json').read_text())
@pytest.fixture
def rim_case():
 points=1e-5*np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0.]])
 faces=np.array([[0,1,2],[1,3,2]]);w=WallGeometry(points[faces]);w.global_node_ids=faces
 # Two caps jointly form the same square perimeter; diagonal is internal.
 w.open_boundary_topology=build_rim_topology(faces,{'INLET':faces})
 def build(x):
  shapes={1:Sphere(x,1e-7)};free={1:np.array([.002,0,0,0,0,0.])}
  base=assemble_v1(shapes,free,.00345312,w)
  new=augment_planar_system(base,shapes,.00345312,w,lambda x:np.zeros((3,3)),wall_gap)
  return base,new,wall_gap(shapes[1],w)
 return build
