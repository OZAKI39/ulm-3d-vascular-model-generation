import json
from pathlib import Path
from functools import lru_cache
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPORT=ROOT/'reports/stage01_6';OUT=ROOT/'outputs/stage01_6'
CANDIDATES=('sparse_A','sparse_B','sparse_C')
def read(p): return json.loads(p.read_text())
POLICY=read(REPORT/'acceptance_policy.json');CONTRACT=read(REPORT/'planar_port_contract_v2.json')
@lru_cache(None)
def source(): return dict(np.load(ROOT/'inputs/stage01/tagged_surface_si.npz'))
@lru_cache(None)
def candidate(name): return dict(np.load(OUT/name/'surface/tagged_surface_si.npz'))
def square():
 return np.array([[0.,0.,0.],[1,0,0],[1,1,0],[0,1,0],[.5,.5,0]]),np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4]])
@lru_cache(None)
def audit_candidate(name):
 from fem3d.planar_port import validate_port
 d=candidate(name); results={}
 for n,p in CONTRACT['ports'].items():
  tri=d['triangles'][d['facet_tags']==p['entity_id']]
  ids=np.setdiff1d(np.unique(tri),p['rim_vertex_ids'])
  results[n]=validate_port(source()['points_m'],d['points_m'],tri,p,ids,POLICY,origin=p['plane_origin_m'],normal=p['outward_normal'],basis=p['basis'])
 return results
@lru_cache(None)
def baseline_volume_audit():
 from fem3d.sparse_volume import audit_volume
 return audit_volume(dict(np.load(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz')),source(),source(),CONTRACT,POLICY)
def good_volume():
 return {'zero_volume':0,'negative_volume':0,'nonfinite_volume':0},{'min_sicn':{'P1':.5,'P5':.6,'median':.8},'cap_adjacent_below_0_1':10,'total_below_0_1':20},dict(POLICY['cost']['baseline'])
