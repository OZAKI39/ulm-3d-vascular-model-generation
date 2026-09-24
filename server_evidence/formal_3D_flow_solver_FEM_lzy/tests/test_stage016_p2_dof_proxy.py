from stage016_helpers import *
import pytest
from fem3d.planar_port import p2_proxy
from fem3d.audit import sha256

def test_proxy_unique_edges_with_shared_face_and_sparse_numbering():
 tetra=np.array([[10,20,30,40],[10,30,20,50]])
 expected={(min(a,b),max(a,b)) for t in tetra for i,a in enumerate(t) for b in t[i+1:]}
 p=p2_proxy(tetra)
 assert p['N_vertex']==5 and p['N_edge']==len(expected)==9
 assert p['N_P2_scalar_proxy']==14 and p['N_P2_velocity_proxy']==42 and p['N_P1_pressure_proxy']==5

def test_baseline_proxy_is_recomputed_and_mesh_hash_locked():
 actual=baseline_volume_audit()['proxy']
 assert actual==POLICY['cost']['baseline']
 assert actual['N_vertex']==42968 and actual['N_edge']==224167 and actual['N_P2_velocity_proxy']==801405
 assert sha256(ROOT/'outputs/stage01/medium/mesh/volume_mesh.npz')==POLICY['cost']['baseline_mesh_sha256']
