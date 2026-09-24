from stage017_helpers import *
import pytest
from fem3d.planar_port import p2_proxy
from fem3d.adaptive_port import volume_acceptance

def test_unique_edges_shared_face_and_renumbering():
 p=p2_proxy(np.array([[10,20,30,40],[20,10,30,50]]))
 assert p['N_vertex']==5 and p['N_edge']==9 and p['N_P2_velocity_proxy']==42

def test_relative_cost_is_measured_not_cap_count():
 m=selected_recomputed();a=volume_acceptance(m,BASE,POLICY,True)
 assert a['C_P2']==m['proxy']['N_P2_velocity_proxy']/BASE['proxy']['N_P2_velocity_proxy']
 assert a['C_tetra']==m['proxy']['N_tetra']/BASE['proxy']['N_tetra']
