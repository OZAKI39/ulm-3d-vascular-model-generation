from stage016_helpers import *
import pytest
from fem3d.cap_remesh import check_rim,rim_loop

@pytest.mark.parametrize('name',CANDIDATES)
def test_frozen_3d_rim_ids_edges_and_coordinates(name):
 for r in audit_candidate(name).values():assert r['rim']['maximum_rim_displacement_m']==0 and r['rim']['rim_edge_set_identical']

def test_malicious_one_rim_coordinate_moved():
 p,t=square();_,edges=rim_loop(t);m=p.copy();m[0,2]=1e-15
 with pytest.raises(ValueError,match='Rim vertex moved'):check_rim(p,edges,m,t)

def test_malicious_rim_edge_split_rejected():
 p,t=square();_,edges=rim_loop(t);m=np.vstack([p,[.5,0,0]]);new=np.vstack([[[0,5,4],[5,1,4]],t[1:]])
 with pytest.raises(ValueError,match='Rim edge set changed'):check_rim(p,edges,m,new)
