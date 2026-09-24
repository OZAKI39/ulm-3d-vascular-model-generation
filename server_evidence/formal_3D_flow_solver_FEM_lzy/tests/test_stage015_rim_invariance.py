import numpy as np
import pytest
from fem3d.cap_remesh import check_rim,rim_loop
from stage015_helpers import square


def test_original_rim_is_bitwise_unchanged():
    p,t=square();_,edges=rim_loop(t)
    assert check_rim(p,edges,p.copy(),t)['maximum_rim_displacement_m']==0


def test_moving_one_rim_vertex_fails():
    p,t=square();_,e=rim_loop(t);changed=p.copy();changed[0,0]=1e-15
    with pytest.raises(ValueError,match='moved'): check_rim(p,e,changed,t)


def test_splitting_an_original_rim_edge_fails_even_with_exact_coverage():
    p,t=square();_,e=rim_loop(t);new=np.vstack([p,[.5,0,0]])
    tri=np.vstack([[[0,5,4],[5,1,4]],t[1:]])
    with pytest.raises(ValueError,match='edge set'): check_rim(p,e,new,tri)


def test_deleting_a_rim_edge_fails():
    p,t=square();_,e=rim_loop(t)
    with pytest.raises(ValueError): check_rim(p,e,p,t[1:])
