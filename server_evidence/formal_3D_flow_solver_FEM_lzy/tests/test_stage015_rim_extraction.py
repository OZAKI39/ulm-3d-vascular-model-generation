import numpy as np
import pytest
from fem3d.cap_remesh import rim_loop


def test_single_loop_has_exact_boundary_edges():
    loop,edges=rim_loop([[0,1,2],[0,2,3]])
    assert set(loop)=={0,1,2,3}
    assert {tuple(e) for e in edges}=={(0,1),(1,2),(2,3),(0,3)}


@pytest.mark.parametrize('triangles',[
    [[0,1,2],[3,4,5]],
    [[0,1,2],[0,3,4]],
    [[0,1,2],[1,0,3],[0,1,4]],
    [[0,1,1]],
])
def test_multiple_branching_nonmanifold_or_degenerate_rims_fail(triangles):
    with pytest.raises(ValueError): rim_loop(np.asarray(triangles))
