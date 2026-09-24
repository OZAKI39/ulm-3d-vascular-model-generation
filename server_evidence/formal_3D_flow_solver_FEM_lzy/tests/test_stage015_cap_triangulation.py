import numpy as np
import pytest
from fem3d.cap_remesh import coverage_check
from stage015_helpers import square


def test_exact_polygon_partition_passes():
    p,t=square();result=coverage_check(p[:,:2],t,np.arange(4))
    assert result['holes']==result['overlaps']==result['triangles_outside_polygon']==0


@pytest.mark.parametrize('bad',['outside','hole','overlap','flipped_normal'])
def test_adversarial_coverage_failures(bad):
    p,t=square()
    if bad=='outside': p[4]=[1.5,.5,0]
    elif bad=='hole': t=t[:-1]
    elif bad=='overlap': t=np.vstack([t,[0,1,2]])
    elif bad=='flipped_normal': t=t[:,[0,2,1]]
    with pytest.raises(ValueError): coverage_check(p[:,:2],t,np.arange(4))


def test_pairwise_clipping_detects_positive_overlap_and_ignores_shared_edge():
    from fem3d.cap_remesh import convex_intersection_area
    t=np.array([[0.,0.],[1,0],[0,1]])
    assert convex_intersection_area(t,t)==.5
    assert convex_intersection_area(t,np.array([[1.,0.],[1,1],[0,1]]))==0
