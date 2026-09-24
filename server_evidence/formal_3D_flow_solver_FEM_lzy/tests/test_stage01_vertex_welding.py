from pathlib import Path
import numpy as np
import pytest
from fem3d.mesh_input import exact_weld,adapt_surface

ROOT=Path(__file__).resolve().parents[1]

def test_actual_welding_is_deterministic_and_keeps_every_triangle():
    a=adapt_surface(ROOT);b=adapt_surface(ROOT)
    for x,y in zip(a[:3],b[:3]):np.testing.assert_array_equal(x,y)
    assert a[3]['welding']==b[3]['welding']
    w=a[3]['welding']
    assert w['points_before_welding']==73416
    assert w['triangle_count_before']==w['triangle_count_after']==67262
    assert w['maximum_merge_distance_m']==w['welding_tolerance_m']==0
    p,t=a[:2]
    assert (np.linalg.norm(np.cross(p[t[:,1]]-p[t[:,0]],p[t[:,2]]-p[t[:,0]]),axis=1)>0).all()

def test_near_but_distinct_points_are_not_merged():
    p=np.array([[1e-4,0,0],[np.nextafter(1e-4,np.inf),0,0],[1e-4,1e-6,0],[1e-4,0,1e-6],[1e-4,0,0]])
    t=np.array([[0,2,3],[1,2,3],[4,2,3]])
    u,w,m=exact_weld(p,t)
    assert len(u)==4 and w[0,0]!=w[1,0] and w[0,0]==w[2,0]
    np.testing.assert_array_equal(u[w],p[t])
    for tol in (1e-15,1e-6,1e-5):
        with pytest.raises(ValueError,match='Tolerance welding'):exact_weld(p,t,tol)

def test_degenerate_and_nonfinite_sources_are_rejected():
    with pytest.raises(ValueError,match='Degenerate'):exact_weld([[0,0,0],[0,0,0],[0,1,0]],[[0,1,2]])
    with pytest.raises(ValueError,match='coordinates'):exact_weld([[0,0,0],[np.nan,0,0],[0,1,0]],[[0,1,2]])
