import numpy as np
import pyvista as pv
from particle_3d.validation_boundary import ValidationBoundaryClassifier


def plane(x):
    return pv.PolyData(np.array([[x,-1.,-1.],[x,1.,-1.],[x,0.,1.]]),np.array([3,0,1,2]))


def test_first_surface_on_whole_segment_even_with_far_endpoint():
    classifier=ValidationBoundaryClassifier({'WALL':plane(.4),'OUTLET_01':plane(.8),'INLET':plane(-1.)})
    hit=classifier.first_event([0.,0.,0.],[2.,0.,0.])
    assert hit.role=='WALL' and hit.segment_fraction==.2
    np.testing.assert_array_equal(hit.position_m,[.4,0.,0.])
    assert classifier.first_event([0.,2.,0.],[2.,2.,0.]) is None


def test_exact_endpoint_and_shared_rim_tie_are_deterministic():
    classifier=ValidationBoundaryClassifier({'WALL':plane(1.),'OUTLET_01':plane(1.)})
    for _ in range(3):
        hit=classifier.first_event([0.,0.,0.],[1.,0.,0.])
        assert hit.role=='WALL' and hit.segment_fraction==1.
        assert set(hit.simultaneous_roles)=={'WALL','OUTLET_01'}


def test_coplanar_contact_and_stationary_segment():
    classifier=ValidationBoundaryClassifier({'WALL':plane(1.)})
    hit=classifier.first_event([1.,-2.,0.],[1.,0.,0.])
    assert hit.role=='WALL'
    np.testing.assert_allclose(hit.position_m,[1.,-.5,0.],rtol=0,atol=32*np.finfo(float).eps)
    assert classifier.first_event([0.,0.,0.],[0.,0.,0.]) is None


def test_wall_event_stops_validation_as_fail_without_boundary_response():
    from particle_3d.field import FrozenFEMField
    from particle_3d.particle1_cases import real_trajectory
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    field=FrozenFEMField(points,np.array([[0,1,2,3]]),np.tile([1.,0.,0.],(4,1)),np.zeros(4))
    wall=pv.PolyData(points[[1,2,3]],np.array([3,0,1,2]))
    classifier=ValidationBoundaryClassifier({'WALL':wall})
    initial={'initial_position_m':[.1,.1,.1],'horizon_s':2.}
    rows,summary=real_trajectory(field,classifier,initial,1e-6,1.)
    assert summary['passed'] is False and summary['wall_crossing'] is True
    assert len(rows)==2 and rows[-1]['boundary_event']=='WALL'
    np.testing.assert_allclose(summary['event']['unmodified_trial_endpoint_m'],[1.1,.1,.1],atol=0,rtol=0)
    np.testing.assert_allclose(summary['event']['position_m'],[.8,.1,.1],atol=8*np.finfo(float).eps,rtol=0)
    assert initial['initial_position_m']==[.1,.1,.1]
