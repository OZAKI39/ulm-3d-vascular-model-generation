from dataclasses import dataclass
import numpy as np
import pytest
from particle_3d.particle1_cases import AffineValidationField
from particle_3d.particle3_motion import initial_wall_state,contact_trial
from particle_3d.physical_time_refinement import refine_interval,PhysicalTimeRefinementError,TrialNeedsSubdivision,capsule_plane_minimum,swept_clearance_certificate
from particle_3d.particle_shapes import Capsule,Sphere


@pytest.mark.parametrize('pieces',[1,2,4])
def test_no_tunneling_and_complete_physical_time(plane_wall,pieces):
    field=AffineValidationField(np.array([1e-6,0,-4e-6]),np.zeros((3,3)))
    state=initial_wall_state([0,0,2e-6],[1,0,0,0],field,plane_wall,radius=.5e-6)
    ledger=[];trial=contact_trial(field,plane_wall)
    for i in range(pieces):state,_=refine_interval(state,(i+1)/pieces,trial,ledger=ledger)
    assert state.time_s==1.
    np.testing.assert_allclose(state.shape.center_m,[1e-6,0,.5e-6],atol=2e-17,rtol=0)
    assert abs(sum(r['dt_s'] for r in ledger)-1)<32*np.finfo(float).eps
    assert all(a['t1_s']==b['t0_s'] for a,b in zip(ledger[:-1],ledger[1:]))
    assert any(r['depth']>0 for r in ledger)
    clear,_,_=swept_clearance_certificate(Sphere([0,0,2e-6],.5e-6),Sphere([1e-6,0,-2e-6],.5e-6),plane_wall)
    assert not clear  # both endpoint gaps positive, but wall lies between them


def test_refinement_guard_reports_unconsumed_interval_and_state():
    @dataclass
    class State:
        time_s:float=0.
        shape_mode:str='SPHERE_MB'
    def fails(state,dt):raise TrialNeedsSubdivision('deliberate software guard test',-1e-6,3)
    s=State()
    with pytest.raises(PhysicalTimeRefinementError) as e:refine_interval(s,1.,fails,max_depth=3)
    assert s.time_s==0 and e.value.record['interval']==[0.,.125]
    assert e.value.record['reason']=='DEPTH_LIMIT' and e.value.record['triangle_id']==3


def test_capsule_swept_plane_checks_interior_axis_maximum(plane_triangle):
    # Endpoints are below the middle support height: checking endpoints alone fails.
    a=Capsule([0,0,2e-6],[-1,0,2],.4e-6,3.3e-6)
    b=Capsule([0,0,2e-6],[1,0,2],.4e-6,3.3e-6)
    value=capsule_plane_minimum(a,b,[0,0,1],plane_triangle)
    assert abs(value-(-.05e-6))<1e-18
    assert value<0


def test_continuous_capsule_stationary_minimum_against_dense_independent_paths(plane_triangle):
    rng=np.random.default_rng(335)
    for _ in range(40):
        axis=rng.normal(size=3);other=axis+.5*rng.normal(size=3)
        a=Capsule(rng.normal(size=3)*1e-6,axis,.5e-6,3e-6)
        radius=rng.uniform(.4,.6)*1e-6
        b=Capsule(a.center_m+rng.normal(size=3)*.3e-6,other,radius,a.volume_m3/(np.pi*radius**2)-4*radius/3)
        normal=np.array([0,0,1.]);minimum=capsule_plane_minimum(a,b,normal,plane_triangle)
        t=np.linspace(0,1,4001);r=a.radius_m+t*(b.radius_m-a.radius_m)
        directions=a.axis_world+t[:,None]*(b.axis_world-a.axis_world);directions/=np.linalg.norm(directions,axis=1)[:,None]
        centers=a.center_m+t[:,None]*(b.center_m-a.center_m)
        length=a.volume_m3/(np.pi*r*r)-4*r/3
        dense=centers@normal-r-length/2*np.abs(directions@normal)
        assert minimum<=dense.min()+512*np.finfo(float).eps*max(a.bounding_radius_m,b.bounding_radius_m)
