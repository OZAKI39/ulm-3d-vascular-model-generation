from itertools import permutations
import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere,Ellipsoid,Capsule,EPS
from particle_3d.pair_geometry import pair_gap
from particle_3d.pair_broadphase import broadphase,all_pairs
from particle_3d.particle4_motion import initial_world,InitialParticleOverlap,simulate_world,pair_interval_certificate
from particle_3d.particle4_cases import shapes_for_geometry,mixed_scene


@pytest.mark.parametrize('dt',[1.,.5,.25])
def test_pair_no_tunneling_and_actual_time_coverage(dt):
    shapes={3:Sphere([-2e-6,0,0],.5e-6),7:Sphere([2e-6,0,0],.5e-6)}
    velocity=lambda i,s,t:(np.array([4e-6 if i==3 else -4e-6,0,0]),np.zeros(3))
    final,rows,ledger=simulate_world(shapes,velocity,dt,1.)
    assert final.time_s==1.
    assert abs(sum(r['dt_s'] for r in ledger)-1.)<=16*EPS
    touching=[r for r in rows if r['pair_gaps'][0]['state']=='TOUCHING']
    assert touching and abs(touching[0]['time_s']-.375)<=16*EPS
    assert all(g['gap_m']>=-g['roundoff_budget_m'] for r in rows for g in r['pair_gaps'])
    assert all(a['t1_s']==b['t0_s'] for a,b in zip(ledger,ledger[1:]))


def test_initial_overlap_rejected_without_moving_particles():
    a=Sphere([0,0,0],1e-6);b=Sphere([1e-6,0,0],1e-6)
    with pytest.raises(InitialParticleOverlap) as error:initial_world({1:a,2:b},lambda *args:(np.zeros(3),np.zeros(3)))
    assert error.value.record['status']=='INITIAL_PARTICLE_OVERLAP'
    assert error.value.record['pair']['gap_m']<0
    np.testing.assert_array_equal(b.center_m,[1e-6,0,0])


def test_glancing_spheres_integrate_contact_ode_without_chattering_or_projection():
    shapes={1:Sphere([-1e-6,0,0],1e-6),2:Sphere([1e-6,0,0],1e-6)}
    free={1:np.array([1.,.5,0])*1e-6,2:np.array([-1.,-.5,0])*1e-6}
    provider=lambda i,s,t:(free[i],np.zeros(3))
    final,rows,ledger=simulate_world(shapes,provider,.125,.5)
    assert len(ledger)==4 and final.time_s==.5
    for row in rows:
        g=row['pair_gaps'][0];assert abs(g['gap_m'])<=g['roundoff_budget_m']
        centers=[np.array(p['center_m']) for p in row['particles']]
        assert np.linalg.norm(sum(centers))<=g['roundoff_budget_m']
    # Independent high-accuracy ODE integration of relative center motion.
    from scipy.integrate import solve_ivp
    w=free[1]-free[2]
    def rhs(t,r):
        n=r/np.linalg.norm(r);return w-min(w@n,0.)*n
    reference=solve_ivp(rhs,[0,.5],shapes[1].center_m-shapes[2].center_m,rtol=128*EPS,atol=128*EPS*2e-6)
    actual=final.shapes[1].center_m-final.shapes[2].center_m
    assert np.linalg.norm(actual-reference.y[:,-1])<=4096*EPS*2e-6


def test_broadphase_no_false_negative_and_all_pairs_equivalence():
    rng=np.random.default_rng(2026092004);base=shapes_for_geometry()
    shapes={100+i:base[i%3].moved(rng.normal(size=3)*3e-6) for i in range(18)}
    exhaustive={ij:pair_gap(shapes[ij[0]],shapes[ij[1]],*ij) for ij in all_pairs(shapes)}
    candidates=set(broadphase(shapes));contacts={ij for ij,g in exhaustive.items() if g.state!='SEPARATED'}
    assert contacts<=candidates
    found={ij:pair_gap(shapes[ij[0]],shapes[ij[1]],*ij) for ij in sorted(candidates)}
    assert {ij for ij,g in found.items() if g.state!='SEPARATED'}==contacts
    for ij in contacts:
        assert found[ij].gap_m==exhaustive[ij].gap_m
        np.testing.assert_array_equal(found[ij].normal_j_to_i,exhaustive[ij].normal_j_to_i)


def test_swept_broadphase_includes_endpoint_separated_collision():
    shapes={1:Sphere([-2e-6,0,0],.5e-6),2:Sphere([2e-6,0,0],.5e-6)}
    end={i:s.moved(-s.center_m) for i,s in shapes.items()}
    assert broadphase(shapes)==[] and broadphase(end)==[]
    assert broadphase(shapes,end_shapes=end)==[(1,2)]
    clear,_=pair_interval_certificate(shapes[1],shapes[2],end[1],end[2],np.zeros(3),np.zeros(3),1.)
    assert not clear


def test_mode_change_and_capsule_axis_change_cannot_reuse_old_bound():
    base=shapes_for_geometry();a=base[0];b=base[1].moved([20e-6,0,0]);c=base[2].moved(b.center_m)
    assert not pair_interval_certificate(a,b,a,c,np.zeros(3),np.zeros(3),1.)[0]
    changed=Capsule(c.center_m,[1,0,0],c.radius_m,c.cylindrical_length_m)
    assert not pair_interval_certificate(a,c,a,changed,np.zeros(3),np.zeros(3),1.)[0]


def test_refinement_failure_is_explicit_and_does_not_consume_time():
    from particle_3d.physical_time_refinement import refine_interval,PhysicalTimeRefinementError
    from particle_3d.particle4_motion import world_trial
    shapes={1:Sphere([-2e-6,0,0],.5e-6),2:Sphere([2e-6,0,0],.5e-6)}
    provider=lambda i,s,t:(np.array([4e-6 if i==1 else -4e-6,0,0]),np.zeros(3))
    state=initial_world(shapes,provider);ledger=[]
    with pytest.raises(PhysicalTimeRefinementError) as error:refine_interval(state,1.,world_trial(provider),max_depth=0,ledger=ledger)
    assert ledger==[] and state.time_s==0.
    assert error.value.record['reason']=='DEPTH_LIMIT' and error.value.record['interval']==[0.,1.]


def test_rotating_pair_certificate_agrees_with_dense_independent_path_samples():
    from particle_3d.particle4_motion import move_shape
    from particle_3d.particle4_cases import place_pair
    a,b=place_pair(shapes_for_geometry()[1],shapes_for_geometry()[0],[1,.3,.4],.5e-6)
    w=np.array([.1,-.2,.3]);v=np.array([.02,-.01,.03])*1e-6;dt=.1
    end=move_shape(a,v,w,dt)
    clear,_=pair_interval_certificate(a,b,end,b,w,np.zeros(3),dt)
    assert clear
    for t in np.linspace(0,dt,41):
        current=move_shape(a,v,w,t) if t else a;g=pair_gap(current,b)
        assert g.gap_m>=-g.roundoff_budget_m


def test_capsule_translation_preserves_axis_bits_and_dimensions_across_contact_steps():
    from particle_3d.particle4_motion import move_shape
    from particle_3d.particle4_cases import place_pair
    rng=np.random.default_rng(9)
    for _ in range(30):
        shape=Capsule([0,0,0],rng.normal(size=3),.4e-6,1e-6)
        moved=move_shape(shape,np.array([1e-6,0,0]),np.ones(3),.1)
        np.testing.assert_array_equal(moved.axis_world,shape.axis_world)
    a,b=place_pair(shape,Sphere([0,0,0],.5e-6),[1,.4,-.1]);n=pair_gap(a,b).normal_j_to_i
    common=np.array([.2,.3,.1])*1e-6
    provider=lambda i,s,t:(common+(-n if i==1 else n)*1e-6,np.array([1.,2,3]))
    final,rows,ledger=simulate_world({1:a,2:b},provider,.025,.1)
    np.testing.assert_array_equal(final.shapes[1].axis_world,a.axis_world)
    assert final.shapes[1].radius_m==a.radius_m and final.shapes[1].cylindrical_length_m==a.cylindrical_length_m
    np.testing.assert_array_equal(final.omegas[1],[1.,2,3])


@pytest.mark.parametrize('dt',[.2,.1,.05])
def test_mixed_four_particle_trajectory(dt):
    shapes,provider=mixed_scene();final,rows,ledger=simulate_world(shapes,provider,dt,1.)
    assert final.time_s==1 and abs(sum(r['dt_s'] for r in ledger)-1)<=16*EPS
    contact_types=set()
    for row in rows:
        for gap in row['pair_gaps']:
            assert gap['gap_m']>=-gap['roundoff_budget_m']
            if gap['state']=='TOUCHING':contact_types.add(tuple(sorted([gap['shape_i'],gap['shape_j']])))
    assert contact_types=={('SPHERE_MB','SPHERE_MB'),('FREE_OBLATE','SPHERE_MB'),('FREE_OBLATE','FREE_OBLATE')}


def test_particle_permutation_preserves_trajectory_positions_and_contact_set():
    shapes,provider=mixed_scene();baseline=None
    for order in [list(shapes),list(shapes)[::-1],[203,409,101,307]]:
        final,rows,ledger=simulate_world({i:shapes[i] for i in order},provider,.2,.5)
        result=np.concatenate([np.r_[final.shapes[i].center_m,final.velocities[i],final.omegas[i]] for i in sorted(shapes)])
        contacts=[g.canonical_pair_id for g in final.pair_gaps if g.state=='TOUCHING']
        if baseline is None:baseline=result,contacts
        else:np.testing.assert_array_equal(result,baseline[0]);assert contacts==baseline[1]
