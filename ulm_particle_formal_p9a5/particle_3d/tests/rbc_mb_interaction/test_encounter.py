from pathlib import Path
import json,sys
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'particle_3d/src'))
sys.path.insert(0,str(ROOT/'particle_3d/scripts'))
from particle_3d.particle_shapes import Sphere,Capsule,unit
from particle_3d.rbc_mb_encounter import advance_encounter,capsule_sphere_gap
from particle_3d.pair_geometry import pair_gap
from particle_3d.rbc_capillary_surrogate import capsule_area
from simulate_rbc_mb_interaction import OUT,CASE,sha,initial_shapes,setup


def test_capsule_gap_agrees_with_existing_independent_convex_kernel():
    rng=np.random.default_rng(27)
    cap=Capsule([100e-6,80e-6,90e-6],unit([1,2,.3]),1.3e-6,6.4e-6)
    for _ in range(12):
        n=unit(rng.normal(size=3));bubble=Sphere(cap.support(n)+(0.588e-6+1e-8)*n,.588e-6)
        a=capsule_sphere_gap(cap,bubble);b=pair_gap(cap,bubble,101,203)
        assert abs(a.gap_m-b.gap_m)<2e-16
        np.testing.assert_allclose(a.normal_j_to_i,b.normal_j_to_i,rtol=0,atol=2e-8)


def test_exact_contact_matches_independent_ode_and_preserves_p4_metric():
    cap=Capsule([0,0,0],[1,0,0],1.3e-6,6.4e-6)
    direction=unit([-1,.5,.15]);radius=cap.radius_m+.588e-6;tip=cap.center_m-3.2e-6*cap.axis_world
    mb=Sphere(tip+radius*direction,.588e-6)
    va=np.array([.1,.03,.01])*1e-3;vb=np.array([.22,.07,.015])*1e-3;w=va-vb;dt=.001
    a,b,v,rec=advance_encounter(cap,mb,va,vb,dt)
    relative=tip-mb.center_m
    def rhs(t,r):
        n=r/np.linalg.norm(r)
        return w-min(float(w@n),0.)*n
    ref=solve_ivp(rhs,[0,dt],relative,rtol=2e-12,atol=1e-19)
    actual=(a.center_m-3.2e-6*cap.axis_world)-b.center_m
    assert np.linalg.norm(actual-ref.y[:,-1])<1e-16
    np.testing.assert_allclose(a.center_m+b.center_m,cap.center_m+mb.center_m+(va+vb)*dt,rtol=0,atol=1e-19)
    g=capsule_sphere_gap(a,b)
    assert abs(g.normal_j_to_i@(v[101]-v[203]))<1e-15
    assert rec['contact_duration_s']==dt and not rec['position_projection']


def test_free_flight_contact_event_and_unsupported_regime_rejection():
    import pytest
    a=Capsule([0,0,0],[1,0,0],1e-6,4e-6);b=Sphere([-3.6e-6,0,0],.5e-6)
    aa,bb,v,r=advance_encounter(a,b,[0,0,0],[1e-4,0,0],.002)
    assert abs(r['hit_time_s']-.001)<1e-14
    assert abs(capsule_sphere_gap(aa,bb).gap_m)<1e-18
    np.testing.assert_allclose(v[101],v[203],rtol=0,atol=1e-18)
    with pytest.raises(ValueError,match='HEMISPHERE'):
        advance_encounter(a,Sphere([0,2e-6,0],.5e-6),[0,0,0],[0,0,0],.001)


def saved(name):return json.loads((OUT/f'data/{name}.json').read_text())


def test_saved_run_time_geometry_wall_and_actual_interaction():
    states=saved('fine_states');ledger=saved('fine_ledger');init=saved('initialization');summary=saved('fine_summary')
    a,b=initial_shapes(init).values()
    assert np.isclose(a.volume_m3,init['rbc_original_geometry']['volume_m3'],rtol=1e-13,atol=0)
    assert a.area_m2<=init['area_budget_m2']
    assert sha(CASE/'frozen_flow/flow_arrays_si.npz')==summary['source_npz_sha256']
    assert sha(CASE/'SV_MESH/mesh-surfaces/WALL.vtp')==summary['source_wall_sha256']
    assert states[0]['time_s']==0 and states[-1]['time_s']==.01
    assert np.all(np.diff([s['time_s'] for s in states])>0)
    assert abs(sum(r['dt_s'] for r in ledger)-.01)<1e-15
    assert all(x['t1_s']==y['t0_s'] for x,y in zip(ledger,ledger[1:]))
    corrections=[]
    for s in states:
        x,y=s['particles'];gap=s['pair_gaps'][0]
        assert gap['gap_m']>=-gap['roundoff_budget_m']
        assert min(x['wall_gap_m'],y['wall_gap_m'])>0
        assert x['radius_m']==a.radius_m and x['cylindrical_length_m']==a.cylindrical_length_m
        np.testing.assert_array_equal(x['axis_world'],a.axis_world)
        if s['projection']:
            r=s['projection'];assert r['wall_clearance_at_start_m']>r['wall_motion_bound_m']
            assert r['hemisphere_margin_m']>r['relative_motion_bound_m']
            corrections.append(np.linalg.norm(np.array(y['velocity_m_s'])-r['mb_free_velocity_m_s']))
    assert max(corrections)>1e-5 and summary['contact_states']>100


def test_free_velocities_are_from_new_fem_and_halved_step_converges():
    fine=saved('fine_states');coarse=saved('coarse_states');field,wall,provider=setup()
    for k in np.linspace(0,len(fine)-2,24,dtype=int):
        for j,key in [(0,'cap_free_velocity_m_s'),(1,'mb_free_velocity_m_s')]:
            sample=field.sample(fine[k]['particles'][j]['center_m'])
            assert sample.inside_lumen
            np.testing.assert_allclose(sample.velocity_m_s,fine[k+1]['projection'][key],rtol=0,atol=1e-15)
    errors=[]
    for a,b in zip(coarse,fine[::2]):
        assert abs(a['time_s']-b['time_s'])<1e-15
        errors.extend(np.linalg.norm(np.array(x['center_m'])-y['center_m']) for x,y in zip(a['particles'],b['particles']))
    assert max(errors)<.01e-6


def test_exact_replay_returns_saved_endpoints_and_never_overlaps():
    from render_rbc_mb_interaction import Replay
    replay=Replay()
    for k in np.linspace(0,len(replay.times)-1,31,dtype=int):
        pos,_,_,g,_=replay.at(replay.times[k])
        np.testing.assert_allclose(pos,replay.centers[k],rtol=0,atol=2e-18)
        assert g.gap_m>=-g.roundoff_budget_m
    for t in np.linspace(0,replay.end,432):
        _,_,_,g,_=replay.at(t)
        assert g.gap_m>=-g.roundoff_budget_m
