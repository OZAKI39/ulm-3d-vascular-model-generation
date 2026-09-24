import numpy as np
import pytest
from particle_3d.particle65_cases import synthetic,records_for,MU
from particle_3d.particle65_motion import Particle65Stepper
from particle_3d.particle81_cache import cached_step_to
from particle_3d.lammps_neighbors import ValidationNeighborPolicy
from particle_3d.injection_population import FluxClock,LinearProfile,ConstantMBConcentrationV0,C_MB
from particle_3d.particle8_replay import read


@pytest.mark.parametrize('dt',[.0001,.0005,.001])
@pytest.mark.parametrize('global_bound',[False,True])
def test_private_query_adapter_matches_original_p65(dt,global_bound):
    shapes,provider,wall,_=synthetic('wall');particles=records_for(shapes,provider)
    policy=ValidationNeighborPolicy(20e-6,.5e-6,'P81_TEST_ONLY')
    original=Particle65Stepper(particles,policy,provider,MU,wall=wall)
    cached=Particle65Stepper(particles,policy,provider,MU,wall=wall)
    advance,info=cached_step_to(cached,global_motion_bound=global_bound)
    for k in range(1,41):
        original.step_to(k*dt);advance(k*dt)
        a=original.read()[0];b=cached.read()[0]
        np.testing.assert_array_equal(a.position,b.position)
        np.testing.assert_array_equal(a.velocity,b.velocity)
        np.testing.assert_array_equal(a.q,b.q)
        assert original.time_s==cached.time_s
    assert info().hits>0


def test_long_acquisition_clock_is_original_p7_cumulative_flux(root):
    b=read(root/'data/birth_ledger.json');clock=FluxClock(LinearProfile([0],[b['Q_in_m3_s']]),ConstantMBConcentrationV0())
    assert b['C_MB_m3']==8.5e12 and b['H_D_feed']==.45
    for i,event in enumerate(b['events'],1):
        assert event['particle_id']==i and event['birth_time_s']==clock.time_at(i)
        assert abs(clock.cumulative(event['birth_time_s'])-i)<=4*np.spacing(float(i))
        assert event['radius_m']==event['diameter_um']*.5e-6
    assert np.all(np.diff([e['birth_time_s'] for e in b['events']])>0)
    assert np.allclose(np.linalg.norm([e['q'] for e in b['events']],axis=1),1.,rtol=0,atol=4e-16)


def test_full_original_reference_not_just_inlet_smoke(root):
    for pid in [1,4]:
        a=np.load(root/f'data/reference_original_run/trajectories/mb_{pid:06d}.npz')['samples']
        b=np.load(root/f'trajectories/mb_{pid:06d}.npz')['samples']
        np.testing.assert_array_equal(a,b)
        assert a.tobytes()==b.tobytes()
        m=read(root/f'trajectories/mb_{pid:06d}.json')
        assert m['completed'] and m['path_length_m']>50e-6 and m['residence_time_s']>.05


def test_cache_v2_preserves_full_successful_and_stopped_physical_records(root):
    record=read(root/'data/cache_v2_parity.json')
    assert record['all_pass'] and {r['particle_id'] for r in record['cases']}=={1,4,7,13}
    for row in record['cases']:
        assert row['raw_arrays_bitwise_identical'] and row['terminal_same'] and row['provider_calls_same']
        pid=row['particle_id']
        actual=np.load(root/f'data/cache_v2_benchmark/trajectories/mb_{pid:06d}.npz')['samples']
        reference=np.load(root/f'trajectories/mb_{pid:06d}.npz')['samples']
        np.testing.assert_array_equal(actual,reference)
        assert actual.tobytes()==reference.tobytes()


def test_saved_recorder_can_use_strict_original_native_queries():
    from particle_3d.particle81_simulation import SavedTrajectoryStepper
    shapes,provider,wall,_=synthetic('wall');particles=records_for(shapes,provider)
    policy=ValidationNeighborPolicy(20e-6,.5e-6,'P81_NATIVE_RECORDER_TEST')
    recorded=SavedTrajectoryStepper(particles,policy,provider,MU,wall=wall,memoize=False)
    native=Particle65Stepper(particles,policy,provider,MU,wall=wall)
    assert recorded.wall is wall
    for k in range(1,31):
        recorded.step_to(k*.0005);native.step_to(k*.0005)
        np.testing.assert_array_equal(recorded.read()[0].position,native.read()[0].position)
        np.testing.assert_array_equal(recorded.read()[0].velocity,native.read()[0].velocity)
    assert recorded.samples[-1][0]==native.time_s
    assert recorded.cache_info().hits==0


def test_safety_stop_keeps_an_exact_prefix_of_original_long_native_run(root):
    original=np.load(root/'data/reference_original_run/trajectories/mb_000013.npz')['samples']
    truncated=np.load(root/'trajectories/mb_000013.npz')['samples']
    assert original[-1,0]==1.5 and len(original)>len(truncated)
    np.testing.assert_array_equal(original[:len(truncated)],truncated)
    assert original[:len(truncated)].tobytes()==truncated.tobytes()
    metadata=read(root/'trajectories/mb_000013.json')
    assert not metadata['completed'] and metadata['end_reason']=='INTEGRATION_SAFETY_STOP'
    assert metadata['failure_detail'].startswith('EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS')
