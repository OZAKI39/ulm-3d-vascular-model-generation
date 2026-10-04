"""Independent force-root benchmark and preserved legacy integrator checks."""
import copy,json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.regularization.equilibrium import make_fixture,FrozenForceEvaluator,equilibrate,diagnostics,NeedPaperInput

ROOT=Path(__file__).resolve().parents[1]


def small_config():
    c=json.loads((ROOT/'configs/mechanical_equilibrium_verification_v1.json').read_text())
    c['clot']['cells']=c['mechanical_equilibrium']['small_fixture_cells']
    return c


def solve(c,method='scipy_reference'):
    cloud=make_fixture(c);g=np.ones(len(cloud.pairs))
    return cloud,equilibrate(cloud,cloud.X,np.zeros_like(cloud.X),g,c,method=method,verification_only=True)


@pytest.fixture(scope='module')
def benchmark():
    c=small_config()
    cloud,reference=solve(c)
    _,legacy=solve(c,'legacy_fixed_damping')
    return c,cloud,reference,legacy


def test_zero_load_undeformed_force_residual_is_small():
    c=small_config();c['pipe']['flow_rate_m3_s']=0.;c['streaming']['streaming_velocity_scale_m_s']=0.
    cloud=make_fixture(c);state=FrozenForceEvaluator(cloud,np.ones(len(cloud.pairs)),c)(cloud.X)
    row=diagnostics(cloud,cloud.X,state,c,np.zeros_like(cloud.X))
    assert row['force_residual_L2_N']<1e-15
    assert row['relative_force_residual']<1e-4
    assert row['physical_kinetic_energy_J']==0


def test_scipy_reference_independently_balances_free_force(benchmark):
    c,cloud,result,_=benchmark
    assert result.converged
    state=FrozenForceEvaluator(cloud,np.ones(len(cloud.pairs)),c)(result.x)
    ratio=np.linalg.norm(state.residual[~cloud.fixed])/max(np.linalg.norm(state.external_force_N[~cloud.fixed]),1e-12)
    assert ratio<=1e-4
    assert result.summary['scipy_success']
    assert result.summary['scipy_termination_reason']
    assert result.summary['free_coordinate_count']==3*int((~cloud.fixed).sum())


@pytest.mark.parametrize('which',[2,3])
def test_fixed_particles_remain_exact_in_both_methods(benchmark,which):
    _,cloud,*_=benchmark;result=benchmark[which]
    np.testing.assert_array_equal(result.x[cloud.fixed],cloud.X[cloud.fixed])
    if result.v is not None:np.testing.assert_array_equal(result.v[cloud.fixed],0.)


def test_legacy_reproduces_saved_erosion_v1_pre_damage_coordinates(benchmark):
    legacy=benchmark[3]
    with np.load(ROOT/'verification/streaming_erosion_v1/tests/erosion_coupled0/baseline/states.npz') as saved:
        np.testing.assert_array_equal(legacy.x,saved['x'][1])
        np.testing.assert_array_equal(legacy.v,saved['v'][1])


def test_scipy_uses_no_physical_velocity_or_pseudo_time(benchmark):
    result=benchmark[2]
    assert result.v is None
    assert result.summary['mechanical_time_s'] is None
    assert result.summary['equilibrium_pseudo_time'] is None
    assert all(row['physical_kinetic_energy_J'] is None for row in result.history)


@pytest.mark.parametrize('which',[2,3])
def test_no_damage_and_no_erosion_clock_in_equilibrium_benchmarks(benchmark,which):
    result=benchmark[which]
    assert not result.summary['damage_advanced']
    for row in result.history:
        assert row['erosion_exposure_time_s']==0
        assert row['maximum_damage']==0


def test_legacy_failure_is_explicit_and_does_not_advance_damage():
    c=small_config();c['mechanical_equilibrium']['legacy_window_s']=2e-6
    _,result=solve(c,'legacy_fixed_damping')
    assert not result.converged
    assert result.summary['status']=='INCOMPLETE_MECHANICAL_EQUILIBRIUM'
    assert not result.summary['damage_advanced']


def test_adr_request_requires_primary_paper_and_preserves_inputs():
    c=small_config();cloud=make_fixture(c);g=np.ones(len(cloud.pairs));x=cloud.X.copy();v=np.zeros_like(x)
    with pytest.raises(NeedPaperInput,match='NEED_PAPER_INPUT'):
        equilibrate(cloud,x,v,g,c,method='adaptive_dynamic_relaxation')
    np.testing.assert_array_equal(x,cloud.X);np.testing.assert_array_equal(g,1.);np.testing.assert_array_equal(v,0.)


def test_scipy_cannot_be_used_as_automatic_production_solver():
    c=small_config();cloud=make_fixture(c)
    with pytest.raises(ValueError,match='verification_only'):
        equilibrate(cloud,cloud.X,np.zeros_like(cloud.X),np.ones(len(cloud.pairs)),c,method='scipy_reference')
    c['clot']['cells']=[10,6,6]
    with pytest.raises(ValueError,match='small cases'):solve(c)


def test_force_callback_recomputes_position_dependent_surface_load():
    c=small_config();cloud=make_fixture(c);evaluator=FrozenForceEvaluator(cloud,np.ones(len(cloud.pairs)),c)
    x=cloud.X.copy();before=evaluator(x);x[~cloud.fixed,0]+=1e-8;after=evaluator(x)
    assert np.linalg.norm(before.external_force_N-after.external_force_N)>1e-20


def test_deltaN_and_frequency_metadata_do_not_change_reference(benchmark):
    c=copy.deepcopy(benchmark[0]);c['damage']['DeltaN']=9123456
    c['simulation']['representative_frequency_Hz']=13456.
    c['streaming']['ultrasound_frequency_Hz']=2.5e9
    _,result=solve(c)
    np.testing.assert_array_equal(result.x,benchmark[2].x)


def test_static_resolved_field_reference_has_nonzero_tensile_driver():
    c=small_config()
    c['streaming'].update(provider='resolved_field',field_path=str(ROOT/'verification/streaming_erosion_v1/tests/test_static_resolved_field_use0/static_field.vtu'),metadata_path=dict(length_unit='m',velocity_unit='m/s',pressure_unit='Pa',time_unit='s',normal_convention='outward_solid',sample_domain='convex_hull'))
    _,result=solve(c)
    assert result.converged
    assert result.summary['final']['maximum_epsilon_plus']>0
    assert result.summary['final']['relative_force_residual']<=1e-4


def test_frozen_integrity_and_physical_density_are_unchanged():
    c=small_config();before=copy.deepcopy(c);cloud=make_fixture(c);g=np.ones(len(cloud.pairs))
    equilibrate(cloud,cloud.X,np.zeros_like(cloud.X),g,c,method='scipy_reference',verification_only=True)
    assert c==before
    np.testing.assert_array_equal(g,1.)
