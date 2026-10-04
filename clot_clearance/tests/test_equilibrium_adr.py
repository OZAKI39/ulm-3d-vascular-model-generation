"""Independent roots, paper updates and author guards for the additive ADR path."""
import copy,json
from pathlib import Path
import numpy as np
import pytest
from pd_clot.regularization.equilibrium_adr import equilibrate,make_fixture,FrozenForceEvaluator,adaptive_damping,tangent_row_sum
from pd_clot.regularization.equilibrium import ForceState
from pd_clot.regularization.erosion import tensile_strain

ROOT=Path(__file__).resolve().parents[1]


def config():
    c=json.loads((ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json').read_text())
    c['clot']['cells']=[4,3,3];c['damage']['enabled']=False
    return c


def solve(c,method=None):
    cloud=make_fixture(c)
    return cloud,equilibrate(cloud,cloud.X,np.zeros_like(cloud.X),np.ones(len(cloud.pairs)),c,method=method,verification_only=True)


@pytest.fixture(scope='module')
def benchmark():
    c=config();cloud,r=solve(c);_,ref=solve(c,'scipy_reference')
    return c,cloud,r,ref


def test_adr_matches_independent_reference(benchmark):
    c,cloud,r,ref=benchmark;assert r.converged and ref.converged
    free=~cloud.fixed;diff=r.x[free]-ref.x[free];umax=np.linalg.norm(ref.x-cloud.X,axis=1).max()
    assert np.sqrt(np.mean(np.sum(diff**2,axis=1)))/umax<=.01
    assert abs(r.summary['final']['maximum_epsilon_plus']/ref.summary['final']['maximum_epsilon_plus']-1)<=.01


def test_adr_residual_reduction_and_independent_recomputation(benchmark):
    c,cl,r,_=benchmark;s=FrozenForceEvaluator(cl,np.ones(len(cl.pairs)),c)(r.x)
    rel=np.linalg.norm(s.residual[~cl.fixed])/np.linalg.norm(s.external_force_N[~cl.fixed])
    assert rel<=1e-4
    assert rel<r.history[0]['relative_force_residual']*1e-4


def test_fixed_exact_and_no_physical_velocity(benchmark):
    _,cl,r,_=benchmark;np.testing.assert_array_equal(r.x[cl.fixed],cl.X[cl.fixed])
    np.testing.assert_array_equal(r.pseudo_half_velocity[cl.fixed],0.)
    assert r.v is None and all(h['physical_kinetic_energy_J'] is None for h in r.history)


def test_mass_bound_and_three_joint_criteria(benchmark):
    c,_,r,_=benchmark;assert r.summary['endpoint_mass_bound_ratio']<1
    assert all(h['all_criteria_passed'] for h in r.checks[-5:])
    assert r.checks[-1]['pseudo_kinetic_ratio']<=1e-6
    assert r.checks[-1]['driver_drift']<=.01


def test_zero_load_undeformed_residual():
    c=config();c['pipe']['flow_rate_m3_s']=0.;c['streaming']['streaming_velocity_scale_m_s']=0.
    cl=make_fixture(c);s=FrozenForceEvaluator(cl,np.ones(len(cl.pairs)),c)(cl.X)
    assert np.linalg.norm(s.residual[~cl.fixed])<1e-15


def test_extra_iterations_do_not_change_driver(benchmark):
    c=copy.deepcopy(benchmark[0]);c['mechanical_equilibrium']['adr']['additional_iterations_after_convergence']=500
    _,r=solve(c);assert r.converged
    assert r.summary['iterations']>=benchmark[2].summary['iterations']+500
    assert abs(r.summary['final']['maximum_epsilon_plus']/benchmark[2].summary['final']['maximum_epsilon_plus']-1)<=.01


def test_adr_static_resolved_field():
    c=config();c['streaming'].update(provider='resolved_field',field_path=str(ROOT/'verification/mechanical_equilibrium_v1/cases/small_scipy_resolved/static_field.vtu'),
        metadata_path=dict(length_unit='m',velocity_unit='m/s',pressure_unit='Pa',time_unit='s',normal_convention='outward_solid',sample_domain='convex_hull'))
    _,r=solve(c);assert r.converged and r.summary['final']['maximum_epsilon_plus']>0


def test_density_integrity_and_input_are_unchanged():
    c=config();before=copy.deepcopy(c);cl=make_fixture(c);g=np.full(len(cl.pairs),.999);x=cl.X.copy()
    r=equilibrate(cl,x,np.zeros_like(x),g,c)
    assert r.converged and c==before
    np.testing.assert_array_equal(g,.999);np.testing.assert_array_equal(x,cl.X)
    assert r.summary['physical_density_kg_m3']==1056.


def test_legacy_dispatch_preserves_saved_fixture():
    _,r=solve(config(),'legacy_fixed_damping')
    z=np.load(ROOT/'verification/streaming_erosion_v1/tests/erosion_coupled0/baseline/states.npz')
    np.testing.assert_array_equal(r.x,z['x'][1]);np.testing.assert_array_equal(r.v,z['v'][1])


def test_author_negative_quotient_guard_and_cap():
    ones=np.ones((1,3));zeros=np.zeros_like(ones)
    c,q,cap=adaptive_damping(ones,ones,zeros,ones,ones,1.)
    assert q==-1 and c==0 and not cap
    c,q,cap=adaptive_damping(ones,-4*ones,zeros,ones,ones,1.)
    assert q==4 and c==1.9 and cap


def test_simple_linear_load_analytic_solution():
    c=config();cl=make_fixture(c);k=5.;load=np.zeros_like(cl.X);load[~cl.fixed,0]=1e-7
    def force(x):
        u=x-cl.X;eps=tensile_strain(np.linalg.norm(x[cl.pairs[:,1]]-x[cl.pairs[:,0]],axis=1)/cl.length)
        return ForceState(-k*u,load,eps,.5*k*float(np.sum(u*u)),1.,np.zeros_like(x),~cl.fixed)
    r=equilibrate(cl,cl.X,np.zeros_like(cl.X),np.ones(len(cl.pairs)),c,force_callback=force)
    assert r.converged
    np.testing.assert_allclose(r.x-cl.X,load/k,atol=1e-11,rtol=0)
    rows=tangent_row_sum(cl,cl.X,force,cl.spacing*1e-5)
    np.testing.assert_allclose(rows,k,rtol=1e-9)


def test_scipy_stays_restricted_to_small_verification():
    c=config();cl=make_fixture(c)
    with pytest.raises(ValueError,match='verification_only'):
        equilibrate(cl,cl.X,np.zeros_like(cl.X),np.ones(len(cl.pairs)),c,method='scipy_reference')
