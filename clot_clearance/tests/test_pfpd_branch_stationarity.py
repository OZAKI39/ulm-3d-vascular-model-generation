import numpy as np
import pytest
from pd_clot.pfpd_candidate.coupling_stability_audit import load_audit_state,read_json,read_npz,AuditProblem,physics_config
from pd_clot.pfpd_candidate.branch_stability_tracking import (
    OUT,settings,solve_stationary_branch_state,stationary_residual,residual_jacobian,BranchDomainExit)


def test_branch_solver_reproduces_A_live():
    p,d=load_audit_state('A');r=solve_stationary_branch_state(p,d['x'],d['s'])
    assert r['accepted'];np.testing.assert_allclose(r['x'],d['x'],rtol=0,atol=settings()['reproduction_tolerances']['x_absolute_m'])
    np.testing.assert_allclose(r['s'],d['s'],rtol=0,atol=settings()['reproduction_tolerances']['s_absolute'])


def test_branch_solver_reproduces_B_live():
    p,d=load_audit_state('B');r=solve_stationary_branch_state(p,d['x'],d['s'])
    assert r['accepted'];np.testing.assert_allclose(r['x'],d['x'],rtol=0,atol=settings()['reproduction_tolerances']['x_absolute_m'])
    np.testing.assert_allclose(r['s'],d['s'],rtol=0,atol=settings()['reproduction_tolerances']['s_absolute'])


def test_all_saved_stationary_force_and_KKT_gates_pass():
    for f in (OUT/'branch_states').glob('*/STATIONARITY.json'):
        r=read_json(f);assert r['accepted'] and all(r['gates'].values())
        assert r['metrics']['force_residual']<=settings()['force_tolerance']
        assert r['metrics']['projected_phase_KKT']<=settings()['phase_KKT_tolerance']


def test_continuation_lower_bound_and_phase_predictor_are_previous_load_exactly():
    for f in (OUT/'branch_states').glob('*/STATIONARITY.json'):
        r=read_json(f)
        if r['initial_kind']!='continuation_predictor':continue
        a=read_npz(f.with_name('STATE.npz'));parent=read_npz(OUT/'branch_states'/r['previous_accepted_state']/'STATE.npz')
        np.testing.assert_array_equal(a['lower'],parent['s'])
        np.testing.assert_array_equal(a['z_initial'][-len(a['s']):],parent['s'])
        assert r['phase_predictor_is_previous_phase']


def test_tracker_stops_at_model_domain_violation():
    p,d=load_audit_state('A');x=d['x'].copy();x[p.free,0]+=1e-5
    assert p.metrics(x,d['s'])['minimum_Green_strain']<-1e-9
    with pytest.raises(BranchDomainExit,match='PURE_TENSION_DOMAIN_VIOLATION'):stationary_residual(p,p.pack(x,d['s']))
    r=solve_stationary_branch_state(p,x,d['s'])
    assert not r['accepted'] and r['exit']['reason']=='PURE_TENSION_DOMAIN_VIOLATION'


def test_tracker_stops_at_sc_gate():
    p,d=load_audit_state('A');s=np.full(p.ns,p.sc)
    r=solve_stationary_branch_state(p,d['x'],s)
    assert not r['accepted'] and r['exit']['reason']=='HIGH_DAMAGE_GATE_REACHED'


def test_residual_jacobian_chain_rule_matches_live_residual_differences():
    p,d=load_audit_state('B');z=p.pack(d['x'],d['s']);J=residual_jacobian(p,z)
    rng=np.random.default_rng(76543)
    for _ in range(5):
        # Use the largest already-declared derivative scale for the independent
        # residual FD reference. Smaller scales amplify subtraction roundoff;
        # all three scales and the original failed assertion are retained.
        v=rng.normal(size=len(z));v/=np.linalg.norm(v);eps=settings()['hessian_epsilon'][0]
        fd=(stationary_residual(p,z+eps*v,False)-stationary_residual(p,z-eps*v,False))/(2*eps)
        np.testing.assert_allclose(J@v,fd,rtol=2e-5,atol=1e-8)
    sweep=read_json(OUT/'tests/RESIDUAL_JACOBIAN_EPSILON_STUDY.json')
    assert len(sweep)==15 and all(r['relative_L2_error']<2e-5 for r in sweep)


def test_fixed_boundaries_and_inputs_remain_exact():
    p,d=load_audit_state('B');before_x=d['x'].copy();before_s=d['s'].copy()
    r=solve_stationary_branch_state(p,d['x'],d['s'])
    np.testing.assert_array_equal(d['x'],before_x);np.testing.assert_array_equal(d['s'],before_s)
    np.testing.assert_array_equal(r['x'][p.cloud.fixed],p.boundary[p.cloud.fixed])
