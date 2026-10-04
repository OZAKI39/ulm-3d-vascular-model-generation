"""Targeted gates for the isolated C5 activation diagnostic."""
import hashlib
import numpy as np
from scripts.c5_stationarity_diagnostics import (OUT,ROOT,read,ca,Problem,classify,q_fields,
    ReducedSystem,full_acceptance,update_active,selections)


def activation_fixture():
    r=read(OUT/'intact/root_offset_1e-6/SUMMARY.json')
    a=ca.read_npz(OUT/'intact/root_offset_1e-6/STATE.npz');p=Problem(5,r['strain'])
    ids=np.flatnonzero(a['q']<0)
    return p,a,ids


def test_old_accepted_and_three_failed_states_reproduce():
    r=read(OUT/'reproduction/ALL.json')
    assert r['passed'] and len(r['states'])==4


def test_s0_q_is_twice_Yc_minus_Psi_plus():
    p,a,_=activation_fixture()
    np.testing.assert_allclose(a['q'],2*(p.yc-a['Y']),rtol=0,atol=1e-16)


def test_strong_weak_violated_and_interior_classification():
    s=np.array([0.,0.,0.,.01,0.]);qs=np.array([2e-7,5e-8,-2e-7,0.,-5e-8])
    masks=classify(s,np.zeros(5),qs,.95)
    assert np.flatnonzero(masks['STRONG_LOWER']).tolist()==[0]
    assert np.flatnonzero(masks['WEAK_LOWER']).tolist()==[1,4]
    assert np.flatnonzero(masks['VIOLATED_LOWER']).tolist()==[2]
    assert np.flatnonzero(masks['INTERIOR']).tolist()==[3]


def test_activation_bracket_contracts_and_preserves_signs():
    b=read(OUT/'activation/BRACKET.json')
    assert b['width']<=1e-7 and len(b['history'])<=12
    assert b['lower']['min_q']>=0 and b['upper']['min_q']<0
    for r in b['history']:
        assert r['after'][0]>=r['before'][0] and r['after'][1]<=r['before'][1]
        assert r['after'][1]-r['after'][0]<r['before'][1]-r['before'][0]


def test_activation_is_two_local_bonds_not_diffuse():
    r=read(OUT/'activation/SPATIAL.json')[1]
    assert r['candidate_count']==2 and r['fraction']<.001
    assert {b['ID'] for b in r['candidates']}=={4488,5059}


def test_projected_residual_piecewise_directional_derivative():
    # Independent scalar perturbations of s and q_scaled, away from kinks.
    for s,q,ds,dq in [(.2,.4,.3,.7),(.2,.1,.3,.7),(.8,-.3,.3,.7)]:
        e=1e-7;arg=s-q
        exact=ds if arg<0 or arg>.95 else dq
        f=lambda a,b:ca.projected_kkt(np.array([a]),np.array([b]),np.array([0.]),.95)[0]
        fd=(f(s+e*ds,q+e*dq)-f(s-e*ds,q-e*dq))/(2*e)
        assert abs(fd-exact)<1e-9


def test_reduced_residual_uses_original_force_and_active_q():
    p,a,ids=activation_fixture();system=ReducedSystem(p,ids);u=system.compress(p.pack(a['x'],a['s']))
    r=system.residual(u);x,s=p.unpack(system.expand(u));state=p.state(x,s)
    np.testing.assert_allclose(r[:p.nx],state.force[p.free].ravel()/p.config['mechanical_equilibrium']['force_floor_N'])
    np.testing.assert_allclose(r[p.nx:],a['q'][ids]/p.yc,rtol=2e-12,atol=1e-14)


def test_inactive_phase_is_exactly_previous_lower_bound():
    p,a,ids=activation_fixture();system=ReducedSystem(p,ids);u=system.compress(p.pack(a['x'],a['s']));u[p.nx:]=.01
    s=system.expand(u)[p.nx:];mask=np.ones(p.ns,bool);mask[ids]=False
    assert np.array_equal(s[mask],p.lower[mask]) and np.all(s[ids]==.01)


def test_full_KKT_rejects_omitted_active_bonds():
    p,a,ids=activation_fixture();r=full_acceptance(p,p.pack(a['x'],a['s']),np.array([],int))
    assert r['gates']['force'] and not r['valid'] and not r['gates']['inactive_inequalities']


def test_active_set_updates_are_deterministic_and_add_missing_bonds():
    p,a,ids=activation_fixture();z=p.pack(a['x'],a['s'])
    first=update_active(p,z,np.array([],int));second=update_active(p,z,np.array([],int))
    for x,y in zip(first,second):assert np.array_equal(x,y)
    assert np.array_equal(first[0],ids) and np.array_equal(first[1],ids)


def test_A3_selection_preserves_both_initial_sets():
    p,a,_=activation_fixture();sets=selections(p,a['x'],a['s'])
    assert set(sets['A1'])<=set(sets['A3']) and set(sets['A2'])<=set(sets['A3'])


def test_reduced_jacobian_three_step_FD_consistency():
    r=read(OUT/'reduced/first_A1/JACOBIAN_CHECK.json')
    assert r['passed'] and {v['epsilon'] for v in r['rows']}=={1e-5,3e-6,1e-6}
    assert max(v['relative_error'] for v in r['rows'])<5e-5


def test_full_residual_blocks_match_existing_operator():
    for name in ['A_last_accepted','B_first_weak','C_failed_coupled']:
        assert read(OUT/f'jacobian/{name}/VALIDATION.json')['passed']


def test_singular_estimates_require_original_operator_residuals():
    for r in read(OUT/'jacobian/ALL.json'):
        if r['condition_resolved']:
            assert r['sigma_max_relative_residual']<=1e-6 and r['sigma_min_relative_residual']<=1e-6
            assert r['condition_estimate']>=1
        else:assert r['condition_estimate'] is None


def test_zero_step_termination_is_not_physical_convergence():
    r=read(OUT/'original_solver/SUMMARY.json');last=r['termination_checks'][-1]
    assert not r['accepted'] and last['status']==4
    assert last['step_norm']==0 and last['dF']==0
    assert r['metrics']['force_residual']>1e-9 and r['metrics']['projected_phase_KKT']>1e-8


def test_preexisting_source_tests_and_reports_are_unchanged():
    pre=read(OUT/'provenance/PRE_RUN_MANIFEST.json')
    for name,digest in pre['files'].items():
        source=name.endswith('.py') and name.split('/')[0] in {'pd_clot','scripts','tests'}
        report=name.startswith(('verification/','results/')) and 'REVIEW' in name and name.endswith('.md')
        if source or report:assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name


def test_positive_root_passes_every_full_KKT_gate():
    rows=read(OUT/'reduced/INITIAL_RESULTS.json')
    assert len(rows)==3
    for r in rows:
        assert r['valid'] and r['positive_phase'] and all(r['gates'].values())
        assert r['inactive_scaled_violation_max']==0


def test_five_step_prefix_preserves_irreversibility():
    rows=read(OUT/'reduced/CONTINUATION.json')
    assert len(rows)==6 and all(r['valid'] for r in rows)
    previous=None
    for r in rows:
        a=ca.read_npz(OUT/f'reduced/{r["label"]}/STATE.npz')
        assert np.all(a['s']>=a['lower'])
        if previous is not None:
            assert np.array_equal(a['lower'],previous) and np.all(a['s']>=previous)
        previous=a['s']


def test_same_load_resolves_root_without_changing_load_step():
    old=read(OUT/'original_solver/SUMMARY.json');new=read(OUT/'reduced/SAME_LOAD.json')
    assert new['strain']==old['strain'] and new['valid'] and new['positive_phase']
    assert new['metrics']['force_residual']<=1e-9 and new['metrics']['projected_phase_KKT']<=1e-8
