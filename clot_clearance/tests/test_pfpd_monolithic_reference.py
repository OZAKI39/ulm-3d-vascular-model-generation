import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import OUT,audit_config,read_json,read_npz


def test_independent_optimizers_and_all_predeclared_starts_retained():
    a=audit_config(); rows=read_json(OUT/'monolithic/ALL_OPTIMIZERS.json')
    failure=[r for r in rows if r['strain']==a['subdivision_endpoint']]
    assert {(r['method'],r['start']) for r in failure}=={(m,s) for m in a['optimizer_methods'] for s in a['monolithic_starts']}


def test_first_order_candidate_meets_bounds_and_original_force_and_domain_gates():
    a=audit_config(); rows=read_json(OUT/'monolithic/ALL_OPTIMIZERS.json'); passed=0
    for r in rows:
        if not r['first_order_candidate_passed']:continue
        passed+=1; b=read_npz(OUT/'monolithic'/r['method']/r['start']/'STATE.npz'); met=r['metrics']
        assert r['optimizer_success'] and r['exit'] is None
        assert np.all(b['s']>=b['lower']-a['optimizer_bounds_tolerance']) and np.all(b['s']<.95)
        assert met['force_residual']<=1e-9 and met['minimum_Green_strain']>=-1e-9
        assert met['projected_phase_KKT']<=a['phase_KKT_tolerance']
    assert passed>0


def test_saddle_points_and_domain_exit_trials_are_not_valid_minima():
    rows=read_json(OUT/'monolithic/ACCEPTANCE_AND_COMPARISON.json')['rows']
    for r in rows:
        if r['domain_exit']:
            assert not r['valid_local_minimum'] and r['recorded_point_kind']=='last_domain_valid_evaluation_not_a_solution'
        if r['curvature']=='NEGATIVE_CURVATURE_FOUND':assert not r['valid_local_minimum']


def test_unreached_subdivision_endpoint_has_no_invented_metrics():
    rows=read_json(OUT/'continuation/ENDPOINT_COMPARISON.json')['rows']
    assert {r['substeps'] for r in rows}=={2,4,8}
    for r in rows:
        if not r['endpoint_reached']:
            assert all(v is None for k,v in r.items() if k.startswith('endpoint_') and k!='endpoint_reached')
