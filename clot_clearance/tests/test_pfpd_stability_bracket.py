import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import read_json
from pd_clot.pfpd_candidate.branch_stability_tracking import OUT,settings


def test_stable_endpoint_classification_is_reproducible_from_curvature_gates():
    for k in ['A','B']:
        r=read_json(OUT/f'critical_bracket/SCHEDULE_{k}.json')['lower']
        assert r['stationary'] and r['classification']=='STABLE_SAMPLED_STATE'
        assert all(v>t for v,t in zip(r['reduced_eigenvalues_min'],r['curvature_thresholds']))
        assert r['active']['counts']['LOWER_WEAK_NEAR_SWITCH']==0


def test_unstable_endpoint_has_resolved_critical_cone_negative_curvature():
    for k in ['A','B']:
        r=read_json(OUT/f'critical_bracket/SCHEDULE_{k}.json')['upper'];w=r['critical_cone_negative_witness']
        assert r['stationary'] and r['classification']=='UNSTABLE_SAMPLED_STATE'
        assert w['critical_cone'] and w['feasible'] and all(v<-t for v,t in zip(w['curvatures'],r['curvature_thresholds']))


def test_refinement_shrinks_bracket_and_uses_only_lower_forward_predictors():
    for k in ['A','B']:
        b=read_json(OUT/f'critical_bracket/SCHEDULE_{k}.json')
        assert b['width']<=settings()['refinement_width']
        for r in b['refinement']:
            assert r['state']['strain']==sum(r['before'])/2
            assert np.isclose(r['after'][1]-r['after'][0],.5*(r['before'][1]-r['before'][0]),rtol=1e-12)
            stat=read_json(OUT/'branch_states'/r['state']['state']/'STATIONARITY.json')
            parent=read_json(OUT/'branch_states'/stat['previous_accepted_state']/'STATIONARITY.json')
            assert parent['strain']==r['before'][0]


def test_schedule_consistency_matches_predeclared_tolerances():
    r=read_json(OUT/'comparisons/SCHEDULE_CONSISTENCY.json');c=settings()['schedule_consistency']
    assert r['status']=='PASSED' and r['endpoint_distance']<=c['critical_interval_tolerance']
    assert all(x['max_phase_difference']<=c['same_load_max_s_tolerance'] and x['active_set_identical'] for x in r['matching_loads'])


def test_critical_mode_energy_curvature_independently_matches_largest_valid_probes():
    r=read_json(OUT/'energy_landscape/nearest_interior_unstable.json');points={p['alpha']:p for p in r['points']}
    alpha=max(a for a in points if a>0 and points[a]['feasible'] and points[-a]['feasible'])
    second=(points[alpha]['delta_energy_J']+points[-alpha]['delta_energy_J'])/(alpha**2*r['energy_scale_J'])
    assert second<0 and r['curvature']<0
    np.testing.assert_allclose(second,r['curvature'],rtol=settings()['hessian_epsilon_relative_tolerance'])
    stable=read_json(OUT/'energy_landscape/last_stable.json')
    assert all(p['delta_energy_J']>0 for p in stable['points'] if p['feasible'] and p['alpha']!=0)


def test_invalid_landscape_directions_are_not_published_as_valid_energies():
    for label in ['last_stable','onset_cone_direction','nearest_interior_unstable']:
        r=read_json(OUT/f'energy_landscape/{label}.json')
        assert set(settings()['landscape_alpha'])<={p['alpha'] for p in r['points']}
        for p in r['points']:
            if not p['feasible']:assert p['delta_energy_J'] is None
