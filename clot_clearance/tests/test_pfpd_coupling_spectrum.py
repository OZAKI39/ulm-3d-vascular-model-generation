import numpy as np
from scipy.linalg import eigvals
from pd_clot.pfpd_candidate.coupling_stability_audit import OUT,audit_config,read_json,read_npz


def test_stable_control_epsilon_stability_without_assuming_contraction():
    a=audit_config(); r=read_json(OUT/'spectrum/EPSILON_STABILITY.json')['A']
    assert r['ordinary_status']=='RESOLVED'
    values=np.array([v['selected_branch_rho'] for v in r['rows']])
    assert np.ptp(values)/max(abs(values))<=a['spectral_epsilon_relative_tolerance']
    # Classification follows the measured result, not an imposed rho < 1.
    rho=r['ordinary_rho']
    assert ('contractive' if rho<1-a['rho_margin'] else 'unstable' if rho>1+a['rho_margin'] else 'marginal')=='marginal'


def test_eigenpairs_reproduce_stored_matrices():
    for key in ['A','B','C']:
        b=read_npz(OUT/f'spectrum/{key}/eps_3e-06/MATRICES.npz')
        J=b['frozen_selected_branch_J']; vals=b['eigenvalues']; vecs=b['eigenvectors']
        np.testing.assert_allclose(J@vecs,vecs*vals,atol=1e-10,rtol=1e-9)
        np.testing.assert_allclose(max(abs(eigvals(J))),max(abs(vals)),rtol=1e-10)


def test_switched_columns_never_reported_as_ordinary_jacobian():
    summaries=read_json(OUT/'spectrum/EPSILON_STABILITY.json')
    for key in ['B','C']:
        assert summaries[key]['ordinary_rho'] is None
        for r in summaries[key]['rows']:
            assert r['switched_columns']>0
            assert r['conditional_frozen_branch_only'] and not r['ordinary_Phi_spectrum_valid']


def test_frozen_trial_derivative_matches_independent_energy_blocks():
    for r in read_json(OUT/'comparisons/HESSIAN_MAP_CROSSCHECK.json'):
        assert r['trial_derivative_Hessian_vs_FD_relative_error']<audit_config()['spectral_epsilon_relative_tolerance']


def test_nonsmooth_directional_evidence_preserved():
    r=read_json(OUT/'spectrum/MODE_ALIGNMENT.json')['directional_probes']
    assert len(r)==6
    assert all(x['selector_switches']>0 for x in r)
    # Both sides are retained: agreement along one side cannot erase the other.
    assert max(x['relative_JVP_error'] for x in r if x['sign']<0)>0.5
    assert max(x['relative_JVP_error'] for x in r if x['sign']>0)<1e-4
