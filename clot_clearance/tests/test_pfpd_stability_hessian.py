import numpy as np
from scipy.linalg import solve,eigh
from pd_clot.pfpd_candidate.coupling_stability_audit import OUT,audit_config,read_json,read_npz


def test_hessian_symmetry_measured_before_symmetrization():
    a=audit_config()
    for key in ['A','B','C']:
        for eps in a['hessian_epsilon']:
            H=read_npz(OUT/f'hessian/{key}/eps_{eps:g}/BLOCKS.npz')['H_raw']
            assert np.linalg.norm(H-H.T)/np.linalg.norm(H)<=a['hessian_symmetry_relative_tolerance']


def test_phase_hessian_matches_positive_diagonal_analytic_energy():
    for key in ['A','B','C']:
        b=read_npz(OUT/f'hessian/{key}/eps_3e-06/BLOCKS.npz'); expected=b['H_ss_analytic']
        assert np.all(np.diag(expected)>0)
        assert np.linalg.norm(b['H_ss']-expected)/np.linalg.norm(expected)<=audit_config()['hss_relative_tolerance']


def test_reduced_hessian_excludes_strong_active_bounds_and_has_valid_eigenpair():
    for key in ['A','B','C']:
        b=read_npz(OUT/f'hessian/{key}/eps_3e-06/BLOCKS.npz'); nx=b['H_xx'].shape[0]
        kept=np.r_[np.arange(nx),nx+b['interior_phase_indices']]
        assert b['reduced'].shape==(len(kept),len(kept))
        v=b['minimum_mode_scaled'][kept]
        np.testing.assert_allclose(b['reduced']@v,b['reduced_eigenvalues'][0]*v,atol=1e-10)
        assert len(kept)==(nx if key=='A' else nx+len(b['s']))


def test_schur_uses_well_conditioned_mechanical_block_and_matches_saved_result():
    for key in ['B','C']:
        b=read_npz(OUT/f'hessian/{key}/eps_3e-06/BLOCKS.npz'); H=b['H_symmetric']; nx=b['H_xx'].shape[0]
        xx=H[:nx,:nx]; xs=H[:nx,nx:]; ss=H[nx:,nx:]
        assert np.linalg.cond(xx)<audit_config()['hxx_condition_limit'] and eigh(xx,eigvals_only=True)[0]>0
        np.testing.assert_allclose(b['Schur'],ss-xs.T@solve(xx,xs,assume_a='pos'),atol=1e-11)


def test_negative_mode_classification_requires_epsilon_stability_and_feasibility():
    summary=read_json(OUT/'hessian/EPSILON_STABILITY.json')
    for key,row in summary.items():
        if row['status']!='NEGATIVE_CURVATURE_FOUND':continue
        assert row['epsilon_stable'] and row['feasible_negative_mode']
        for eps in audit_config()['hessian_epsilon']:
            r=read_json(OUT/f'hessian/{key}/eps_{eps:g}/SUMMARY.json')
            assert r['minimum_reduced_eigenvalue']<-r['curvature_threshold']
            assert any(all(p.get('bounds') and p.get('in_tensile_domain') for p in probe['points']) for probe in r['minimum_mode_probes'])
