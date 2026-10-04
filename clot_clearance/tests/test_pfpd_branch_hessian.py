import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import read_json,read_npz,AuditProblem,physics_config
from pd_clot.pfpd_candidate.branch_stability_tracking import OUT,settings
from scripts.run_monolithic_stability_audit import full_hessian


def test_new_audit_reuses_identical_hessian_assembly():
    b=read_npz(OUT/'hessian/A_grid_00/eps_3e-06/BLOCKS.npz')
    p=AuditProblem(physics_config(),.009,b['lower']);H,_=full_hessian(p,p.pack(b['x'],b['s']),3e-6)
    np.testing.assert_array_equal(H,b['H_raw'])


def test_all_hessian_scales_have_consistent_signs_and_original_block_gates():
    for f in (OUT/'hessian').glob('*/BRANCH_SUMMARY.json'):
        r=read_json(f);vals=np.array(r['reduced_eigenvalues_min']);t=np.array(r['curvature_thresholds'])
        sign=np.where(vals>t,1,np.where(vals<-t,-1,0))
        assert r['curvature_sign_resolved'] and np.all(sign==sign[0])
        for eps in settings()['hessian_epsilon']:
            a=read_json(f.parent/f'eps_{eps:g}/SUMMARY.json')
            assert a['Hessian_symmetry_relative_error']<=1e-5 and a['Hss_analytic_relative_error']<=1e-6


def test_strong_bound_components_are_absent_from_critical_negative_witnesses():
    for f in (OUT/'hessian').glob('*/BRANCH_SUMMARY.json'):
        r=read_json(f);data=read_npz(f.parent/'CONE_AND_MODE.npz');nx=len(data['critical_mode'])-len(data['s'])
        strong=np.array(r['active']['bond_IDs']['LOWER_STRONG_ACTIVE'],dtype=int)
        for probe in read_json(f.parent/'CONE_PROBES.json'):
            if probe['critical_cone']:
                assert abs(probe['linear_energy_variation'])<=settings()['critical_linear_tolerance']
                assert np.max(abs(data['directions'][probe['index'],nx+strong]),initial=0)<1e-14


def test_weak_bound_negative_curvature_is_a_feasible_witness_not_released_eigenvalue():
    bracket=read_json(OUT/'critical_bracket/SCHEDULE_B.json');r=bracket['upper'];w=r['critical_cone_negative_witness']
    assert r['active']['counts']['LOWER_WEAK_NEAR_SWITCH']>0
    assert r['indicator_kind']=='feasible_critical_cone_witness'
    assert w['critical_cone'] and w['feasible']
    assert all(v<-t for v,t in zip(w['curvatures'],r['curvature_thresholds']))
    assert all(v>0 for v in r['reduced_eigenvalues_min'])
