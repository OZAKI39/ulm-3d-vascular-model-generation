from stage03_helpers import *

def test_primary_p2_p1_and_derived_qc_outputs_remain_distinct():
    require_solution()
    for name in ('solution','qc','metadata','checkpoints','logs'):assert (BASE/name).is_dir()
    primary=np.load(BASE/'checkpoints/primary.npz');visual=np.load(BASE/'solution/fields_for_visualization.npz')
    assert primary['velocity_values'].shape==(preflight()['velocity_dofs']//3,3)
    assert primary['pressure_values'].shape==(preflight()['pressure_dofs'],1)
    assert visual['cells'].shape==(preflight()['mesh']['cells'],10)
    assert visual['velocity_gradient_s_inv'].shape==(preflight()['mesh']['cells'],9)
    assert visual['strain_rate_tensor_s_inv'].shape==(preflight()['mesh']['cells'],9)
    assert visual['vorticity_s_inv'].shape==(preflight()['mesh']['cells'],3)
    assert not qc()['wss_computed'] and 'DG0' in qc()['derived_definitions']['representation']
    s=solve();assert s['config_sha256'] and s['code_sha256'] and s['mesh']['source_sha256']


def test_failed_attempt_never_exports_fake_solution_or_flow_integrals():
    f=read('outputs/stage03/reference/metadata/failure.json');u=read('outputs/stage03/reference/qc/unavailable.json')
    assert f['status']=='FAIL' and not f['pde_solution_available'] and f['roundtrip']=='NOT_EXECUTED_NO_SOLUTION'
    assert u['no_target_flow_substitution'] and u['no_zero_filled_fields']
    assert not (BASE/'checkpoints/primary.npz').exists() and not (BASE/'qc/solution.json').exists()
