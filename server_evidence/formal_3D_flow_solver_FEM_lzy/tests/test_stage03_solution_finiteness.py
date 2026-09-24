from stage03_helpers import *

def test_original_coefficients_and_derived_samples_are_all_finite():
    require_solution()
    for path in ('checkpoints/primary.npz','solution/fields_for_visualization.npz'):
        with np.load(BASE/path) as data:
            assert all(np.isfinite(data[k]).all() for k in data.files)
    assert qc()['all_fields_finite'] and solve()['solver']['physical_fields_finite']
    assert solve()['solver']['converged_reason']>0 and np.isfinite(qc()['lambda_pa'])
    assert not qc()['divergence_hard_gate']
