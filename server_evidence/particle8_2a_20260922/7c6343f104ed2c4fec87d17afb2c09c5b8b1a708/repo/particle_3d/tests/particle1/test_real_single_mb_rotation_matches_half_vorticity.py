import numpy as np


def test_all_saved_angular_velocities_match_fresh_field_samples(trajectory_cases,trajectory_samples):
    for (rows,summary),samples in zip(trajectory_cases,trajectory_samples):
        omega=np.array([[float(r[f'Omega_{a}_s_inv']) for a in 'xyz'] for r in rows])
        vorticity=np.array([[float(r[f'vorticity_{a}_s_inv']) for a in 'xyz'] for r in rows])
        np.testing.assert_array_equal(vorticity,samples.vorticity_s_inv)
        np.testing.assert_array_equal(omega,.5*samples.vorticity_s_inv)
        assert summary['max_angular_relation_error_s_inv']==0
