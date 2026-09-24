import numpy as np


def test_all_saved_positions_requery_original_p0_velocity(trajectory_cases,trajectory_samples):
    for (rows,summary),samples in zip(trajectory_cases,trajectory_samples):
        velocity=np.array([[float(r[f'V_{a}_m_s']) for a in 'xyz'] for r in rows])
        background=np.array([[float(r[f'u_{a}_m_s']) for a in 'xyz'] for r in rows])
        assert samples.inside_lumen.all()
        np.testing.assert_array_equal(velocity,samples.velocity_m_s)
        np.testing.assert_array_equal(background,samples.velocity_m_s)
        np.testing.assert_array_equal([int(r['tetra_id']) for r in rows],samples.tetra_id)
        np.testing.assert_array_equal([float(r['pressure_pa']) for r in rows],samples.pressure_pa)
        assert summary['max_velocity_relation_error_m_s']==0
