from particle_3d.particle1_cases import SYNTHETIC_ATOL,SYNTHETIC_RTOL


def test_rigid_rotation_recovers_analytic_omega_and_zero_strain(rotation_rows):
    assert len(rotation_rows)==49
    for row in rotation_rows:
        assert row['strain_norm_s_inv']==0.
        for axis in 'xyz':
            expected=row[f'analytic_omega_{axis}_s_inv']
            assert abs(row[f'sampled_omega_{axis}_s_inv']-expected)<=SYNTHETIC_ATOL+SYNTHETIC_RTOL*abs(expected)
    print('rotation max error (1/s)',max(r['angular_velocity_error_s_inv'] for r in rotation_rows))
