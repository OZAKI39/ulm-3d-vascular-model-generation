import numpy as np


def test_three_validation_dts_match_exact_straight_line(uniform_rows):
    dts=sorted({r['dt_s'] for r in uniform_rows})
    np.testing.assert_array_equal(np.array(dts[1:])/np.array(dts[:-1]),[2.,2.])
    for row in uniform_rows:
        assert row['timestep_role']=='VALIDATION_ONLY'
        assert row['position_error_m']<=row['position_error_bound_m']
        assert row['velocity_error_m_s']==0
    print('uniform max position error (m)',max(r['position_error_m'] for r in uniform_rows))
