from sv12_support import accepted
def test_final_actual_wall_velocity():
    q=accepted();assert q['wall_noslip_pass']
    assert q['wall_velocity_max_m_s']<=q['wall_noslip_tolerance_m_s']
    assert q['wall_velocity_P95_m_s']<=q['wall_noslip_tolerance_m_s']
