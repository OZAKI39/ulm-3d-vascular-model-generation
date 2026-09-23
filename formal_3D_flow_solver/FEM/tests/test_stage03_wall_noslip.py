from stage03_helpers import *

def test_actual_wall_coefficients_are_within_scaled_roundoff():
    q=qc();assert q['max_abs_wall_velocity_m_s']<=q['wall_roundoff_tolerance_m_s']
    expected=128*np.finfo(float).eps*config()['physics']['inlet_volume_flow_m3_s']/config()['inlet_geometry']['projected_area_m2']
    close(q['wall_roundoff_tolerance_m_s'],expected)
