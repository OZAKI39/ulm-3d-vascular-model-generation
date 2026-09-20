import pytest
from sv_validation.sv11 import load,REPORT

def test_accepted_wall_noslip():
    if not (REPORT/'accepted_solution.json').exists():pytest.skip('No accepted steady solution')
    q=load('accepted_solution');assert q['wall_noslip_pass']
    assert q['wall_velocity_max_m_s']<=q['wall_noslip_tolerance_m_s']
