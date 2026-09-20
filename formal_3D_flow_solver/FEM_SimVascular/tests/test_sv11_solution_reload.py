import pytest
from sv_validation.sv11 import load,REPORT

def test_fresh_process_reload():
    if not (REPORT/'solution_reload.json').exists():pytest.skip('No vascular solution to reload')
    assert load('solution_reload')['status']=='PASS'
    assert load('solution_reload')['fresh_process'] is True
