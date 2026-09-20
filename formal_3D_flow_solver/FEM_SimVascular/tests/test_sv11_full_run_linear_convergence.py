import pytest
from sv_validation.sv11 import load,REPORT,linear_gate

def test_full_run_every_linear_solve():
    if not (REPORT/'petsc_full_execution.json').exists():pytest.skip('Full run not authorized by short gate')
    assert linear_gate(load('petsc_full_execution'))
    assert load('petsc_full_qc')['completed_steps']>=400
