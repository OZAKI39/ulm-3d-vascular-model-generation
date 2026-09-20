import pytest
from sv_validation.sv11 import load,REPORT,linear_gate

def test_actual_short_linear_convergence():
    if not (REPORT/'petsc_short_execution.json').exists():pytest.skip('Short run blocked by preceding gate')
    assert linear_gate(load('petsc_short_execution'))
    assert load('petsc_short_qc')['completed_steps']==10
