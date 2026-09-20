import copy
import pytest
from sv_validation.sv11 import parse_solver_log,load,nonlinear_gate
from sv_validation.validation import ValidationError

def test_nonlinear_failure_recognized():
    text='[svMultiPhysics] WARNING: The number of nonlinear iterations (13) has exceeded the maximum number set by the value of the Add_equation/Max_iterations parameter'
    assert parse_solver_log(text,1)['nonlinear_failure_messages']

def test_smoke_nonlinear_residual_decreases():
    rows=load('petsc_smoke')['run']['history']['linear_solves']
    assert rows[-1]['nonlinear_Ri_over_R0']<rows[0]['nonlinear_Ri_over_R0']

def test_silent_iteration_cap_does_not_mean_convergence():
    history=copy.deepcopy(load('petsc_smoke')['run']['history'])
    row=history['linear_solves'][-1];row.update(nonlinear_iteration=12,nonlinear_Ri_over_R0=.1,nonlinear_Ri_over_R1=.1)
    with pytest.raises(ValidationError,match='cap'):nonlinear_gate(history)
