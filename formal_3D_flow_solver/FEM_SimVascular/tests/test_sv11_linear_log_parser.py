import copy
import pytest
from sv_validation.sv11 import load,parse_solver_log,linear_gate
from sv_validation.validation import ValidationError

def test_zero_exit_does_not_override_divergence():
    data=copy.deepcopy(load('petsc_smoke')['run']);data['history']['petsc_reasons'][0]['diverged']=True
    assert data['exit_code']==0
    with pytest.raises(ValidationError,match='divergence'):linear_gate(data)

def test_nonfinite_summary_residual_rejected():
    line=' NS      1-1   1.000e+00 [0  1.000e+00  1.000e+00  nan] [   21 -233   84]'
    history=parse_solver_log(line,1.)
    assert len(history['linear_solves'])==1 and history['failed_linear_solves']==1
    assert history['linear_solves'][0]['linear_residual_relative'] is None
    with pytest.raises(ValidationError):linear_gate({'exit_code':0,'history':history},False)

def test_full_precision_monitors_associated_with_each_solve():
    h=load('petsc_smoke')['run']['history']
    assert len(h['petsc_reasons'])==len(h['linear_solves'])
    for row in h['linear_solves']:
        assert row['petsc_reason']['iterations']==row['linear_iterations']
        assert row['petsc_monitor'][0]['iteration']==0
        assert row['petsc_monitor'][-1]['iteration']==row['linear_iterations']
        assert row['linear_residual_absolute']==row['petsc_monitor'][-1]['residual_norm']

def test_nonfinite_true_residual_rejects_finite_summary():
    text=' 1 KSP unpreconditioned resid norm 1e-12 true resid norm nan ||r(i)||/||b|| nan\nLinear solve converged due to CONVERGED_RTOL iterations 1\n NS 1-1 1.0 [0 1.0 1.0 1e-12] [1 -200 99]'
    h=parse_solver_log(text,1.)
    assert h['failed_linear_solves']==1
    assert not h['linear_solves'][0]['residual_trustworthy']
    with pytest.raises(ValidationError):linear_gate({'exit_code':0,'history':h})
