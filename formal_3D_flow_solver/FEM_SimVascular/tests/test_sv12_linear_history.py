import copy
import pytest
from sv_validation.sv11 import linear_gate
from sv_validation.sv12 import load
from sv_validation.validation import ValidationError
from sv12_support import artifact

def test_every_current_linear_solve_success():
    flow=artifact('flow_execution')
    for name in flow['runs']:
        run=artifact(name+'_execution')
        assert run['linear_failures']==0 and run['ill_conditioned_warnings']==0 and run['exit_code']==0
        assert run['linear_solves'][0]['step']==(11 if name=='production' else 401)
        for row in run['linear_solves']:
            assert row['linear_converged'] and row['residual_finite']
            assert row['petsc_reason']['reason'].startswith('CONVERGED_')
            assert row['petsc_monitor'] and all(r['finite'] for r in row['petsc_monitor'])

def test_exit_zero_with_divergence_rejected():
    run=copy.deepcopy(load('petsc_short_execution','sv1_1'));run['history']['petsc_reasons'][0]['diverged']=True
    assert run['exit_code']==0
    with pytest.raises(ValidationError):linear_gate(run)
