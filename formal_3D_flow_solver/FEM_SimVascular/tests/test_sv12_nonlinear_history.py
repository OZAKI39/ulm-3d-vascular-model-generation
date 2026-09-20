import copy
import pytest
from sv_validation.sv11 import nonlinear_gate
from sv_validation.sv12 import load
from sv_validation.validation import ValidationError
from sv12_support import artifact

def test_current_nonlinear_history():
    flow=artifact('flow_execution')
    for name in flow['runs']:
        data=artifact(name+'_execution');assert data['nonlinear_failures']==0
        assert all(r['converged'] and 2<=r['iterations']<=12 for r in data['nonlinear_history'])

def test_iteration_cap_without_residual_convergence_rejected():
    history=copy.deepcopy(load('petsc_short_execution','sv1_1')['history'])
    history['linear_solves'][-1].update(nonlinear_iteration=12,nonlinear_Ri_over_R0=.1,nonlinear_Ri_over_R1=.1)
    with pytest.raises(ValidationError):nonlinear_gate(history)
