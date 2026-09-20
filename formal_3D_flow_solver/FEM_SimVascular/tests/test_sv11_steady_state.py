import pytest
from sv_validation.sv11 import load,REPORT
from sv_validation.validation import steady_state,ValidationError

def test_one_saved_state_never_steady():
    with pytest.raises(ValidationError):steady_state([],[])

def test_five_consecutive_intervals_required():
    assert steady_state([1e-6]*5,[1e-7]*5)
    assert not steady_state([1e-6]*4+[1e-4],[1e-7]*5)

def test_actual_full_steady_state():
    if not (REPORT/'petsc_full_qc.json').exists():pytest.skip('No full run; steady solution unavailable')
    assert load('petsc_full_qc')['steady_reached']
