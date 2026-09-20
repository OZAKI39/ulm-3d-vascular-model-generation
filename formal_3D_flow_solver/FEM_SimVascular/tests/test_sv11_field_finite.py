import numpy as np
import pytest
from sv_validation.sv11 import load,REPORT,ROOT
from sv_validation.validation import parse_result_vtu,fields_finite,ValidationError

def test_nonfinite_field_rejected():
    with pytest.raises(ValidationError):fields_finite(np.array([[0.,np.nan,0.]]),np.array([0.]))

def test_short_actual_fields_finite():
    if not (REPORT/'petsc_short_qc.json').exists():pytest.skip('No vascular saved state')
    states=load('petsc_short_qc')['states']
    if not states:pytest.skip('Stopped before first scheduled VTU')
    for state in states:parse_result_vtu(ROOT/state['path'])
