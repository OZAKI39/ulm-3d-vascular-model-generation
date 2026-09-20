import pytest
from sv_validation.sv11 import load,REPORT

def test_actual_step_ten_mass_strictly_improves():
    if not (REPORT/'petsc_short_qc.json').exists():pytest.skip('Short run unavailable')
    states=load('petsc_short_qc')['states']
    if not states or states[-1]['step']!=10:pytest.skip('Stopped before step 10; no comparable saved VTU')
    assert states[-1]['source']=='actual_vtu_surface_integration'
    assert states[-1]['epsilon_mass']<load('flow_qc','sv1')['epsilon_mass']
