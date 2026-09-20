import pytest
from sv_validation.sv11 import load,REPORT,measured_mass_gate
from sv_validation.validation import ValidationError

def test_mass_from_target_instead_of_measured_solution_rejected():
    with pytest.raises(ValidationError,match='measured'):measured_mass_gate({'source':'Qtarget'})

def test_fabricated_zero_closure_rejected():
    data={'source':'actual_vtu_surface_integration','Q_target_m3_s':1.,'Q_in_m3_s':1.,
          'outlet_flows_m3_s':{'OUTLET_01':.2,'OUTLET_02':.2,'OUTLET_03':.2},'epsilon_mass':0.}
    with pytest.raises(ValidationError,match='inconsistent'):measured_mass_gate(data)

def test_accepted_solution_actual_mass():
    if not (REPORT/'accepted_solution.json').exists():pytest.skip('No accepted steady solution')
    assert measured_mass_gate(load('accepted_solution'))
