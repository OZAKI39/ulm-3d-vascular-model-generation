import csv,math
import pytest
from sv_validation.sv12 import OUTPUT
from sv_validation.validation import mass_balance,ValidationError
from sv12_support import artifact

def test_saved_mass_history_uses_every_actual_state():
    states=artifact('saved_state_qc')['states']
    rows=list(csv.DictReader((OUTPUT/'qc/mass_history.csv').open()))
    assert len(rows)==len(states)
    for row,state in zip(rows,states):
        assert row['source']=='actual_vtu_surface_integration'
        assert math.isclose(float(row['epsilon_mass']),abs(float(row['Qout_total'])-float(row['Qin']))/float(row['Qtarget']),rel_tol=1e-12,abs_tol=1e-16)
        assert int(row['step'])==state['step']

def test_nan_outlet_flow_rejected():
    with pytest.raises(ValidationError):mass_balance(1.,-1.,{'OUTLET_01':float('nan'),'OUTLET_02':.4,'OUTLET_03':.6})
