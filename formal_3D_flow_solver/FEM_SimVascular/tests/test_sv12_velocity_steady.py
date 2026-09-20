from sv_validation.sv12 import steady_gate
from sv12_support import artifact

def test_mass_pass_does_not_override_velocity_unsteady():
    rows=[{'E_u':1e-3,'E_Q':1e-8,'epsilon_mass':0.,'step':350+10*i} for i in range(6)]
    assert not steady_gate(rows)

def test_actual_velocity_last_five_intervals():
    state=artifact('flow_execution');qc=artifact('saved_state_qc')
    assert len(qc['intervals'])>=5 and all(r['E_u']<=1e-5 for r in qc['intervals'][-5:])
