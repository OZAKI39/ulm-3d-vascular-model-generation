from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_proof_equivalence():equivalence_gate(actual('gpu_cpu_equivalence')['metrics'],json.loads((ROOT/'configs/sv1_3g/policy.json').read_text()))
def test_changed_solution_rejected():
    p=json.loads((ROOT/'configs/sv1_3g/policy.json').read_text())
    m=dict(velocity_relative_L2=.01,pressure_relative_L2=0,Qin_relative=0,mass_error_difference=0,Qout_relative={'OUTLET_01':0})
    with pytest.raises(GateError,match='NUMERICAL_DIFFERENCE'):equivalence_gate(m,p)
