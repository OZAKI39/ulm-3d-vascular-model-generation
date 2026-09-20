from sv_validation.sv12 import steady_gate
from sv12_support import artifact

def test_velocity_pass_flow_fail_is_not_steady():assert not steady_gate([{'E_u':1e-8,'E_Q':1e-3}]*5)
def test_step400_with_three_intervals_not_steady():assert not steady_gate([{'step':400,'E_u':1e-8,'E_Q':1e-8}]*3)
def test_only_one_passing_interval_not_steady():assert not steady_gate([{'E_u':1.,'E_Q':1.}]*5+[{'E_u':0.,'E_Q':0.}])
def test_actual_five_joint_intervals():
    artifact('flow_execution');qc=artifact('saved_state_qc')
    assert steady_gate(qc['intervals']) and qc['consecutive_passing_intervals']>=5
