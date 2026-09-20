from sv13l_support import *
def test_actual_matching_cpu_proof():
    d=actual('cpu_proof');proof_gate(d);assert d['initial_state']=='t=0'
def test_rank_mismatched_checkpoint_rejected():
    with pytest.raises(GateError):restart_gate('native',4,1)
    restart_gate('t=0',None,1)

