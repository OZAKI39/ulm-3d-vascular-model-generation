from sv13h_support import *
def test_actual_matched_cpu_twenty_steps():
    d=actual('cpu_proof');proof_gate(d);assert d['initial_state']=='t=0'
def test_four_rank_checkpoint_cannot_seed_one_rank_proof():
    with pytest.raises(GateError,match='INCOMPATIBLE_NATIVE_RESTART'):restart_gate('native',4,1)
    restart_gate('t=0',None,1)

