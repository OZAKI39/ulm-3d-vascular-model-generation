from sv13_support import artifact
def test_no_unbounded_gpu_pc_exploration():
    d=artifact('gpu_decision');assert d['gpu_B_attempt_count'] in (0,1)
    assert d['gpu_C_or_later_attempt_count']==0
    if d['gpu_B_attempt_count']:assert d['gpu_A_scientific_pass'] and d['gpu_A_PC_or_transfer_dominated']
