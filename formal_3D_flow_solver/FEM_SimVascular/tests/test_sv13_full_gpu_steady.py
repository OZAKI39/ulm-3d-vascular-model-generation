from sv13_support import gpu_evidence
def test_actual_gpu_full_run_gates():
    d=gpu_evidence('gpu_full_run')
    assert d['from_t0'] and d['automatic_steady_stop'] and d['linear_failures']==d['nonlinear_failures']==0
    assert d['mass_pass'] and d['steady_pass'] and d['wall_pass'] and d['field_finite']
