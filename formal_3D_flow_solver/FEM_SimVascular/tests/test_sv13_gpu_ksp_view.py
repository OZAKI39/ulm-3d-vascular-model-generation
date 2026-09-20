from sv13_support import gpu_evidence
def test_actual_ksp_pc_view_retained():
    r=gpu_evidence('gpu_ksp_view')
    assert r['ksp_type']=='gmres' and r['rtol']==1e-10 and r['atol']==1e-24
    assert r['view_log_sha256'] and r['actual_pc_type']
