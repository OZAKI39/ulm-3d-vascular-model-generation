from sv13_support import gpu_evidence
def test_actual_residency_evidence():
    r=gpu_evidence('gpu_residency_audit')
    assert r['matrix_resident'] and r['vectors_resident']
    assert r['transfer_scales_with']!='KSP_iterations'
