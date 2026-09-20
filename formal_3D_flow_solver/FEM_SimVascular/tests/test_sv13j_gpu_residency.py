from sv13j_support import *
def test_actual_vascular_gpu_residency():
    d=actual('gpu_residency');assert d['measured']
def test_standalone_smoke_is_not_vascular_residency_evidence():
    assert load('petsc_gpu_smoke')['status']=='PASS'
    d=load('gpu_residency');assert not d['measured'] and d['matrix_resident'] is None and d['ILU2_location'] is None

