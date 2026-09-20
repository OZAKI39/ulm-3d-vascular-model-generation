from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_transfer_scaling():
    d=actual('gpu_transfer');assert d['measured'] and not d['large_transfer_per_iteration']
def test_large_transfer_every_iteration_rejected():
    d=dict(measured=True,mat_type='seqaijcusparse',vec_type='seqcuda',matrix_resident=True,vectors_resident=True,pc_apply_device='GPU',large_transfer_per_iteration=True)
    with pytest.raises(GateError,match='EXCESSIVE_TRANSFER'):residency_gate(d)
def test_cpu_pc_apply_not_called_fully_gpu_resident():
    d=dict(measured=True,mat_type='seqaijcusparse',vec_type='seqcuda',matrix_resident=True,vectors_resident=True,pc_apply_device='CPU',large_transfer_per_iteration=False)
    with pytest.raises(GateError,match='PC_NOT_RESIDENT'):residency_gate(d)
