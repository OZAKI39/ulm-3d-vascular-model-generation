import pytest
from sv13_support import gpu_evidence
from sv_validation.sv13 import gpu_gate
@pytest.mark.parametrize('types',[('seqaij','seqcuda'),('seqaijcusparse','seq'),('mpiaij','mpi')])
def test_cpu_types_rejected_even_with_gpu_utilization(types):
    assert gpu_gate({'mat_type':types[0],'vec_type':types[1],'gpu_utilization':99})=='GPU_BACKEND_NOT_REACHABLE'
def test_actual_mat_and_vec_are_cuda():
    d=gpu_evidence('cuda_types')
    assert d['mat_type'] in ('seqaijcusparse','mpiaijcusparse')
    assert d['vec_type'] in ('seqcuda','mpicuda')
