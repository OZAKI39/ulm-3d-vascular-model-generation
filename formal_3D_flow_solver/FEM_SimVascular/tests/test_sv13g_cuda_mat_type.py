from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_runtime_cuda_mat():
    d=actual('cuda_types');cuda_types_gate(d['mat_type'],d['vec_type'])
def test_cpu_aij_rejected_despite_configure_success_and_gpu_load():
    evidence={'configure_success':True,'gpu_utilization':90,'mat_type':'seqaij','vec_type':'seqcuda'}
    with pytest.raises(GateError,match='Mat'):cuda_types_gate(evidence['mat_type'],evidence['vec_type'])
