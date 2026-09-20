from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_runtime_cuda_vec():
    d=actual('cuda_types');cuda_types_gate(d['mat_type'],d['vec_type'])
def test_cpu_vector_rejected():
    with pytest.raises(GateError,match='Vec'):cuda_types_gate('seqaijcusparse','seq')
