from sv13j_support import *
def test_actual_vec_is_cuda():
    d=actual('petsc_cuda_types');assert d['vec_type']=='seqcuda'
def test_cpu_vector_rejected():
    with pytest.raises(GateError,match='Vec'):cuda_types_gate('seqaijcusparse','seq')

