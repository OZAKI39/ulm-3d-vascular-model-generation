from sv13j_support import *
def test_actual_mat_is_cuda():
    d=actual('petsc_cuda_types');cuda_types_gate(d['mat_type'],d['vec_type'])
def test_gpu_utilization_does_not_rescue_cpu_matrix():
    d={'mat_type':'seqaij','vec_type':'seqcuda','GPU_utilization':90}
    with pytest.raises(GateError,match='Mat'):cuda_types_gate(d['mat_type'],d['vec_type'])

