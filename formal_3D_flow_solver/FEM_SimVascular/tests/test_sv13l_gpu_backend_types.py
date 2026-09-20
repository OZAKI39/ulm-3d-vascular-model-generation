from sv13l_support import *
def test_actual_solver_cuda_types():
    d=actual('runtime_backend');cuda_types_gate(d['mat_type'],d['vec_type'])
@pytest.mark.parametrize('mat,vec',[('seqaij','seqcuda'),('seqaijcusparse','seq')])
def test_standalone_gpu_success_cannot_hide_cpu_solver_types(mat,vec):
    with pytest.raises(GateError,match='NOT_ACTIVE'):cuda_types_gate(mat,vec)
