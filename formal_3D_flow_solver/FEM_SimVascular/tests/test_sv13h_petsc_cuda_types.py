from sv13h_support import *
def test_actual_runtime_mat_vec_types():
    d=actual('petsc_cuda_types');cuda_types_gate(d['mat_type'],d['vec_type'])
@pytest.mark.parametrize('mat,vec',[('seqaij','seqcuda'),('seqaijcusparse','seq')])
def test_cpu_type_rejected(mat,vec):
    with pytest.raises(GateError,match='PETSC_GPU_TYPES_NOT_ACTIVE'):cuda_types_gate(mat,vec)

