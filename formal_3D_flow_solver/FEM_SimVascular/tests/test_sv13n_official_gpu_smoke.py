from sv13n_support import *
def test_official_gpu_complete():flow_gate(accepted('svmp_gpu_smoke'))
@pytest.mark.parametrize('field,value',[('mat_type','seqaij'),('vec_type','seq'),('PETSc_error_detected',True),('MPI_error_detected',True)])
def test_exit_zero_hard_error_or_cpu_backend_rejected(field,value):
 d=good_flow();d[field]=value
 with pytest.raises(GateError):flow_gate(d)
