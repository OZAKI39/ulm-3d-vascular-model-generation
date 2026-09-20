from sv13n_support import *
def test_cpu_mpi_ghost():ghost_run_gate(ghost('CPU_MPI'))
def test_optional_cuda_mpi_failure_is_retained_and_rejected():
 d=read('ghost_probe_new');r=ghost('CUDA_MPI')
 assert d['status']=='FAIL' and not r['accepted']
 with pytest.raises(GateError):ghost_run_gate(r)
 assert read('gpu_ghost_target')['MPI_CUDA_supported'] is False
