from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_solver_cuda_linkage():
    d=actual('svmp_gpu_build');linkage_gate(d['resolved_petsc'],d['expected_petsc']);assert d['source_unmodified']
def test_cpu_petsc_link_rejected():
    with pytest.raises(GateError,match='WRONG_PETSC'):linkage_gate('/cpu/libpetsc.so','/cuda/libpetsc.so')
