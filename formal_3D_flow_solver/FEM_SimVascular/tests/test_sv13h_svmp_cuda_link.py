from sv13h_support import *
def test_actual_svmp_cuda_link():
    d=actual('svmp_cuda_link');linkage_gate(d['resolved_PETSc'],d['expected_PETSc'])
def test_cpu_petsc_link_rejected():
    with pytest.raises(GateError,match='LINKED_WRONG_PETSC'):
        linkage_gate('/cpu/lib/libpetsc.so','/cuda12/lib/libpetsc.so')

