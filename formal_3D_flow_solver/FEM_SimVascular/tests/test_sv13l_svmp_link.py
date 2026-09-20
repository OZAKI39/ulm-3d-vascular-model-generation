from sv13l_support import *
def test_actual_solver_links_rebuilt_stack():
    d=actual('svmp_gpu_link')
    stack_link_gate(d['resolved_MPI'],d['expected_MPI'],'SVMP_MPI')
    stack_link_gate(d['resolved_PETSc'],d['expected_PETSc'],'SVMP_PETSC')
    assert d['source_unmodified']
@pytest.mark.parametrize('label',['SVMP_MPI','SVMP_PETSC'])
def test_wrong_stack_rejected(label):
    with pytest.raises(GateError,match=label):stack_link_gate('/old/lib.so','/new/lib.so',label)
