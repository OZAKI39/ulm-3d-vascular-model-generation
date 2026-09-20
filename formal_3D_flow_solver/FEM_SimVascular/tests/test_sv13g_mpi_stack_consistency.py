from sv13g_support import *
from sv_validation.sv13g import *
def test_actual_compiler_runtime_stack():
    d=load('mpi_stack_consistency');assert d['status']=='PASS'
    assert 'libmpi.so' in d['binary_ldd']['stdout']
    assert 'not found' not in d['binary_ldd']['stdout']
    assert 'libopen-rte' in d['launcher_ldd']['stdout']
def test_different_implementation_rejected():
    with pytest.raises(GateError,match='MPI_STACK_MISMATCH'):mpi_stack_gate('MPICH','Open MPI','/x','/x')
def test_different_prefix_rejected():
    with pytest.raises(GateError,match='MPI_STACK_MISMATCH'):mpi_stack_gate('Open MPI','Open MPI','/a','/b')
