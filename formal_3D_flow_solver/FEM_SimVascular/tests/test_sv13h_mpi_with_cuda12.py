from sv13h_support import *
def test_actual_five_rank1_plus_rank2_with_unchanged_mpi():
    d=actual('mpi_with_cuda12')
    assert len(d['rank1'])==5 and not d['MPI_rebuilt'] and not d['MPI_wrapper_modified']
    for r in d['rank1']:mpi_process_gate(r,1)
    mpi_process_gate(d['rank2'],2)
    a=actual('petsc_failure_audit')
    assert not a['MPI_changes'] and a['MPI_wrapper_unchanged']
@pytest.mark.parametrize('change',[{'exit_code':127},{'timeout':True},{'stdout':'libmpi.so: not found'}])
def test_cuda_library_path_breaking_mpi_is_rejected(change):
    r=copy.deepcopy(load('mpi_with_cuda12')['rank1'][0]);r.update(change)
    with pytest.raises(GateError):mpi_process_gate(r,1)

