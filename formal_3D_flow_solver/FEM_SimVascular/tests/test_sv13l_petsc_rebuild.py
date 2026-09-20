from sv13l_support import *
def test_actual_clean_petsc_rebuild():
    d=actual('petsc_rebuild');assert d['version']=='3.19.6' and d['CUDA_version']=='12.3.2'
    assert d['configure']['exit_code']==d['make']['exit_code']==0 and d['source_unmodified']
    assert len(d['configuration_diff_from_J'])==5
    source_integrity_gate(load('petsc_source_integrity_before'))
    assert all(actual('petsc_self_test')['observed'].values())
def test_petsc_links_new_mpi():
    d=actual('petsc_mpi_link');stack_link_gate(d['resolved_MPI'],d['expected_MPI'],'PETSC_MPI')
def test_petsc_old_mpi_rejected():
    with pytest.raises(GateError,match='PETSC_MPI'):stack_link_gate('/old/gpu_mpi/libmpi.so','/new/gpu_mpi_fortran/libmpi.so','PETSC_MPI')
