from sv13j_support import *
def test_actual_official_fluid_gpu_smoke_acceptance():
    # Actual observed failure remains a failing acceptance test, never xfailed.
    d=actual('svmp_gpu_smoke')
    assert d['execution']['exit_code']==0 and d['VTU_count']>0
def test_native_datatype_probe_confirms_application_gap():
    d=load('mpi_datatype_diagnosis');assert d['status']=='CONFIRMED'
    assert not d['MPI_changed'] and not d['solver_changed']
    for run in d['runs']:
        assert not run['timeout']
        with pytest.raises(GateError,match='DATATYPE_UNAVAILABLE'):application_mpi_gate(run['rows'])
        assert all(r['bcast_rc']==0 for r in run['rows'] if r['name'] in ('MPI_INT','MPI_DOUBLE','MPI_CHAR','MPI_CXX_BOOL'))
def test_failed_smoke_did_not_reach_backend_or_generate_field():
    d=load('svmp_gpu_smoke')
    assert d['runtime_cuda_types']=='NOT_REACHED' and d['VTU_count']==0
    assert not d['linear_solves_reached'] and not d['nonlinear_solves_reached']

