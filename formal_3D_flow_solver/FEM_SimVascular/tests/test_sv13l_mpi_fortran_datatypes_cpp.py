from sv13l_support import *
def test_actual_cpp_three_complete_probes():
    d=actual('mpi_fortran_datatypes');assert len(d['repetitions'])==3
    for repeat in d['repetitions']:
        assert repeat['accepted'] and len(repeat['runs'])==2
        for ranks,run in enumerate(repeat['runs'],1):datatype_rows_gate(run['rows'],ranks,d['native_Fortran_sizes'])
def test_zero_size_rejected():
    d=load('mpi_fortran_datatypes');rows=d['repetitions'][0]['runs'][0]['rows']
    next(r for r in rows if r['name']=='MPI_INTEGER')['size']=0
    with pytest.raises(GateError,match='SIZE_INVALID'):datatype_rows_gate(rows,1,d['native_Fortran_sizes'])
