from sv13l_support import *
def test_actual_application_gate():mpi_application_gate(actual('mpi_application_gate'))
@pytest.mark.parametrize('field,value,reason',[('bcast_rc',3,'BCAST_FAILED'),('data_ok',0,'PAYLOAD_MISMATCH')])
def test_type_size_success_not_sufficient(field,value,reason):
    d=load('mpi_fortran_datatypes');rows=d['repetitions'][0]['runs'][0]['rows']
    next(r for r in rows if r['name']=='MPI_INTEGER')[field]=value
    with pytest.raises(GateError,match=reason):datatype_rows_gate(rows,1,d['native_Fortran_sizes'])
