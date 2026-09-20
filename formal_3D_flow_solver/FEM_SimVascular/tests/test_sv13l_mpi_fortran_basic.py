from sv13l_support import *
def test_actual_native_fortran():
    d=actual('mpi_fortran_basic');assert len(d['runs'])==2
    for i,r in enumerate(d['runs'],1):
        assert r['exit_code']==0 and len(r['rows'])==i
        assert all(x['failures']==0 and all(s>0 for s in x['sizes']) for x in r['rows'])
def test_c_success_cannot_hide_fortran_failure():
    d=load('mpi_application_gate');d['Fortran_rank2']='FAIL'
    with pytest.raises(GateError,match='FORTRAN_BASIC_FAILED'):mpi_application_gate(d)
