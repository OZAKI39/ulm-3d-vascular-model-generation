from sv13g_support import *
from sv_validation.sv13g import *
def test_five_consecutive_actual_rank1_runs():
    d=actual('mpi_hard_gate');assert len(d['rank1'])==5
    for r in d['rank1']:mpi_process_gate(r,1)
    assert len({tuple(r['command']) for r in d['rank1']})==1
def test_mpiexec_timeout_is_failure():
    with pytest.raises(GateError,match='MPI_RUNTIME_UNUSABLE'):mpi_process_gate({'exit_code':0,'timeout':True,'stdout':'rank=0 size=1'},1)
def test_actual_hang_was_bounded_and_rejected():
    d=load('system_mpi_tests');assert d['status']=='FAIL'
    for r in d['tests']:
        if r['timeout']:
            assert r['timeout_s']==10 and r['wall_time_s']<12
            with pytest.raises(GateError,match='MPI_RUNTIME_UNUSABLE'):mpi_process_gate(r,1)
