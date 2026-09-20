from sv13g_support import *
from sv_validation.sv13g import *
def test_verified_singleton():mpi_process_gate(actual('mpi_hard_gate')['singleton'],1)
def test_system_singleton_failure_retained():
    d=load('system_mpi_tests');assert len(d['tests'])>=6
    r=next(x for x in d['tests'] if x['name']=='direct_singleton')
    assert r['timeout'] and r['wall_time_s']<12
