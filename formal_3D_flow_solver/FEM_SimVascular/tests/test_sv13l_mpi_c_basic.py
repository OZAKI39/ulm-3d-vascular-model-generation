from sv13l_support import *
def test_actual_c_rank_tests():
    d=actual('mpi_c_basic');assert len(d['rank1'])==5
    for r in d['rank1']:mpi_process_gate(r,1)
    mpi_process_gate(d['rank2'],2)
