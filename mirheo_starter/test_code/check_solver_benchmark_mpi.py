"""CPU MPI transport test, not a fluid measurement and never imports Mirheo."""
import json
import time
from py_scripts.solver_benchmark.native_mpi import NativeMPI


def main():
    m=NativeMPI();m.check(m.lib.MPI_Init(None,None))
    try:
        m.Barrier();t=time.perf_counter()
        if m.rank==1:time.sleep(.025)
        elapsed=m.allreduce(time.perf_counter()-t)
        assert elapsed>=.02
        assert m.allreduce(m.rank+1,op=m.SUM)==3
        assert m.bcast(True if m.rank==0 else None)
        assert not m.bcast(False if m.rank==0 else None)
        m.Barrier()
        if m.rank==0:print(json.dumps({'status':'PASS','ranks':2,'maximum_rank_elapsed_s':elapsed,'sum':3,'category':'CPU_TRANSPORT_TEST_NOT_FLUID_DATA'}))
    finally:m.check(m.lib.MPI_Finalize())


if __name__=='__main__':main()
