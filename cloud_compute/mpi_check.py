"""Actual two-rank collectives through the archived, locally tested OpenMPI ABI wrapper."""
import sys,json
sys.path.insert(0,sys.argv[1])
from py_scripts.solver_benchmark.native_mpi import NativeMPI
mpi=NativeMPI();mpi.check(mpi.lib.MPI_Init(None,None))
assert mpi.allreduce(mpi.rank+1,mpi.SUM)==3.0,'ALLREDUCE_FAILED'
assert mpi.bcast(True if mpi.rank==0 else None),'BROADCAST_FAILED'
mpi.Barrier()
print(json.dumps(dict(rank=mpi.rank,allreduce_sum=3.0,broadcast=True,barrier=True,status='MPI_COMMUNICATION_PASS')),flush=True)
mpi.check(mpi.lib.MPI_Finalize())
