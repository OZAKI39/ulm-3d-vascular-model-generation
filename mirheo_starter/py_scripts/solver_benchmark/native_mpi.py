"""Three OpenMPI collectives through the already installed C ABI; no mpi4py.

OpenMPI-specific exported handles are checked by a real two-rank CPU test.
Mirheo retains MPI initialization/finalization ownership in the fluid worker.
"""
import ctypes as C
import os


class NativeMPI:
    MAX='max'
    SUM='sum'
    def __init__(self):
        self.lib=C.CDLL('/usr/lib/x86_64-linux-gnu/libmpi.so.40')
        self.rank=int(os.environ.get('OMPI_COMM_WORLD_RANK','0'))
        def handle(name):return C.c_void_p(C.addressof(C.c_byte.in_dll(self.lib,name)))
        self.world=handle('ompi_mpi_comm_world');self.double=handle('ompi_mpi_double');self.integer=handle('ompi_mpi_int')
        self.ops={self.MAX:handle('ompi_mpi_op_max'),self.SUM:handle('ompi_mpi_op_sum')}
        self.lib.MPI_Allreduce.argtypes=[C.c_void_p,C.c_void_p,C.c_int,C.c_void_p,C.c_void_p,C.c_void_p]
        self.lib.MPI_Bcast.argtypes=[C.c_void_p,C.c_int,C.c_void_p,C.c_int,C.c_void_p]
        self.lib.MPI_Barrier.argtypes=[C.c_void_p]
    def check(self,rc):
        if rc:raise RuntimeError('OPENMPI_C_ABI_ERROR '+str(rc))
    def Barrier(self):self.check(self.lib.MPI_Barrier(self.world))
    def allreduce(self,value,op=MAX):
        source=C.c_double(value);dest=C.c_double()
        self.check(self.lib.MPI_Allreduce(C.byref(source),C.byref(dest),1,self.double,self.ops[op],self.world))
        return dest.value
    def bcast(self,value,root=0):
        buf=C.c_int(bool(value) if value is not None else 0)
        self.check(self.lib.MPI_Bcast(C.byref(buf),1,self.integer,root,self.world))
        return bool(buf.value)
