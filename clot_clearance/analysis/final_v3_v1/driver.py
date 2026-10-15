from common import *
import os,sys,time,resource,traceback,ctypes as ct,numpy as np
from mpi4py import MPI
from lammps import lammps
from runtime_audit import E,Flow,Audit
comm=MPI.COMM_WORLD;rank=comm.rank;assert comm.size==2
assert json.loads((O/'STATIC_PREFLIGHT.json').read_text())['status']=='PASS'
for name,entry in json.loads((O/'V3_INPUT_IDENTITY_AUDIT.json').read_text())['files'].items():assert sha(P/name)==entry['sha256']==sha(OLD/name)
os.chdir(P);E.OUT=P;t=time.perf_counter();L=lammps(comm=comm,cmdargs=['-log','log.lammps','-screen','screen.lammps']);wall=flow=audit=None
try:
 dump(P/f'LOADED_LIBRARY_rank{rank}.json',loaded(L))
 prefix,suffix=(P/'in.lammps').read_text().split('\nrun ',1);L.commands_string(prefix)
 wall=E.AuditedWall(L,comm,P,1927);flow=Flow(L,comm,P,wall);audit=Audit(L,comm,flow,wall)
 L.set_fix_external_callback('nativewall',wall.callback);L.set_fix_external_callback('earlyflow',flow.callback)
 L.command('fix audit_readonly all audit/readonly/finalv3')
 dump(P/f'RANK_READY_{rank}.json',{'rank':rank,'pid':os.getpid(),'library_sha256':EXPECTED,'observer_sha256':sha(O/'observer.so'),'readonly_fix_command':'fix audit_readonly all audit/readonly/finalv3','input_in_lammps_byte_identical':True})
 comm.Barrier();L.commands_string('run '+suffix)
 assert L.extract_global('ntimestep')==4660 and audit.lastcomplete==4660
 x=np.ctypeslib.as_array(L.gather_atoms('x',1,3)).reshape(-1,3).copy();v=np.ctypeslib.as_array(L.gather_atoms('v',1,3)).reshape(-1,3).copy()
 if rank==0:np.savez_compressed(P/'terminal_full_precision.npz',position_m=x,velocity_m_s=v)
 flow.close();wall.close();audit.close();L.close()
 q=resource.getrusage(resource.RUSAGE_SELF);dump(P/f'RANK_FINISHED_{rank}.json',{'rank':rank,'accepted_step':4660,'time_s':4660*5e-7,'wall_s':time.perf_counter()-t,'cpu_user_s':q.ru_utime,'cpu_system_s':q.ru_stime,'max_RSS_KiB':q.ru_maxrss})
 comm.Barrier();MPI.Finalize()
except BaseException:
 dump(P/f'DRIVER_ERROR_rank{rank}.json',{'last_completed':audit.lastcomplete if audit else None,'traceback':traceback.format_exc()});comm.Abort(61)
