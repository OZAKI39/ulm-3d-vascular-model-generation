from common import *
import ctypes,sys,time
from mpi4py import MPI
from lammps import lammps
COMM=MPI.COMM_WORLD;rank=COMM.rank
class Native:
 def __init__(self,name,sc=.24,restart=None,legacy=False,contact=False):
  self.legacy=legacy;self.l=lammps(cmdargs=['-log',str(O/(name+'.log')),'-screen','none'],comm=COMM);self.V=DX**3
  self.lib=ctypes.CDLL(str(LIB));self.lib.lpsvar_action.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.c_longlong,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_action.restype=ctypes.c_int
  self.lib.lpsvar_output.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_output.restype=ctypes.c_int
  if restart:
   self.l.commands_string(f'atom_modify map array sort 0 0\nread_restart {restart}')
  else:
   self.l.commands_string(f'''units si
dimension 3
boundary f f f
atom_style bond
atom_modify map array sort 0 0
newton on off
processors {COMM.size} 1 1
read_data {O/'provenance/body.data'} extra/bond/per/atom 150 extra/special/per/atom 200
''')
  self.l.commands_string(f'''lattice sc {DX:.17g}
neighbor {DX:.17g} bin
neigh_modify every 1 delay 0 check no
comm_modify cutoff {3*DELTA:.17g}
fix props all property/atom d_vfrac d_theta ghost yes
set group all d_vfrac {self.V:.17g}
''')
  if not restart:
   if contact:self.l.commands_string(f'pair_style bpm/peri contact_distance_factor 1.0\npair_coeff * * {18*K/(np.pi*DELTA**4):.17g} {DELTA:.17g}')
   else:self.l.commands_string(f'pair_style zero {DELTA:.17g}\npair_coeff * *')
   if legacy:self.l.commands_string(f'bond_style bpm/peri break no\nbond_coeff 1 lps {K:.17g} {G:.17g} {DELTA:.17g} .15 0')
   else:self.l.commands_string(f'bond_style bpm/peri/variational {O/"provenance/reference.native"}\nbond_coeff 1 {K:.17g} {G:.17g} {DELTA:.17g} .15 {sc:.17g}')
  self.l.commands_string('special_bonds lj 0 1 1 coul 1 1 1\nthermo 100000\nthermo_style custom step atoms bonds pe ebond\nrun 0')
  maps=sorted(set(s.split()[-1] for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in s));assert maps and all(sha(p)==sha(LIB) for p in maps)
 def ptr(self,a):return np.ascontiguousarray(a,dtype=np.float64).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
 def field(self,i):
  size={0:1023,1:1023,2:12047,3:12047,4:341,5:24094,6:24094,7:341,8:10,9:12047,10:12047}[i];out=np.empty(size);assert self.lib.lpsvar_output(self.l.lmp,i,self.ptr(out))==size
  return out
 def owned(self):
  n=self.l.extract_global('nlocal');ids=self.l.numpy.extract_atom('id')[:n].copy();return ids
 def gather_atom(self,key):
  ids=self.owned();n=len(ids);val=self.l.numpy.extract_atom(key)[:n].copy();parts=COMM.allgather((ids,val));ii=np.concatenate([a[0] for a in parts]);vv=np.concatenate([a[1] for a in parts]);assert len(np.unique(ii))==341;return vv[np.argsort(ii)]
 def snapshot(self):
  st=dict(x=self.gather_atom('x'),force=self.gather_atom('f')*self.V,owned=self.owned(),bonds=int(self.l.get_thermo('bonds')),ebond_J=float(self.l.get_thermo('ebond'))*self.V)
  if not self.legacy:st.update(kappa=self.field(2),d=self.field(3),m=self.field(4),v=self.field(5),tdir=self.field(6),theta=self.field(7),scalars=self.field(8),bond_ids=self.field(9),failed=self.field(10),force_cpp=self.field(0).reshape(-1,3))
  return st
 def eval(self,x=None):
  if x is not None:self.l.numpy.extract_atom('x')[:len(self.owned())]=x[self.owned()-1]
  self.l.command('run 0');return self.snapshot()
 def import_state(self,kappa,d):
  kappa=np.ascontiguousarray(kappa);d=np.ascontiguousarray(d);z=self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(kappa),self.ptr(d));assert z==0,z
 def commit(self,token):
  z=self.lib.lpsvar_action(self.l.lmp,2,token,None,None);assert z in [0,1],z;self.l.command('run 0');return self.snapshot()
 def close(self):self.l.close()
