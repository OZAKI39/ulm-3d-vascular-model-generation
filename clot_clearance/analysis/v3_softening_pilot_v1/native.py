from common import *
import ctypes,time
from mpi4py import MPI
from lammps import lammps
from external import External,FIELDS
class Native:
 def __init__(self,name,restart=None,restart_state=None,dt=BASE,query=None):
  self.comm=MPI.COMM_WORLD;assert self.comm.size==1
  self.ref=dict(np.load(O/'provenance/reference.npz'));self.N=len(self.ref['X']);self.n=int(self.ref['clot_count']);self.NB=len(self.ref['bonds']);self.dt=dt;self.dynamic=False;self.restart_state=restart_state;self.query=query
  self.l=lammps(cmdargs=['-log',str(O/(name+'.log')),'-screen','none'],comm=self.comm)
  self.lib=ctypes.CDLL(str(LIB));self.lib.lpsvar_action.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.c_longlong,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_action.restype=ctypes.c_int
  self.lib.lpsvar_output.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_output.restype=ctypes.c_int
  # Setting special_bonds BEFORE read_data avoids the unused three-hop table.
  # It does not alter original graph/m or the active 0/1/1 contact exclusions.
  if restart:self.l.commands_string(f'atom_modify map array sort 0 0\nread_restart {restart}')
  else:self.l.commands_string(f'''units si
boundary m m m
atom_style bond
atom_modify map array sort 0 0
newton on off
special_bonds lj 0 1 1 coul 1 1 1
read_data {O/'provenance/body.data'} extra/bond/per/atom 160 extra/special/per/atom 180
''')
  self.l.commands_string(f'''lattice sc {DX:.17g}
neighbor 0.0001 bin
neigh_modify every 1 delay 0 check no
comm_modify cutoff {3*HORIZON:.17g}
fix props all property/atom d_vfrac d_theta ghost yes
set group all d_vfrac {V:.17g}
''')
  if not restart:
   self.l.commands_string(f'''pair_style bpm/peri contact_distance_factor 1.0
pair_coeff * * {18*K/(np.pi*HORIZON**4):.17g} {HORIZON:.17g}
bond_style bpm/peri/variational {O/'provenance/reference.native'}
''')
   for t in [1,2,3]:self.l.command(f'bond_coeff {t} {K:.17g} {G:.17g} {HORIZON:.17g} {S0[t]:.17g} {SC[t]:.17g}')
  self.l.commands_string(f'''special_bonds lj 0 1 1 coul 1 1 1
timestep {dt:.17g}
thermo 100000
thermo_style custom step time atoms bonds pe ebond epair ke
thermo_modify lost error norm no
''')
  self.maps=sorted(set(s.split()[-1] for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in s));assert self.maps and all(sha(p)==sha(LIB) for p in self.maps)
 def ptr(self,a):return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
 def field(self,i):
  size={0:3*self.N,1:3*self.N,2:self.NB,3:self.NB,4:self.N,5:2*self.NB,6:2*self.NB,7:self.N,8:10,9:self.NB,10:self.NB,11:5,12:self.NB,13:self.NB,14:self.NB}[i];out=np.empty(size);assert self.lib.lpsvar_output(self.l.lmp,i,self.ptr(out))==size;return out
 def gather(self,key):
  n=self.l.extract_global('nlocal');ids=self.l.numpy.extract_atom('id')[:n].copy();vals=self.l.numpy.extract_atom(key)[:n].copy();assert n==self.N and np.array_equal(np.sort(ids),np.arange(1,self.N+1));return vals[np.argsort(ids)]
 def state(self):
  if self.dynamic and self.ext.error:raise RuntimeError(self.ext.error)
  d=self.field(3);td=self.field(6);failed=d==1
  z=dict(step=int(self.l.get_thermo('step')),time=float(self.l.get_thermo('time')),x=self.gather('x'),velocity=self.gather('v'),force=self.gather('f')*(1 if self.dynamic else V),material_force=self.field(0).reshape(-1,3),kappa=self.field(2),d=d,scalars=self.field(8),counters=self.field(11),native_bonds=int(self.l.get_thermo('bonds')),Uc=float(self.l.get_thermo('epair'))*V,K=float(self.l.get_thermo('ke')),native_U=float(self.l.get_thermo('ebond'))*V,failed_tension_max_N=float(np.max(abs(td[np.r_[failed,failed]]),initial=0)),atom_ids=self.gather('id').astype(np.int64),box_lo=np.asarray(self.l.extract_box()[0]),box_hi=np.asarray(self.l.extract_box()[1]))
  if self.dynamic:z.update({k:v.copy() for k,v in self.ext.last.items()})
  return z
 def static(self,x=None):
  if x is not None:
   ids=self.l.numpy.extract_atom('id')[:self.N];self.l.numpy.extract_atom('x')[:self.N]=x[ids-1]
  self.l.command('run 0');return self.state()
 def dynamics(self):
  self.l.commands_string('group clot type 1\ngroup anchors type 2\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC '+format(RHO,'.17g')+'\nfix v3forces clot external pf/callback 1 1')
  self.ext=External(self,self.restart_state,self.query);self.l.set_fix_external_callback('v3forces',self.ext.callback);self.dynamic=True;self.l.command('run 0');return self.state()
 def close(self):self.l.close()
class Contact:
 def __init__(self,ref):
  from scipy.spatial import cKDTree
  self.tree=cKDTree;self.b=ref['bonds'];self.n=len(ref['X']);self.codes=self.b[:,0]*self.n+self.b[:,1];self.c=15*(18*K/(np.pi*HORIZON**4))*V**2/HORIZON
 def evaluate(self,x,d):
  pairs=self.tree(x).query_pairs(DX,output_type='ndarray');pairs=pairs[~np.isin(pairs[:,0]*self.n+pairs[:,1],self.codes[d<1])];f=np.zeros_like(x)
  if not len(pairs):return 0.,f,pairs
  dv=x[pairs[:,1]]-x[pairs[:,0]];r=np.linalg.norm(dv,axis=1);dr=r-DX;ff=self.c*dr[:,None]*dv/r[:,None]
  for c in range(3):f[:,c]=np.bincount(pairs[:,0],weights=ff[:,c],minlength=self.n)-np.bincount(pairs[:,1],weights=ff[:,c],minlength=self.n)
  return float(.5*self.c*np.dot(dr,dr)),f,pairs
