from common import *
import ctypes,sys,time,os
from mpi4py import MPI
from lammps import lammps
sys.path.insert(0,str(R/'analysis/softening_feasibility_v1'))
from model import SofteningLPS
COMM=MPI.COMM_WORLD
class Native:
 def __init__(self,name,root=None,sc=SC,restart=None,dt=.5e-6,box=BOX):
  self.root=Path(root) if root else O/'provenance';self.ref=np.load(self.root/'reference.npz');self.N=len(self.ref['X']);self.NB=len(self.ref['bonds']);self.dt=dt;self.box=box;self.deferred_box=bool(restart)
  self.l=lammps(cmdargs=['-log',str(O/(name+'.log')),'-screen','none'],comm=COMM)
  self.lib=ctypes.CDLL(str(LIB));self.lib.lpsvar_action.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.c_longlong,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_action.restype=ctypes.c_int
  self.lib.lpsvar_output.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_output.restype=ctypes.c_int
  if restart:self.l.commands_string(f'atom_modify map array sort 0 0\nread_restart {restart}')
  else:self.l.commands_string(f'''units si
boundary f f f
atom_style bond
atom_modify map array sort 0 0
newton on off
processors {COMM.size} 1 1
read_data {self.root/'body.data'} extra/bond/per/atom 150 extra/special/per/atom 200
''')
  if not restart:self.set_box()
  self.l.commands_string(f'''lattice sc {DX:.17g}
neighbor {DX:.17g} bin
neigh_modify every 1 delay 0 check no
comm_modify cutoff {3*HORIZON:.17g}
fix props all property/atom d_vfrac d_theta ghost yes
set group all d_vfrac {V:.17g}
''')
  if not restart:self.l.commands_string(f'''pair_style bpm/peri contact_distance_factor 1.0
pair_coeff * * {18*K/(np.pi*HORIZON**4):.17g} {HORIZON:.17g}
bond_style bpm/peri/variational {self.root/'reference.native'}
bond_coeff 1 {K:.17g} {G:.17g} {HORIZON:.17g} {S0:.17g} {sc:.17g}
''')
  self.l.commands_string(f'''special_bonds lj 0 1 1 coul 1 1 1
timestep {dt:.17g}
thermo 100000
thermo_style custom step time atoms bonds pe ebond epair ke
thermo_modify lost error norm no
''')
  maps=sorted(set(s.split()[-1] for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in s));assert maps and all(sha(p)==sha(LIB) for p in maps),maps
  if not restart:self.import_precrack()
 def ptr(self,a):return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
 def field(self,i):
  size={0:3*self.N,1:3*self.N,2:self.NB,3:self.NB,4:self.N,5:2*self.NB,6:2*self.NB,7:self.N,8:10,9:self.NB,10:self.NB,11:5}[i];out=np.empty(size);assert self.lib.lpsvar_output(self.l.lmp,i,self.ptr(out))==size;return out
 def gather(self,key):
  n=self.l.extract_global('nlocal');ids=self.l.numpy.extract_atom('id')[:n].copy();vals=self.l.numpy.extract_atom(key)[:n].copy();chunks=COMM.allgather((ids,vals));ii=np.concatenate([a[0] for a in chunks]);vv=np.concatenate([a[1] for a in chunks]);assert len(np.unique(ii))==self.N;return vv[np.argsort(ii)]
 def state(self):
  d=self.field(3);td=self.field(6);failed=d==1
  return dict(step=int(self.l.get_thermo('step')),time=float(self.l.get_thermo('time')),x=self.gather('x'),velocity=self.gather('v'),force=self.gather('f'),material_force=self.field(0).reshape(-1,3),kappa=self.field(2),d=d,scalars=self.field(8),counters=self.field(11),native_bonds=int(self.l.get_thermo('bonds')),Uc=float(self.l.get_thermo('epair'))*V,K=float(self.l.get_thermo('ke')),native_U=float(self.l.get_thermo('ebond'))*V,failed_tension_max_N=float(np.max(abs(td[np.r_[failed,failed]]),initial=0)),atom_ids=self.gather('id').astype(np.int64),box_lo=np.asarray(self.l.extract_box()[0]),box_hi=np.asarray(self.l.extract_box()[1]))
 def import_precrack(self):
  cut=self.ref['initial_cut'];d=cut.astype(float);k=SC*d
  assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(k),self.ptr(d))==0
 def static(self,x=None,kappa=None,d=None):
  if x is not None:self.l.numpy.extract_atom('x')[:self.N]=x
  if d is not None:
   assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(np.ascontiguousarray(kappa)),self.ptr(np.ascontiguousarray(d)))==0
  self.l.command('run 0');z=self.state();z['force']*=V;return z
 def dynamics(self,bottom,top,free):
  for name,ids in [('bottom',bottom),('top',top),('free',free)]:self.l.command('group '+name+' id '+' '.join(map(str,ids)))
  self.l.commands_string(f'''variable q equal ternary(time<{RAMP:.17g},time/{RAMP:.17g},1)
variable disp equal {AMPLITUDE/2:.17g}*(1-cos(PI*v_q))
variable speed equal {AMPLITUDE*np.pi/(2*RAMP):.17g}*sin(PI*v_q)
variable negdisp equal -v_disp
variable negspeed equal -v_speed
fix integrate free nve
fix motion_bottom bottom move variable NULL NULL v_negdisp NULL NULL v_negspeed
fix motion_top top move variable NULL NULL v_disp NULL NULL v_speed
fix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}
run 0
''')
  if self.deferred_box:
   changed=self.set_box();self.deferred_box=False
   if changed:self.l.command('run 0')
 def set_box(self):
  lo,hi,*_=self.l.extract_box()
  if self.box is None or (np.array_equal(lo,-self.box) and np.array_equal(hi,self.box)):return False
  before=self.gather('x');tags=self.gather('id')
  self.l.command('change_box all '+ ' '.join(f'{axis} final {-q:.17g} {q:.17g}' for axis,q in zip('xyz',self.box))+' units box')
  assert np.array_equal(before,self.gather('x')) and np.array_equal(tags,self.gather('id'))
  return True
 def close(self):self.l.close()

class Contact:
 def __init__(self,ref):
  from scipy.spatial import cKDTree
  self.tree=cKDTree;self.b=ref['bonds'];self.n=len(ref['X']);self.codes=self.b[:,0]*self.n+self.b[:,1];self.c=15*(18*K/(np.pi*HORIZON**4))*V**2/HORIZON
 def evaluate(self,x,d):
  pairs=self.tree(x).query_pairs(DX,output_type='ndarray');pairs=pairs[~np.isin(pairs[:,0]*self.n+pairs[:,1],self.codes[d<1])];f=np.zeros_like(x)
  if not len(pairs):return 0.,f,0
  dv=x[pairs[:,1]]-x[pairs[:,0]];r=np.linalg.norm(dv,axis=1);dr=r-DX;ff=self.c*dr[:,None]*dv/r[:,None]
  for c in range(3):f[:,c]=np.bincount(pairs[:,0],weights=ff[:,c],minlength=self.n)-np.bincount(pairs[:,1],weights=ff[:,c],minlength=self.n)
  return float(.5*self.c*np.dot(dr,dr)),f,len(pairs)
