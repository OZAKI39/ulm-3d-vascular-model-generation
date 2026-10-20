"""Reuse the frozen V3 native bridge; only synthetic geometry/BCs are replaced."""
from common import *
import ast,ctypes,time,types
from mpi4py import MPI
from lammps import lammps
BASE=.25e-6;RAMP=.002;END=.0025;AMPLITUDE=240e-6
SOURCE=R/'analysis/v3_softening_pilot_v1/native.py'

def native_factory(case,axis):
 tree=ast.parse(SOURCE.read_text());classes=[x for x in tree.body if isinstance(x,ast.ClassDef) and x.name in ['Native','Contact']]
 namespace=dict(globals(),O=case)
 exec(compile(ast.Module(body=classes,type_ignores=[]),str(SOURCE),'exec'),namespace)
 Base=namespace['Native'];Contact=namespace['Contact']
 class SpecimenNative(Base):
  def __init__(self,name='native'):
   super().__init__(name,dt=BASE)
   self.axis=axis;self.top=self.ref['top_ids']-1;self.released=False
   d=self.ref['initial_cut'].astype(float);k=SC[self.ref['bond_types']]*d
   assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(k),self.ptr(d))==0
  def dynamics(self):
   top=self.top+1;free=np.setdiff1d(np.arange(1,self.n+1),top)
   self.l.command('group top id '+' '.join(map(str,top)));self.l.command('group free id '+' '.join(map(str,free)))
   disp=['NULL']*3;speed=['NULL']*3;disp[axis]='v_disp';speed[axis]='v_speed'
   self.l.commands_string(f'''variable q equal ternary(time<{RAMP:.17g},time/{RAMP:.17g},1)
variable disp equal {AMPLITUDE/2:.17g}*(1-cos(PI*v_q))
variable speed equal {AMPLITUDE*np.pi/(2*RAMP):.17g}*sin(PI*v_q)
fix integrate free nve
fix motion_top top move variable {' '.join(disp+speed)}
fix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}
''')
   self.ext=types.SimpleNamespace(error=None,last={});self.dynamic=True;self.l.command('run 0');return self.state()
  def release(self):
   before=self.state();self.l.commands_string('unfix motion_top\nunfix integrate\ngroup clot type 1\nfix integrate clot nve\nrun 0');self.released=True;after=self.state()
   for k in ['x','velocity','d','kappa']:assert np.array_equal(before[k],after[k]),k
   assert np.max(abs(before['force']-after['force']))<1e-15
 return SpecimenNative,Contact

def boundary(t):
 q=min(t/RAMP,1);return AMPLITUDE/2*(1-np.cos(np.pi*q)),AMPLITUDE*np.pi/(2*RAMP)*np.sin(np.pi*q)

def update(model,ref,x,kappa):
 e,_=model.geometry(x);k=np.maximum(kappa,np.maximum(e/ref['r0'],0));on=S0[ref['bond_types']];off=SC[ref['bond_types']];d=np.zeros(len(k));part=(k>on)&(k<off);d[k>=off]=1;d[part]=off[part]*(k[part]-on[part])/(k[part]*(off[part]-on[part]));return k,d
