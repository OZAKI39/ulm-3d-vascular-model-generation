"""Native VV with the accepted physical-N external force callback."""
from common import *
from shear_loading import ShearLoad
import traceback
import ast,ctypes,types
from mpi4py import MPI
from lammps import lammps
def make_case(case):
 case=Path(case);(case/'provenance').mkdir(parents=True,exist_ok=True);(case/'keys').mkdir(exist_ok=True)
 for name in ['reference.npz','reference.native','body.data']:
  p=case/'provenance'/name
  if not p.exists():p.symlink_to(os.path.relpath(O/'provenance'/name,p.parent))
 return case
def factory(case):
 source=R/'analysis/v3_softening_pilot_v1/native.py';tree=ast.parse(source.read_text());classes=[x for x in tree.body if isinstance(x,ast.ClassDef) and x.name in ['Native','Contact']];ns=dict(globals(),O=Path(case));exec(compile(ast.Module(body=classes,type_ignores=[]),str(source),'exec'),ns);return ns['Native'],ns['Contact']
def native(case,protocol='stokes',restart=None):
 case=make_case(case);Base,Contact=factory(case)
 class Native(Base):
  def __init__(self):
   super().__init__('native',restart=restart,dt=BASE)
   d=self.ref['initial_cut'].astype(float);k=SC[self.ref['bond_types']]*d
   assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(k),self.ptr(d))==0
   self.tool=ShearLoad(self.ref,protocol);self.ext=types.SimpleNamespace(error=None,last={})
  def dynamics(self):
   self.l.commands_string(f'group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}\nfix v3forces clot external pf/callback 1 1')
   def callback(caller,step,nlocal,tags,x,f):
    f[:]=0
    try:
     ids=tags.astype(int)-1;assert nlocal==self.N
     t=float(self.l.get_thermo('time'));force=self.tool.force(t);f[:]=force[ids]
     self.ext.last=dict(external_force=force,external_time=np.array(t))
    except BaseException:self.ext.error=traceback.format_exc();f[:]=np.nan
   self.callback=callback;self.l.set_fix_external_callback('v3forces',callback);self.dynamic=True;self.l.command('run 0');return self.state()
 return Native(),Contact
def predict(model,ref,x,k,onset):
 e,_=model.geometry(x);kk=np.maximum(k,np.maximum(e/ref['r0'],0));off=SC[ref['bond_types']];d=np.zeros(len(kk));d[kk>=off]=1;part=(kk>onset)&(kk<off);d[part]=off[part]*(kk[part]-onset[part])/(kk[part]*(off[part]-onset[part]));return kk,d
