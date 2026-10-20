"""Native VV and native fix indent; Python independently audits, never applies force."""
from common import *
from indenter import Indenter
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
def native(case,heterogeneous=False,restart=None):
 case=make_case(case);Base,Contact=factory(case)
 class Native(Base):
  def __init__(self):
   super().__init__('native',restart=restart,dt=BASE)
   self.lib.lpsvar_set_onset.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.POINTER(ctypes.c_longlong),ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_set_onset.restype=ctypes.c_int
   if not restart:
    if heterogeneous:
     field=load_npz(O/'HETEROGENEITY_FIELD.npz');ids=field['original_bond_ids'].copy();on=field['s0_by_original_bond_id'].copy()
     # Reverse input ordering deliberately; native identity is original bond ID.
     ids=ids[::-1].copy();on=on[::-1].copy();assert self.lib.lpsvar_set_onset(self.l.lmp,len(ids),ids.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),self.ptr(on))==0
    d=self.ref['initial_cut'].astype(float);k=SC[self.ref['bond_types']]*d;assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(k),self.ptr(d))==0
   self.tool=Indenter(self.ref);self.ext=types.SimpleNamespace(error=None,last={})
  def dynamics(self):
   self.l.commands_string(f'group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}\n'+self.tool.command());self.dynamic=True;self.l.command('run 0');return self.state()
  def state(self):
   s=super().state()
   if self.dynamic:
    contact=self.tool.evaluate(s['x'],s['time']);reaction=np.array([self.l.extract_fix('indenter',0,1,nrow=i) for i in range(3)])
    s.update(indenter_force=contact['force'],indenter_U=np.array(self.l.extract_fix('indenter',0,0)),indenter_U_oracle=np.array(contact['U']),indenter_reaction=reaction,indenter_center=contact['center'],contact_atom_ids=contact['contact_atom_ids'],indenter_overlap=contact['overlap'])
   return s
 return Native(),Contact
def predict(model,ref,x,k,onset):
 e,_=model.geometry(x);kk=np.maximum(k,np.maximum(e/ref['r0'],0));off=SC[ref['bond_types']];d=np.zeros(len(kk));d[kk>=off]=1;part=(kk>onset)&(kk<off);d[part]=off[part]*(kk[part]-onset[part])/(kk[part]*(off[part]-onset[part]));return kk,d
