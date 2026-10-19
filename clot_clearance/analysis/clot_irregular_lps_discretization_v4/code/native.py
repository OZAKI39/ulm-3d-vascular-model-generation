"""Isolated V4 native driver: nonuniform masses, SI forces, unchanged commit order."""
from common import *
import ast,ctypes,types,traceback
from mpi4py import MPI
from lammps import lammps

def make_case(case):
 case=Path(case);(case/'provenance').mkdir(parents=True,exist_ok=True);(case/'keys').mkdir(exist_ok=True)
 for name in ['reference.npz','reference.native','body.data']:
  p=case/'provenance'/name
  if not p.exists():p.symlink_to(os.path.relpath(O/'provenance'/name,p.parent))
 return case

def base_factory(case):
 source=R/'analysis/v3_softening_pilot_v1/native.py';text=source.read_text()
 text=text.replace('14:self.NB}', '14:self.NB,15:self.N}')
 text=text.replace("self.gather('f')*(1 if self.dynamic else V)","self.gather('f')")
 text=text.replace("float(self.l.get_thermo('epair'))*V","float(self.l.get_thermo('epair'))")
 text=text.replace("float(self.l.get_thermo('ebond'))*V","float(self.l.get_thermo('ebond'))")
 marker='  if not restart:\n   self.l.commands_string'
 assert text.count(marker)==1
 text=text.replace(marker,"""  self.l.command('fix nodal_mass all property/atom rmass ghost yes')
  ids=self.l.numpy.extract_atom('id')[:self.N].astype(int)-1
  self.l.numpy.extract_atom('rmass')[:self.N]=RHO*self.ref['V'][ids]
  self.l.numpy.extract_atom('d_vfrac')[:self.N]=self.ref['V'][ids]
"""+marker)
 tree=ast.parse(text);classes=[q for q in tree.body if isinstance(q,ast.ClassDef) and q.name=='Native'];ns=dict(globals(),O=Path(case));exec(compile(ast.Module(body=classes,type_ignores=[]),str(source),'exec'),ns);return ns['Native']

class Contact:
 def __init__(self,ref):
  from scipy.spatial import cKDTree
  self.tree=cKDTree;self.b=ref['bonds'];self.N=len(ref['X']);self.codes=self.b[:,0]*self.N+self.b[:,1];self.vol=ref['V'];self.c=15*(18*K/(np.pi*HORIZON**4))/HORIZON
 def coefficients(self,pairs):return self.c*((self.vol[pairs[:,0]]+self.vol[pairs[:,1]])/2)**2
 def energy_for_pairs(self,x,pairs):
  return float(.5*np.sum(self.coefficients(pairs)*np.minimum(np.linalg.norm(x[pairs[:,1]]-x[pairs[:,0]],axis=1)-DX,0)**2))
 def evaluate(self,x,d):
  pairs=self.tree(x).query_pairs(DX,output_type='ndarray');pairs=pairs[~np.isin(pairs[:,0]*self.N+pairs[:,1],self.codes[d<1])];f=np.zeros_like(x)
  if not len(pairs):return 0.,f,pairs
  delta=x[pairs[:,1]]-x[pairs[:,0]];r=np.linalg.norm(delta,axis=1);dr=r-DX;c=self.coefficients(pairs);ff=c[:,None]*dr[:,None]*delta/r[:,None]
  for k in range(3):f[:,k]=np.bincount(pairs[:,0],weights=ff[:,k],minlength=self.N)-np.bincount(pairs[:,1],weights=ff[:,k],minlength=self.N)
  return float(.5*np.dot(c,dr*dr)),f,pairs

def native(case,protocol='biphasic',restart=None,install_properties=True,surface=True):
 case=make_case(case);Base=base_factory(case)
 class Native(Base):
  def __init__(self):
   super().__init__('native',restart=restart,dt=BASE)
   self.vol=self.ref['V'];self.mass=RHO*self.vol
   if restart is None:
    on=S0[self.ref['bond_types']].copy();off=SC[self.ref['bond_types']].copy()
    if install_properties:
     if protocol=='biphasic':field=load_npz(O/'BIPHASIC_FIELD.npz');on=field['bond_s0'];off=field['bond_sc']
     ids=np.arange(1,self.NB+1,dtype=np.int64)[::-1].copy();self.lib.lpsvar_set_thresholds.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.POINTER(ctypes.c_longlong),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)];self.lib.lpsvar_set_thresholds.restype=ctypes.c_int
     assert self.lib.lpsvar_set_thresholds(self.l.lmp,self.NB,ids.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),self.ptr(on[::-1].copy()),self.ptr(off[::-1].copy()))==0
    d=self.ref['initial_cut'].astype(float);k=.30*d;assert self.lib.lpsvar_action(self.l.lmp,1,0,self.ptr(k),self.ptr(d))==0
   assert np.array_equal(self.field(15),self.vol)
   if surface:
    from current_surface import CurrentSurface
    self.tool=CurrentSurface(self.ref,'uniform')
   self.ext=types.SimpleNamespace(error=None,last={});self.force_stage_dt=0.
  def step(self,h):
   if h!=self.dt:self.l.command(f'timestep {h:.17g}');self.dt=h
   self.force_stage_dt=h
   try:self.l.command('run 1 pre no post no')
   finally:self.force_stage_dt=0.
   return self.state()
  def dynamics(self):
   self.l.commands_string(f'group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}\nfix v3forces clot external pf/callback 1 1')
   def callback(caller,step,nlocal,tags,x,f):
    f[:]=0
    try:
     ids=tags.astype(int)-1;assert nlocal==self.N and np.array_equal(np.sort(ids),np.arange(self.N))
     current=np.empty((self.N,3));velocity=np.empty_like(current);current[ids]=x;velocity[ids]=self.l.numpy.extract_atom('v')[:self.N]
     t=float(self.l.get_thermo('time'));d=self.field(3);drift=velocity.copy();iterations=0
     if self.force_stage_dt:
      internal=np.empty_like(current);internal[ids]=self.l.numpy.extract_atom('f')[:self.N];factor=self.force_stage_dt/(2*self.mass[:self.n,None])
      for iterations in range(1,13):
       z=self.tool.evaluate(current,velocity,d,t);updated=drift.copy();updated[:self.n]+=factor*(internal[:self.n]+z['force'][:self.n]);err=np.max(abs(updated-velocity));velocity=updated
       if err<2e-13*max(1.,float(np.max(abs(velocity)))):break
      else:raise RuntimeError('Follower slip half-kick did not converge')
     z=self.tool.evaluate(current,velocity,d,t);f[:]=z['force'][ids]
     self.ext.last=dict(external_force=z['force'],external_drive_force=z['drive_force'],external_drag_force=z['drag_force'],external_accessible_force=np.asarray(z['W']@(z['face_force']*z['accessible'][:,None])),external_L2_force=np.asarray(z['W']@(z['face_force']*(self.tool.method==1)[:,None])),external_time=np.array(t),force_stage_velocity=velocity,verlet_drift_velocity=drift,slip_iterations=np.array(iterations),surface_centroid=z['centroid'],surface_normal=z['normal'],surface_area=z['area'],surface_active=z['active'],surface_velocity=z['face_velocity'],surface_outer_velocity=z['outer_velocity'],surface_traction=z['traction'],surface_force=z['face_force'],surface_weight=z['weight'],surface_original=z['original'],surface_new_open=z['new_open'],surface_accessible=z['accessible'],surface_unresolved=z['unresolved'],surface_aperture_m=z['aperture_m'],surface_fluid_region=z['fluid_region'],surface_occluded=z['occluded'],surface_method=self.tool.method,surface_class=z['face_class'],field_amplitude_Pa=np.array(self.tool.field.amplitude(t)))
    except BaseException:self.ext.error=traceback.format_exc();f[:]=np.nan
   self.callback=callback;self.l.set_fix_external_callback('v3forces',callback);self.dynamic=True;self.l.command('run 0');return self.state()
 return Native(),Contact

def predict(model,ref,x,k,onset,off):
 e,_=model.geometry(x);kk=np.maximum(k,np.maximum(e/ref['r0'],0));d=np.zeros(len(kk));d[kk>=off]=1;part=(kk>onset)&(kk<off);d[part]=off[part]*(kk[part]-onset[part])/(kk[part]*(off[part]-onset[part]));return kk,d
