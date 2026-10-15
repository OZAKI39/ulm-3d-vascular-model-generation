"""Task-isolated static native evaluator. No time integrator and no production changes."""
from common import *
import numpy as np,ctypes as ct,traceback
from mpi4py import MPI
from lammps import lammps
from pd_clot.lammps_fragmentation.lammps_writer import force_conversion
from pd_clot.lammps_fragmentation.native_y_wall_contact import WallCallback

class Native:
 def __init__(self,out,comm=MPI.COMM_WORLD):
  self.comm=comm;self.rank=comm.rank;self.out=out;out.mkdir(parents=True,exist_ok=True);self.ref=np.load(O/'REFERENCE_AND_LOAD.npz');self.X=self.ref['X'];self.load=self.ref['force_N'];self.enabled=False;self.calls=0;self.eval_index=0
  self.L=lammps(comm=comm,cmdargs=['-log',str(out/'log.lammps'),'-screen','none']);dump(out/f'LOADED_LIBRARY_rank{self.rank}.json',loaded(self.L))
  self.probe=ct.CDLL(str(O/'observer.so'));self.probe.observe.argtypes=[ct.c_void_p,ct.c_int];self.probe.observe.restype=ct.POINTER(ct.c_double);self.probe.counter.argtypes=[ct.c_void_p,ct.c_int];self.probe.counter.restype=ct.c_longlong;self.probe.permutation.argtypes=[ct.c_void_p];self.probe.permutation.restype=ct.POINTER(ct.c_int);self.probe.endpoint_records.argtypes=[ct.c_void_p,ct.c_void_p,ct.c_int];self.probe.endpoint_records.restype=ct.c_int;self.probe.contact_readout.argtypes=[ct.c_void_p,ct.c_void_p,ct.c_void_p,ct.c_int];self.probe.contact_readout.restype=ct.c_int
  original=(O/'provenance/original_C1.in.lammps').read_text().splitlines();elastic='\n'.join(l for l in original if l.startswith(('pair_style ','pair_coeff ','bond_coeff ')))
  deck=f'''# Static evaluator; original V3 geometry and elastic/contact/critical parameters.
units si
dimension 3
boundary m m m
atom_style bond
processors {comm.size} 1 1
atom_modify map array sort 1 {DX:.17g}
newton on off
special_bonds lj 0 1 1 coul 1 1 1
lattice sc {DX:.17g}
read_data {O/'source_inputs/body.data'} extra/bond/per/atom 160 extra/special/per/atom 180
neighbor {2*DX:.17g} bin
neigh_modify every 1 delay 0 check no
comm_modify cutoff {3*DELTA:.17g}
fix props all property/atom d_vfrac d_theta d_damage ghost yes
set group all d_vfrac {V:.17g}
group clot type 1
group anchors type 2
bond_style bpm/peri
{elastic}
velocity all set 0 0 0 units box
{force_conversion()}
fix wallprops all property/atom d_wclear d_wpen d_wtri ghost yes
fix nativewall clot external pf/callback 1 1
fix cfdload all external pf/callback 1 1
timestep 5e-7
compute contactPE all pe pair
compute bondPE all pe bond
thermo 1
thermo_style custom step atoms bonds ke c_contactPE c_bondPE
thermo_modify lost error lost/bond error
# No nve, no damping, no streaming or drag; anchor coordinates constrained by static optimizer.
# No anchor setforce: physical internal force is retained for reaction recovery.
'''
  if self.rank==0:(out/'in.lammps').write_text(deck+'# Python callbacks are registered before run 0. Subsequent run 1 pre no post no are force evaluations, not dynamic steps.\n')
  self.L.commands_string(deck);self.wall=WallCallback(self.L,comm,out,NC);self.L.set_fix_external_callback('nativewall',self.wall.callback);self.L.set_fix_external_callback('cfdload',self.callback);self.L.command('run 0')
 def callback(self,caller,step,nlocal,tags,x,fext):
  try:
   ids=np.asarray(tags[:nlocal],dtype=int);assert np.all((ids>=1)&(ids<=N));fext[:]=0
   self.before=self.L.numpy.extract_atom('f')[:nlocal].copy();self.callback_ids=ids.copy();self.applied=self.load[ids-1].copy() if self.enabled else np.zeros((nlocal,3));fext[:nlocal]=self.applied;self.calls+=1
  except BaseException:
   (self.out/f'CALLBACK_ERROR_rank{self.rank}.txt').write_text(traceback.format_exc());self.comm.Abort(47)
 def evaluate(self,x,label,enabled=True,initial=False):
  L=self.L;self.enabled=bool(enabled);b=self.ref['bonds'];r=self.ref['reference_lengths_m'];s=(np.linalg.norm(x[b[:,1]]-x[b[:,0]],axis=1)-r)/r;crit=np.choose(self.ref['bond_types']-1,[.15,.140624859375,1.]);assert np.max(s/crit)<1,'SUBCRITICAL_GUARD'
  assert np.array_equal(x[NC:],self.X[NC:]),'ANCHOR_MOVEMENT'
  if not initial:
   xx=np.ascontiguousarray(x,dtype=np.float64);L.scatter_atoms('x',1,3,xx.ctypes.data_as(ct.POINTER(ct.c_double)));L.command('run 1 pre no post no');self.eval_index+=1
  nl=L.extract_global('nlocal');ng=L.extract_global('nghost');nt=nl+ng;allids=L.numpy.extract_atom('id')[:nt].astype(int).copy();ids=allids[:nl];assert np.array_equal(ids,self.callback_ids)
  mm=np.ctypeslib.as_array(self.probe.observe(L.lmp,0),shape=(nt,)).copy();theta=np.ctypeslib.as_array(self.probe.observe(L.lmp,1),shape=(nt,)).copy();refm=self.ref['reference_m'][allids-1];merr=float(np.max(abs(mm-refm)/refm));assert merr<2e-12
  f=L.numpy.extract_atom('f')[:nl].copy();actual_x=L.numpy.extract_atom('x')[:nl].copy();v=L.numpy.extract_atom('v')[:nl].copy();assert np.array_equal(actual_x,x[ids-1]);assert not np.any(v)
  endpoints=np.empty((2*len(b),3),np.int32);ne=self.probe.endpoint_records(L.lmp,endpoints.ctypes.data,len(endpoints));assert ne>=0;endpoints=endpoints[:ne].copy()
  pairs=np.sort(endpoints[:,:2].astype(np.int64)-1,axis=1);keys=b[:,0]*N+b[:,1];ki=np.searchsorted(keys,pairs[:,0]*N+pairs[:,1]);assert np.all(ki<len(b));assert np.array_equal(keys[ki],pairs[:,0]*N+pairs[:,1]);assert np.array_equal(self.ref['bond_types'][ki],endpoints[:,2])
  deg=np.bincount(b.ravel(),minlength=N);ndeg=np.bincount(endpoints[:,0]-1,minlength=N);assert np.array_equal(ndeg[ids-1],deg[ids-1]);assert len(np.unique((endpoints[:,0].astype(np.int64)-1)*N+endpoints[:,1]-1))==ne
  cf=np.zeros((nt,3));ce=np.empty((len(b),2),np.int32);npair=self.probe.contact_readout(L.lmp,cf.ctypes.data,ce.ctypes.data,len(ce));assert npair>=0;cglobal=np.zeros((N,3));np.add.at(cglobal,allids-1,cf*V)
  nbin=int(self.probe.counter(L.lmp,3));perm=np.ctypeslib.as_array(self.probe.permutation(L.lmp),shape=(nl,)).copy() if nbin>1 else np.arange(nl)
  packets=self.comm.allgather({'ids':ids,'x':actual_x,'force':f,'before':self.before,'applied':self.applied,'m':mm[:nl],'theta':theta[:nl],'merr':merr,'endpoints':ne,'contact_force':cglobal,'contact_edges':ce[:npair].copy(),'nonidentity_permutation':int(np.count_nonzero(perm!=np.arange(nl))),'ghost_count':ng})
  assert np.array_equal(np.sort(np.concatenate([p['ids'] for p in packets])),np.arange(1,N+1));assert sum(p['endpoints'] for p in packets)==2*len(b)
  arrays={k:np.zeros((N,3)) for k in ['x','force','before','applied']};arrays.update(m=np.zeros(N),theta=np.zeros(N),owner=np.zeros(N,int));contact=np.zeros((N,3))
  for rank,p in enumerate(packets):
   for k in ['x','force','before','applied','m','theta']:arrays[k][p['ids']-1]=p[k]
   arrays['owner'][p['ids']-1]=rank;contact+=p['contact_force']
  arrays['contact_force_N']=contact;arrays['wall_force_N']=np.zeros((N,3));arrays['wall_clearance_m']=np.zeros(NC)
  wr=self.comm.allgather((ids[ids<=NC],self.wall.clearance[ids[ids<=NC]-1].copy(),self.wall.normal[ids[ids<=NC]-1].copy()))
  from pd_clot.lammps_fragmentation.native_y_wall_contact import penalty
  for ii,clear,normal in wr:arrays['wall_clearance_m'][ii-1]=clear;arrays['wall_force_N'][ii-1]=penalty(clear,normal,np.zeros((len(ii),3)))
  expected=self.load if enabled else np.zeros_like(self.load);derr=float(np.linalg.norm(arrays['force']-arrays['before']-expected,axis=1).max());assert derr<1e-18
  diag={'label':label,'native_force_evaluation_counter':int(L.extract_global('ntimestep')),'physical_dynamic_steps':0,'CFD_enabled':bool(enabled),'MPI_ranks':self.comm.size,'owned_counts':[len(p['ids']) for p in packets],'ghost_counts':[p['ghost_count'] for p in packets],'nonidentity_permutation_counts':[p['nonidentity_permutation'] for p in packets],'m_max_relative':max(p['merr'] for p in packets),'post_force_increment_error_N':derr,'bonds':len(b),'new_breaks':0,'contact_pairs':sum(len(p['contact_edges']) for p in packets),'pair_contact_PE_J':float(L.extract_compute('contactPE',0,0))*V,'native_bond_PE_times_V_J':float(L.extract_compute('bondPE',0,0))*V,'wall_active_nodes':int(np.sum(np.linalg.norm(arrays['wall_force_N'],axis=1)>0)),'wall_min_clearance_m':float(arrays['wall_clearance_m'].min()),'KE_J':float(L.get_thermo('ke'))}
  if self.rank==0:np.savez_compressed(self.out/(label+'.npz'),**arrays);dump(self.out/(label+'.json'),diag)
  return arrays,diag
 def close(self):self.wall.close();self.L.close()
