"""Exactly one V3 replay, unchanged input deck/callback physics; scoped outputs."""
from pathlib import Path
import os,json,hashlib,shutil,time,resource,ctypes as ct
import numpy as np
from mpi4py import MPI
from lammps import lammps
from pd_clot.lammps_fragmentation import early_flow_coupling as E
C=json.loads((Path(__file__).parent/'config.json').read_text());R=Path(C['server']);O=R/'verification'/C['task']/C['run'];P=O/'V3_replay';OLD=R/'verification/y_clot_point_cloud_and_early_flow_coupling_v1/cases/C1_GAUSSIAN_A200';comm=MPI.COMM_WORLD
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
gates=json.loads((O/'SMALL_TEST_GATES.json').read_text());assert gates['V3_eligible']
if comm.rank==0:
 P.mkdir(exist_ok=False)
 for name in ['body.data','request.json','LOAD_PROVIDER.json','in.lammps']:shutil.copyfile(OLD/name,P/name)
 identity={'input_exact':{name:sha(P/name)==sha(OLD/name) for name in ['body.data','request.json','LOAD_PROVIDER.json','in.lammps']},'input_sha256':{name:sha(P/name) for name in ['body.data','request.json','LOAD_PROVIDER.json','in.lammps']},'non_target_python_sha256':{name:sha(R/'pd_clot/lammps_fragmentation'/name) for name in ['early_flow_coupling.py','native_y_wall_contact.py','native_y_wall_query.py','y_fragmentation_cases.py','current_position_provider.py']},'library_sha256':sha(O/'build/liblammps.so'),'baseline_path':str(OLD),'original_baseline_rerun':False,'max_steps':4660,'mpi_ranks':2,'instrumentation':'read-only state identity observer after existing wall callback; output root E.OUT redirected to new run; full exact deck unchanged','sort_frequency':'unchanged LAMMPS default1000','force_energy_gate':'independent existing C++ potential gradient PASS; native PE tally caveat explicitly retained'}
 (P/'REPLAY_DESIGN_AND_IDENTITY.json').write_text(json.dumps(identity,indent=2));(O/'V3_LAUNCHED.json').write_text(json.dumps({'one_formal_replay':True,'epoch':time.time()}))
comm.Barrier();os.chdir(P);E.OUT=P # EarlyFlow.close writes one shared audit here, not accepted historical OUT.
g=E.load_geometry();X=g['X'];b=g['bonds'];r=np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1);h=np.where(abs(r-150.5e-6)<=25e-6,(175.5e-6-r)/50e-6,1);m=np.bincount(b.ravel(),weights=np.repeat(r*h*(50e-6)**3,2),minlength=len(X));probe=ct.CDLL(str(O/'probe_fixed.so'));probe.observe.argtypes=[ct.c_void_p,ct.c_int];probe.observe.restype=ct.POINTER(ct.c_double)
class ObserverWall(E.AuditedWall):
 def __init__(self,*args):
  super().__init__(*args);self.last_ids=None;self.identity_rows=[];self.reorders=[]
 def compute(self,step,nlocal,tags,x,f):
  super().compute(step,nlocal,tags,x,f)
  ids=tags.astype(int);mm=np.ctypeslib.as_array(probe.observe(self.lmp.lmp,0),shape=(nlocal,)).copy();err=float(np.max(abs(mm-m[ids-1]))/m.max())
  changed=self.last_ids is not None and not np.array_equal(ids,self.last_ids)
  oldset=set() if self.last_ids is None else set(self.last_ids.tolist());incoming=[] if self.last_ids is None else sorted(set(ids.tolist())-oldset)
  if changed:self.reorders.append({'step':int(step),'rank':comm.rank,'before_ids':self.last_ids.tolist(),'after_ids':ids.tolist(),'incoming_ids':incoming})
  reaction=self.lmp.numpy.extract_atom('f')[:nlocal][ids>1927].sum(0)
  self.identity_rows.append([int(step),nlocal,float(err),int(changed),len(incoming),*reaction.tolist()]);self.last_ids=ids.copy()
 def close(self):
  super().close();np.savez_compressed(P/f'identity_rank{comm.rank}.npz',history=np.asarray(self.identity_rows),columns=['step','nlocal','m_relerror','local_order_changed','incoming_ids','anchor_fx_N','anchor_fy_N','anchor_fz_N']);(P/f'identity_events_rank{comm.rank}.json').write_text(json.dumps(self.reorders))
t0=time.perf_counter();L=lammps(comm=comm,cmdargs=['-log','log.lammps','-screen','screen.lammps']);wall=flow=None
try:
 if comm.rank==0:
  (P/'LOADED_LIBRARY_MAPS.txt').write_text('\n'.join(s for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in s or 'libmpi.so' in s))
 prefix,suffix=(P/'in.lammps').read_text().split('\nrun ',1);L.commands_string(prefix);wall=ObserverWall(L,comm,P,1927);flow=E.EarlyFlow(L,comm,P,wall)
 L.set_fix_external_callback('nativewall',wall.callback);L.set_fix_external_callback('earlyflow',flow.callback);L.commands_string('run '+suffix)
 xx=np.ctypeslib.as_array(L.gather_atoms('x',1,3)).reshape(-1,3).copy();vv=np.ctypeslib.as_array(L.gather_atoms('v',1,3)).reshape(-1,3).copy()
 if comm.rank==0:np.savez_compressed(P/'terminal_full_precision.npz',position_m=xx,velocity_m_s=vv)
finally:
 if flow:flow.close()
 if wall:wall.close()
 L.close();q=resource.getrusage(resource.RUSAGE_SELF);stats=comm.gather({'rank':comm.rank,'wall_s':time.perf_counter()-t0,'cpu_s':q.ru_utime+q.ru_stime,'max_RSS_KiB':q.ru_maxrss},root=0)
 if comm.rank==0:(P/'PERFORMANCE.json').write_text(json.dumps(stats,indent=2))
