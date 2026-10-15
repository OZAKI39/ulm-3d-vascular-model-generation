"""Two authorized tiny dynamics only. Native VV/force semantics unchanged.
Same historical X initializes reference m; exact T0 initializes dynamic x.
Only physics adaptation: prescribed corner motion and native NVE integration.
"""
from common import *
import os,sys,time,resource,traceback,ctypes as ct,numpy as np,csv
from model import LPS
from mpi4py import MPI
from lammps import lammps
import lammps as lammps_module
sys.path.insert(0,str(R))
from pd_clot.lammps_fragmentation.lammps_writer import force_conversion
name=sys.argv[1];case=json.loads((O/name/'case.json').read_text());out=O/name;dt=case['dt_s'];steps=case['steps'];freeze=json.loads((O/'TIMESTEP_DESIGN_FREEZE.json').read_text())
assert MPI.COMM_WORLD.size==1 and sha(B/'liblammps.so')==EXPECTED
for n,h in freeze['diagnostic_source_sha256'].items():assert sha(Path(__file__).parent/n)==h,(n,'executed code changed after freeze')
with (out/'LAUNCHED.json').open('x') as f:json.dump({'epoch':time.time(),'pid':os.getpid(),'case':case},f)
os.chdir(out);g=np.load(O/'historical_static/reference.npz');model=LPS(g);n=len(g['X']);nb=len(g['bonds']);V=model.V;mass=1000*V;free=np.ones(n,bool);free[89]=False
key=g['bonds'][:,0]*n+g['bonds'][:,1];alive=np.ones(nb,bool);snapshots=[];halves={};events=[];checks=[];last=-1;armed=False;beg=time.perf_counter()
L=lammps(comm=MPI.COMM_WORLD,cmdargs=['-log',str(out/'log.lammps'),'-screen',str(out/'screen.lammps')])
paths=sorted({s.split()[-1] for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps.so' in s or 'libmpi.so' in s});libpaths=[p for p in paths if 'liblammps.so' in p]
assert len(libpaths)==1 and sha(libpaths[0])==EXPECTED and Path(libpaths[0]).parent==B
dump(out/'LOADED_LIBRARY.json',{'actual_maps':{p:sha(p) for p in paths},'python':sys.executable,'python_version':sys.version,'lammps_python':lammps_module.__file__,'lammps_version':L.version(),'MPI':MPI.Get_library_version(),'diagnostic_library_sha256':EXPECTED})
lib=ct.CDLL(str(O/'observer.so'));lib.observe.argtypes=[ct.c_void_p,ct.c_int];lib.observe.restype=ct.POINTER(ct.c_double);CB=ct.CFUNCTYPE(None,ct.c_longlong,ct.c_int);lib.register_observer.argtypes=[CB]
header=(O/'historical_static/in.lammps').read_text().replace(str(PRIOR/'final_fracture/body.data'),str(out/'body.data')).replace('timestep 5e-7',f'timestep {dt:.17g}')
(out/'reference_initialization.in').write_text(header)
def fail(stage):
 dump(out/'CALLBACK_ERROR.json',{'stage':stage,'last_accepted':last,'traceback':traceback.format_exc()});os._exit(71)
def post_force(caller,step,nlocal,tags,x,fext):
 try:
  assert nlocal==n and np.array_equal(tags,np.arange(1,n+1));fext[:]=0
  if armed:halves[int(step)]=L.numpy.extract_atom('v')[:n].copy()
 except BaseException:fail('post_force')
stream=(out/'native_scalars.csv').open('w',buffering=1);w=csv.writer(stream);w.writerow(['step','time_s','alive_bonds','new_breaks','vmax_all_m_s','vmax_free_m_s','max_step_displacement_m','K_free_J','K_all_J','max_m_relative_error','native_vs_mixed_max_abs_N','integrator_velocity_error_m_s'])
def accepted(step,kind):
 global alive,last
 try:
  if not armed:return
  step=int(step);assert step==last+1,(step,last,kind)
  ids=L.numpy.extract_atom('id')[:n].copy();assert np.array_equal(ids,np.arange(1,n+1))
  x=L.numpy.extract_atom('x')[:n].copy();v=L.numpy.extract_atom('v')[:n].copy();f=L.numpy.extract_atom('f')[:n].copy()
  m=np.ctypeslib.as_array(lib.observe(L.lmp,0),shape=(n,)).copy();th=np.ctypeslib.as_array(lib.observe(L.lmp,1),shape=(n,)).copy()
  count,raw=L.gather_bonds();rec=np.ctypeslib.as_array(raw,shape=(3*count,)).reshape(-1,3);ends=np.sort(rec[:,1:].astype(np.int64)-1,axis=1);idx=np.searchsorted(key,ends[:,0]*n+ends[:,1]);assert np.all(idx<nb) and np.array_equal(key[idx],ends[:,0]*n+ends[:,1]) and len(np.unique(idx))==count
  now=np.zeros(nb,bool);now[idx]=True;assert np.all(now<=alive) and int(L.extract_global('nbonds'))==int(now.sum())
  broken=np.flatnonzero(alive&~now);tp,_,_=model.eval(x,alive);_,fmix,_=model.eval(x,now,tp)
  ferr=float(np.linalg.norm(f-fmix,axis=1).max());merr=float(np.max(abs(m-g['m'])/g['m']))
  assert merr<=freeze['gates']['m_relative_per_ID'] and ferr<=freeze['gates']['independent_force_relative_to_historical_scale']*freeze['reference_scales']['force_N'],(merr,ferr)
  assert max(abs(th-tp))<1e-12
  half=halves.pop(step);verr=0.;dxmax=0.
  if step:
   prev=snapshots[-1];dx=x-prev['x'];dxmax=float(np.linalg.norm(dx,axis=1).max());prediction=half.copy();prediction[free]+=.5*dt*f[free]/mass;verr=float(np.max(abs(prediction[free]-v[free])))
   assert verr<freeze['gates']['kinematic_velocity_abs_m_s']
   assert np.max(abs(half[free]-prev['v_full'][free]-.5*dt*prev['force_N'][free]/mass))<2e-10
   assert np.max(abs(dx[free]-dt*half[free]))<2e-16
  vmax=float(np.linalg.norm(v,axis=1).max());assert np.isfinite(x).all() and np.isfinite(v).all() and np.isfinite(f).all() and vmax<20 and dxmax<.2*50e-6
  events.extend([step,float(step*dt),int(j+1),int(g['bonds'][j,0]+1),int(g['bonds'][j,1]+1)] for j in broken)
  snapshots.append(dict(step=step,time_s=step*dt,x=x,v_full=v,v_half=half,force_N=f,theta_native=th,alive_packed=np.packbits(now),m=m))
  checks.append(dict(step=step,m_relative=merr,force_error_N=ferr,velocity_error_m_s=verr))
  w.writerow([step,step*dt,int(now.sum()),len(broken),vmax,float(np.linalg.norm(v[free],axis=1).max()),dxmax,.5*mass*float(np.sum(v[free]**2)),.5*mass*float(np.sum(v**2)),merr,ferr,verr])
  if len(broken) or step in [0,steps]:
   np.savez_compressed(out/f'key_{step:04d}.npz',**snapshots[-1],alive=now,alive_before=alive)
  alive=now;last=step
 except BaseException:fail('end_of_step')
callback=CB(accepted);lib.register_observer(callback)
try:
 L.commands_string(header)
 assert np.array_equal(L.numpy.extract_atom('id')[:n],np.arange(1,n+1))
 # Set the approved prestrained initial coordinates only after reference history initialized.
 L.numpy.extract_atom('x')[:n]=np.load(out/'initial_T0.npz')['x']
 tr=freeze['loading']['ramp_s'];amp=freeze['loading']['amplitude_m']
 text=force_conversion()+f'''velocity all set 0 0 0 units box
group driven id 90
group free subtract all driven
variable zero equal 0
variable disp equal {amp:.17g}*0.5*(1-cos(PI*min(time/{tr:.17g},1)))
variable speed equal (time<{tr:.17g})*{amp:.17g}*PI/(2*{tr:.17g})*sin(PI*min(time/{tr:.17g},1))
fix drive driven move variable v_disp v_zero v_zero v_speed v_zero v_zero units box
fix integrate free nve
fix half_read all external pf/callback 1 1
fix full_read all audit/readonly/finalv3
thermo_style custom step time atoms bonds ke pe
thermo_modify flush yes
'''
 (out/'dynamics.in').write_text(text+f'run {steps}\n');L.commands_string(text);L.set_fix_external_callback('half_read',post_force);armed=True;L.command(f'run {steps}');assert last==steps and len(snapshots)==steps+1
 np.savez_compressed(out/'native_states.npz',**{k:np.asarray([s[k] for s in snapshots]) for k in snapshots[0]},bond_count=nb,free_node_mask=free)
 np.savez_compressed(out/'native_break_events.npz',events=np.asarray(events).reshape(-1,5))
 with (out/'native_break_events.csv').open('w') as ef:
  ew=csv.writer(ef);ew.writerow(['step','time_s','bond_id','atom_i','atom_j']);ew.writerows(events)
 dump(out/'RUNTIME_CHECKS.json',{'status':'PASS','accepted_states':len(snapshots),'max_m_relative':max(q['m_relative'] for q in checks),'max_native_vs_mixed_force_N':max(q['force_error_N'] for q in checks),'max_velocity_reconstruction_error_m_s':max(q['velocity_error_m_s'] for q in checks),'events':len(events),'actual_velocity_layer':'Direct END_OF_STEP after native final_integrate','half_velocity_layer':'Post_force prior to final_integrate, not reconstructed','physics_changes_relative_to_static':'Only explicitly approved kinematic adaptation/native integration and existing required N/m3 to N conversion','extra_force_from_audit_N':0,'checks':checks})
 L.close();q=resource.getrusage(resource.RUSAGE_SELF);dump(out/'PERFORMANCE.json',{'status':'COMPLETED','accepted_step':last,'time_s':last*dt,'wall_s':time.perf_counter()-beg,'CPU_user_s':q.ru_utime,'CPU_system_s':q.ru_stime,'max_RSS_KiB':q.ru_maxrss,'MPI_ranks':1,'GPU_used':False,'native_instances':1,'reference_initialization_run0':1,'formal_dynamic_runs':1})
 print(name,'COMPLETED',last,'breaks',len(events),flush=True)
except BaseException:
 dump(out/'DRIVER_ERROR.json',{'last_accepted':last,'traceback':traceback.format_exc()});raise
finally:stream.close()
