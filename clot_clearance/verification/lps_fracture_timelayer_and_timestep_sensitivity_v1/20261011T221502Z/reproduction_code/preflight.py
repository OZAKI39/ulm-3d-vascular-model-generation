from common import *
import subprocess,sys,shutil,ctypes,platform,time
O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
def command(a):
 t=time.perf_counter();r=subprocess.run(a,text=True,capture_output=True)
 return dict(command=a,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr,wall_s=time.perf_counter()-t)
cmds=[['pwd'],['git','status','--short'],['git','branch','--show-current'],['git','rev-parse','HEAD'],['supervisorctl','status'],['ps','-eo','pid,comm,pcpu,pmem,etimes'],['df','-h','.'],['lscpu'],['free','-b'],['nvidia-smi','--query-gpu=name,memory.used,utilization.gpu','--format=csv'],['ldd',str(B/'liblammps.so')]]
rec={'commands':[command(a) for a in cmds],'python':sys.executable,'python_version':sys.version,'platform':platform.platform(),'cgroup':{}}
for p in ['/sys/fs/cgroup/cpu/cpu.cfs_quota_us','/sys/fs/cgroup/cpu/cpu.cfs_period_us','/sys/fs/cgroup/memory/memory.limit_in_bytes','/sys/fs/cgroup/memory/memory.usage_in_bytes']:
 rec['cgroup'][p]=Path(p).read_text().strip()
dump(O/'PREFLIGHT.json',rec)
assert sha(B/'liblammps.so')==EXPECTED
files=[B/'liblammps.so',PROD,Path('/workspace/y_clot/env.sh'),PRIOR/'BUGFIX_SOURCE_DIFF.patch',R/'analysis/material_macro_audit_v1/energy.cpp',R/'analysis/material_macro_audit_v1/energy.so',R/'pd_clot/lammps_fragmentation/macro_lab_v1.py']
files += [S/'src'/n for n in ['BPM/bond_bpm_peri.cpp','BPM/bond_bpm_peri.h','BPM/bond_bpm.cpp','fix_nve.cpp','verlet.cpp','fix_move.cpp','fix_external.cpp']]
identity={str(p):sha(p) for p in files}
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',{'status':'ON_DISK_IDENTITY_PASS_RUNTIME_IDENTITY_PENDING','files_sha256':identity,'expected_final_library_sha256':EXPECTED,'production_library_sha256':sha(PROD),'source_base_commit':'8de817dd79bfe4525d5d39246a212d833e6dee07','new_native_compute_count':0})
pd=O/'provenance';pd.mkdir(exist_ok=True)
for p in files:
 if p.suffix not in ['.so']:
  dest=pd/('native_source/'+p.name if S in p.parents else p.name);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
for p in [R/'AGENTS.md',R/'PROJECT_ACTIVE_CONTEXT.md',Path('/workspace/y_clot/README_MIGRATION_ZH.txt'),Path('/etc/vast-agents-guide.md')]:shutil.copy2(p,pd/p.name)
old=O/'historical_static';old.mkdir(exist_ok=True)
for p in (PRIOR/'final_fracture').iterdir():
 if p.is_file() and p.suffix in ['.npz','.json','.lammps','.data']:shutil.copy2(p,old/p.name)
protected={}
for base in [PRIOR,V3,R/'verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z']:
 for p in base.rglob('*'):
  if p.is_file() and not any(s in p.parts for s in ['build_owner_fix','build_fixed','source_owner_fix','source_fixed']):protected[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
for p in files:protected[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
dump(O/'HISTORICAL_PROTECTION_BEFORE.json',protected)
dump(O/'SMALL_SPECIMEN_IDENTITY.json',{'status':'HISTORICAL_STATIC_SPECIMEN_IDENTIFIED_DYNAMIC_TIMING_ABSENT','source':str(PRIOR/'final_fracture'),'copied_files':{p.name:sha(p) for p in old.iterdir()},'nodes':341,'initial_bonds':12047,'historical_MPI_ranks':1,'historical_integrator':None,'historical_loading':'Direct coordinate overwrite: affine T0, then +0.3 dx at corner atom 90; no continuous physical loading function','historical_initial_velocity':'Default zero; never integrated','historical_contact':'pair zero','historical_damping':'none','baseline_reuse_limitation':'No actual dynamic baseline exists; velocities, impulse and work cannot be recovered from static force evaluations','dynamic_adaptation_authorization':'PENDING_USER_REPLY'})
print(json.dumps({'identity':identity,'protected_files':len(protected),'protected_bytes':sum(x['bytes'] for x in protected.values()),'cgroup':rec['cgroup']}))
