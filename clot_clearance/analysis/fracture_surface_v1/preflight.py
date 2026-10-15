from common import *
import subprocess,os,sys,importlib.metadata,ctypes,time,shutil
start=time.perf_counter();assert not (O/'SOURCE_AND_LIBRARY_IDENTITY.json').exists();O.mkdir(parents=True,exist_ok=True)
commands={}
for name,cmd in [('pwd',['pwd']),('git_status',['git','status','--short']),('git_branch',['git','branch','--show-current']),('git_HEAD',['git','rev-parse','HEAD']),('CPU',['lscpu']),('RAM',['free','-b']),('disk',['df','-B1','/workspace/y_clot']),('GPU',['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu','--format=csv']),('processes',['ps','-eo','pid,ppid,etimes,pcpu,pmem,comm']),('library_deps',['ldd',str(B/'liblammps.so')])]:
 p=subprocess.run(cmd,cwd=R,text=True,capture_output=True);commands[name]={'argv':cmd,'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
assert sha(B/'liblammps.so')==EXPECTED
lib=ctypes.CDLL(str(B/'liblammps.so'));paths=sorted({l.split()[-1] for l in Path('/proc/self/maps').read_text().splitlines() if 'liblammps.so' in l});assert len(paths)==1 and sha(paths[0])==EXPECTED
q={p:Path(p).read_text() for p in ['/sys/fs/cgroup/cpu.max','/sys/fs/cgroup/memory.max','/sys/fs/cgroup/cpuset.cpus.effective'] if Path(p).exists()}
dump(O/'PREFLIGHT.json',{'commands':commands,'python':sys.executable,'version':sys.version,'packages':{n:importlib.metadata.version(n) for n in ['numpy','scipy','matplotlib','psutil','mpi4py']},'cgroup':q,'native_library_handle_loaded_without_LAMMPS_instance':{p:sha(p) for p in paths},'zero_LAMMPS_instances':True,'task_threads':1,'GPU_used':False})
prov=O/'provenance';prov.mkdir(exist_ok=True)
source_paths=[S/'src/BPM/bond_bpm_peri.cpp',S/'src/BPM/bond_bpm_peri.h',S/'src/BPM/bond_bpm.cpp',S/'src/BPM/bond_bpm.h',R/'analysis/dual_track_v1/core_a.py',R/'analysis/dual_track_v1/energy_y.cpp',R/'analysis/dual_track_v1/common.py',R/'analysis/dual_track_v1/preflight_a.py',R/'analysis/material_macro_audit_v1/energy.cpp',R/'pd_clot/lammps_fragmentation/macro_lab_v1.py',R/'pd_clot/lammps_fragmentation/lammps_writer.py',PRIOR/'BUGFIX_SOURCE_DIFF.patch',Path('/workspace/y_clot/github_sync_20261010/isolated_checkout/clot_clearance/github_sync_20261010/native_patches/lammps/SOURCE.json'),R/'AGENTS.md',R/'PROJECT_ACTIVE_CONTEXT.md',Path('/workspace/y_clot/README_MIGRATION_ZH.txt'),Path('/etc/vast-agents-guide.md')]
for p in source_paths:assert p.is_file(),str(p);shutil.copy2(p,prov/p.name)
ident={str(p):{'sha256':sha(p),'bytes':p.stat().st_size} for p in source_paths+[B/'liblammps.so',PROD,R/'analysis/dual_track_v1/energy_y.so',R/'analysis/material_macro_audit_v1/energy.so',Path('/workspace/y_clot/env.sh')]};dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',{'files':ident,'loaded_library':paths,'expected':EXPECTED,'source_upstream':'8de817dd79bfe4525d5d39246a212d833e6dee07','native_CXX_changed':False})
# Protect inherited accepted evidence plus all directly used PVCE data and the latest static task.
old=json.loads((R/C['prior_static']/'HISTORICAL_PROTECTION_BEFORE.json').read_text());bad=[p for p,v in old.items() if not Path(p).is_file() or sha(p)!=v['sha256']];assert not bad,bad
protect=old.copy()
for base in [R/C['prior_static'],R/'results'/Path(C['prior_static']).parts[1]/Path(C['prior_static']).parts[2],W,R/'verification/lps_pvce_energy_release_feasibility_v1/20261011T132103Z',R/'verification/lps_pvce_continuum_and_near_onset_audit_v2/20261011T142003Z',R/'analysis/dual_track_v1',R/'analysis/material_macro_audit_v1']:
 for p in base.rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:protect[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
protect.update(ident);dump(O/'HISTORICAL_PROTECTION_BEFORE.json',protect)
# Source data copied exactly; existing pristine array includes alive masks but bonds itself is complete.
inputs=O/'source_inputs';inputs.mkdir(exist_ok=True)
for name in ['reference_topology.npz','PERIODIC_TOPOLOGY_AUDIT.json','DESIGN_FREEZE.json','PVCE_DESIGN_FREEZE.json','FINAL_STATUS.json','PVCE_ENERGY_RELEASE_RESULTS.json','REVIEW_LPS_WIDE_PURE_SHEAR_PVCE_FEASIBILITY_V1.md']:
 p=W/name;shutil.copy2(p,inputs/name)
for name in ['reference.npz','periodic_nonaffine.npz']:
 p=W/'preflight/full'/name;shutil.copy2(p,inputs/('wide_full_'+name))
rt=np.load(W/'reference_topology.npz');meta={k:{'shape':rt[k].shape,'dtype':str(rt[k].dtype)} for k in rt.files};dump(O/'SOURCE_ARRAY_INVENTORY.json',{'path':str(W/'reference_topology.npz'),'sha256':sha(W/'reference_topology.npz'),'fields':meta,'axis_counts':[len(np.unique(rt['X'][:,k])) for k in range(3)],'axes_min':rt['X'].min(0),'axes_max':rt['X'].max(0),'preflight_s':time.perf_counter()-start,'protected_files':len(protect),'protected_bytes':sum(v['bytes'] for v in protect.values())});print(json.dumps({'reference_topology':meta,'protected_files':len(protect),'seconds':time.perf_counter()-start},indent=2))
