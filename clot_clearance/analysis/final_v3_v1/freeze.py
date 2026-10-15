from common import *
import subprocess,shutil,sys,os,time
O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True);P.mkdir(exist_ok=False)
prior=json.loads((PRIOR/'SOURCE_AND_BINARY_IDENTITY.json').read_text());patch=json.loads((PRIOR/'FINAL_PATCH_IDENTITY.json').read_text())
files={}
for name,q in patch.items():
 if not isinstance(q,dict):continue
 p=S/name;files[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size};assert files[str(p)]['sha256']==q['final_sha256']
for p in [B/'liblammps.so',B/'lmp',B/'CMakeCache.txt',PROD,Path('/workspace/y_clot/env.sh'),PRIOR/'BUGFIX_SOURCE_DIFF.patch']:
 files[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
assert sha(B/'liblammps.so')==EXPECTED==prior['diagnostic_builds']['build_owner_fix']['liblammps.so']
assert sha(PROD)=='d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a'
base=(S/'src/BPM/bond_bpm.cpp').read_text();assert 'if (id >= 0 && id < nlocal)' in base and 'if (jd >= 0 && jd < nlocal)' in base
orientation={'pwd':str(Path.cwd()),'python':sys.executable,'files':files,'commands':{}}
for name,args in {'branch':['git','branch','--show-current'],'head':['git','rev-parse','HEAD'],'git_status':['git','status','--short'],'native_head':['git','-C',str(S),'rev-parse','HEAD'],'cpu':['lscpu'],'memory':['free','-b'],'disk':['df','-B1','.'],'gpu':['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu','--format=csv'],'jobs':['ps','-eo','pid,comm,pcpu,pmem','--sort=-pcpu'],'mpi':['mpiexec','--version'],'ldd_final':['ldd',str(B/'liblammps.so')]}.items():
 q=subprocess.run(args,capture_output=True,text=True);orientation['commands'][name]={'exit':q.returncode,'stdout':q.stdout,'stderr':q.stderr}
orientation['cgroup']={str(p):p.read_text() for p in [Path('/sys/fs/cgroup/cpu/cpu.cfs_quota_us'),Path('/sys/fs/cgroup/cpu/cpu.cfs_period_us'),Path('/sys/fs/cgroup/memory/memory.limit_in_bytes')] if p.exists()}
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',orientation)
# Protect every accepted result/input, the preceding fix provenance (except rebuildable source/build caches), physics and production environment.
protected={}
for dirname in ['verification','results','pd_clot']:
 for p in (R/dirname).rglob('*'):
  if not p.is_file() or O in p.parents or F in p.parents:continue
  if any(x in p.parts for x in ['source','source_owner_fix','build','build_owner_fix','__pycache__']):continue
  protected[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
for p in [PROD,Path('/workspace/y_clot/env.sh'),*list((S/'src/BPM').glob('bond_bpm*')),*list((R/'tools/lammps-contact-onset-v1/src/BPM').glob('bond_bpm*'))]:protected[str(p)]={'sha256':sha(p),'bytes':p.stat().st_size}
dump(O/'PROTECTED_START.json',protected)
audit={}
for name in ['body.data','in.lammps','request.json','LOAD_PROVIDER.json']:
 shutil.copy2(OLD/name,P/name);a=sha(OLD/name);assert sha(P/name)==a;audit[name]={'source':str(OLD/name),'copy':str(P/name),'sha256':a,'byte_identical':True}
oldfreeze=json.loads((PRIOR/'V3_replay/REPLAY_DESIGN_AND_IDENTITY.json').read_text());assert {k:v['sha256'] for k,v in audit.items()}==oldfreeze['input_sha256']
for name,hash0 in oldfreeze['non_target_python_sha256'].items():assert sha(R/'pd_clot/lammps_fragmentation'/name)==hash0
physics={str(p.relative_to(R)):sha(p) for p in (R/'pd_clot/lammps_fragmentation').glob('*.py')}
dump(O/'V3_INPUT_IDENTITY_AUDIT.json',{'status':'PASS','files':audit,'previous_accepted_replay_manifest':str(PRIOR/'V3_replay/REPLAY_DESIGN_AND_IDENTITY.json'),'physics_python_sha256':physics})
dump(O/'PATCH_PROVENANCE.json',{'prior_run':str(PRIOR),'patch':str(PRIOR/'BUGFIX_SOURCE_DIFF.patch'),'patch_sha256':sha(PRIOR/'BUGFIX_SOURCE_DIFF.patch'),'final_source':str(S),'final_build':str(B),'final_source_sha256':patch,'final_lib_sha256':EXPECTED,'new_native_source_fix':False,'production_replaced':False,'owned_guard_present':True})
dump(O/'DESIGN_FREEZE.json',{'task':C['task'],'run':C['run'],'formal_run_limit':1,'steps':4660,'dt_s':5e-7,'mpi_ranks':2,'inputs':'byte-identical actual accepted C1 deck/body/request/provider','physical_changes':'Only existing final two-fix LAMMPS library','library_sha256':EXPECTED,'monitor':'Independent read-only END_OF_STEP fix; original wall and early-flow calculations unchanged. Bulk audit file sinks redirected to /dev/null. CSV/JSONL flush each step; exact FP64 key NPZ states.','m_gate':'per-ID max relative error <=2e-12 using strictly positive independently calculated m; absolute error also recorded; no denominator floor','max_output_budget_bytes':600000000,'memory_budget_bytes':8000000000,'failure_rule':'STOP on input/library mismatch, bad m, topology/ownership inconsistency, nonfinite state, or native crash. Preserve logs and last completed END_OF_STEP. At most one minimal debug if necessary; NO second full V3.','theta_scheme':'UNCHANGED explicit time layer; known limitation','deployment':'Diagnostic-only; production lib and env.sh protected','restart_and_other_unexercised':'NOT_EXERCISED unless direct evidence'})
print(json.dumps({'run':C['run'],'protected_files':len(protected),'protected_bytes':sum(q['bytes'] for q in protected.values()),'source_and_input_identity':'PASS'}))
