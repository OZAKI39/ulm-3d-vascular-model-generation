from common import *
import datetime,subprocess
assert not (HERE/'TASK.json').exists()
# Recheck every protected input at the actual final boundary.
for n,h in read(HERE/'INPUT_SHA256.json').items():assert sha(BASE/n)==h,n
old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json')
assert sha(old['actual_binary'])==old['actual_binary_sha256']
assert sha(old['PETSc']['path'])==old['PETSc']['sha256']
assert all(sha(RETAINED/n)==h for n,h in read(IDENTITY/'EXACT_DEPLOYED_SOURCE_MANIFEST.json').items())
processes=subprocess.check_output(['ps','-eo','pid,ppid,args'],text=True)
active=[s for s in processes.splitlines() if str(HERE/'runtime') in s and ('svmultiphysics_audit' in s or 'pressure_corrector_audit' in s) and 'ps -eo' not in s]
assert not active,active
put('FINAL_INTEGRITY.json',dict(status='PASS',UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),protected_inputs_verified=len(read(HERE/'INPUT_SHA256.json')),production_source_files_verified=2431,production_binary_unchanged=True,PETSc_unchanged=True,own_active_solver_processes=active,own_TASK_exists=False,final_cgroup=managed().cgroup_memory(),GPU=subprocess.check_output(['nvidia-smi','--query-gpu=memory.used,memory.free','--format=csv'],text=True),other_GPU_processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv'],text=True)))
files={str(p.relative_to(HERE)):sha(p) for p in HERE.rglob('*') if p.is_file() and not any(x in {'runtime','monitoring','__pycache__'} for x in p.relative_to(HERE).parts) and p.name not in {'OUTPUT_SHA256.json','TASK.json'}}
put('OUTPUT_SHA256.json',files);print('Closeout PASS',len(files),'deliverables')
