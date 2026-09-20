#!/usr/bin/env python3
"""One-shot freeze of the pre-existing scientific and development evidence."""
import gzip,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,git_state,sha256,write_json,now
R=ROOT/'reports/sv1_3j'
assert not (R/'reference_manifest.json').exists(), 'Frozen manifest already exists'
for kind in ('reports','outputs','logs','configs','benchmarks'):(ROOT/kind/'sv1_3j').mkdir(parents=True,exist_ok=True)
shutil.copyfile('/home/lzy/.codex/attachments/014d479e-7fb9-40ca-84c5-cee2b20c4898/pasted-text.txt',R/'USER_REQUEST.txt')
files={}; scopes=[]
for kind in ('reports','outputs','logs','configs','benchmarks'):
    for p in (ROOT/kind).iterdir():
        if p.is_dir() and p.name!='sv1_3j':scopes.append(str(p.relative_to(ROOT)))
scopes+=['inputs']
for name in scopes:files.update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items()})
for name in ('configs','src','scripts','tests'):
    for k,v in inventory(ROOT/name)['files'].items():
        if 'sv13j' in k or 'sv1_3j' in k or '__pycache__' in k:continue
        files[name+'/'+k]=v
(R/'history_baseline.json.gz').write_bytes(gzip.compress(json.dumps({'timestamp':now(),'files':files,'scopes':scopes}).encode(),mtime=0))
old=inventory(ROOT.parent/'FEM');old['git']=git_state(ROOT.parent/'FEM')
(R/'old_fem_baseline.json.gz').write_bytes(gzip.compress(json.dumps(old).encode(),mtime=0))
ref=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
accepted=json.loads((ROOT/'reports/sv1_3/cpu_accepted_solution.json').read_text())
paths={f['path'] for f in ref['files']}
paths.update(('reports/sv1_3g/mpi_resolution.json','reports/sv1_3g/mpi_hard_gate.json','reports/sv1_3g/cuda_build.json','reports/sv1_3g/stage_result.json','reports/sv1_3g/native_artifact_mirror.json','scripts/run_gpu_mpi.sh'))
paths.update((accepted['path'],'reports/sv1_3/cpu_accepted_solution.json','configs/sv1_3/production_performance.yaml','configs/sv1_3/policy.json','configs/sv1_3/petsc_options.json','reports/sv1_3/gpu_environment.json','reports/sv1_3/cuda_petsc_build.json'))
paths.update(('reports/sv1_3h/stage_result.json','reports/sv1_3h/cuda12_source_manifest.json','reports/sv1_3h/cuda12_runtime.json','reports/sv1_3h/petsc_cuda12_build.json','reports/sv1_3h/final_environment.json','reports/sv1_3h/petsc_failure_audit.json','reports/sv1_3h/delivery_manifest.json','scripts/use_cuda12_gpu_env.sh'))
for f in ref['files']:assert sha256(ROOT/f['path'])==f['sha256'], 'REFERENCE_CHANGED: '+f['path']
write_json(R/'reference_manifest.json',{'status':'PASS','timestamp':now(),'source':'Stage SV1.3 artifacts (read, not prompt constants)','files':[{'path':p,'sha256':sha256(ROOT/p)} for p in sorted(paths)],'SV1.2':ref['accepted_solution'],'SV1.3':accepted,'physics':ref['physics'],'dt':ref['dt'],'BC':ref['BC'],'build':ref['build'],'official_source_git':git_state(ROOT/'external/svMultiPhysics'),'CPU_production_read_only':True,'compatibility_candidates':['12.3.2','12.2.2','12.1.1'],'PETSc_archive_sha256':json.loads((ROOT/'reports/sv1_3h/petsc_cuda12_build.json').read_text())['source_archive_sha256'],'MPI_prefix':json.loads((ROOT/'reports/sv1_3g/mpi_resolution.json').read_text())['prefix'],'existing_CUDA_versions':['13.2.86',json.loads((ROOT/'reports/sv1_3h/cuda12_runtime.json').read_text())['selected_release']],'driver':next(x['stdout'] for x in json.loads((ROOT/'reports/sv1_3h/final_environment.json').read_text())['probes'] if '--query-gpu=name,driver_version,memory.total,compute_cap' in x['command'])})
write_json(ROOT/'configs/sv1_3j/policy.json',{'frozen_before_runs':now(),'proof_steps':20,'initial_state':'t=0','GPU_ranks':1,'GPU_count':1,'OMP_NUM_THREADS':1,'velocity_relative_L2':1e-5,'pressure_relative_L2':1e-5,'Qin_relative':1e-6,'Qout_relative':1e-5,'mass_error_difference':1e-6,'max_velocity_relative':1e-5,'MPI_timeout_s':10,'MPI_rank1_repetitions':5,'MPI_fallback_max':0,'candidates':['12.3.2','12.2.2','12.1.1'],'stop_on_first_make_pass':True,'host_compiler_major':12,'allow_unsupported_compiler':False,'source_patch_allowed':False,'PETSc_upgrade_allowed':False,'PETSc_smoke_repetitions':3,'benchmark_min_repeats':2,'benchmark_reported_time':'median of all 2 or 3 observed repetitions (never fastest)' ,'benchmark_third_repeat_difference':0.1,'GPU_speedup_threshold':1.25,'GPU_B_allowed':False,'production_changes_allowed':False})
print(json.dumps({'historical_entries':len(files),'old_FEM_entries':len(old['files']),'reference_files':len(paths)}))
