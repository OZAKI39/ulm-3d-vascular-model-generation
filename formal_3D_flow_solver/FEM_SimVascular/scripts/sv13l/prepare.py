#!/usr/bin/env python3
"""One-shot freeze of the pre-existing scientific and development evidence."""
import gzip,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,git_state,sha256,write_json,now
R=ROOT/'reports/sv1_3l'
assert not (R/'reference_manifest.json').exists(), 'Frozen manifest already exists'
for kind in ('reports','outputs','logs','configs','benchmarks'):(ROOT/kind/'sv1_3l').mkdir(parents=True,exist_ok=True)
shutil.copyfile('/home/lzy/.codex/attachments/099d4753-1297-458a-ba44-29c70ec30192/pasted-text.txt',R/'USER_REQUEST.txt')
files={}; scopes=[]
for kind in ('reports','outputs','logs','configs','benchmarks'):
    for p in (ROOT/kind).iterdir():
        if p.is_dir() and p.name!='sv1_3l':scopes.append(str(p.relative_to(ROOT)))
scopes+=['inputs']
for name in scopes:files.update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items()})
for name in ('configs','src','scripts','tests'):
    for k,v in inventory(ROOT/name)['files'].items():
        if 'sv13l' in k or 'sv1_3l' in k or '__pycache__' in k:continue
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
paths.update(str(p.relative_to(ROOT)) for p in (ROOT/'reports/sv1_3j').glob('*.json'))
paths.update(('reports/sv1_3j/REPORT.md','reports/sv1_3g/mpi_fallback_build.json','scripts/use_cuda123_env.sh'))
for f in ref['files']:assert sha256(ROOT/f['path'])==f['sha256'], 'REFERENCE_CHANGED: '+f['path']
write_json(R/'reference_manifest.json',{'status':'PASS','timestamp':now(),'source':'Stage SV1.3 artifacts (read, not prompt constants)','files':[{'path':p,'sha256':sha256(ROOT/p)} for p in sorted(paths)],'SV1.2':ref['accepted_solution'],'SV1.3':accepted,'physics':ref['physics'],'dt':ref['dt'],'BC':ref['BC'],'build':ref['build'],'official_source_git':git_state(ROOT/'external/svMultiPhysics'),'CPU_production_read_only':True,'CUDA_winner':'12.3.2', 'MPI_source':json.loads((ROOT/'reports/sv1_3g/mpi_fallback_build.json').read_text())['source'],'PETSc_archive_sha256':json.loads((ROOT/'reports/sv1_3h/petsc_cuda12_build.json').read_text())['source_archive_sha256'],'MPI_prefix':json.loads((ROOT/'reports/sv1_3g/mpi_resolution.json').read_text())['prefix'],'existing_CUDA_versions':['13.2.86',json.loads((ROOT/'reports/sv1_3h/cuda12_runtime.json').read_text())['selected_release']],'driver':next(x['stdout'] for x in json.loads((ROOT/'reports/sv1_3h/final_environment.json').read_text())['probes'] if '--query-gpu=name,driver_version,memory.total,compute_cap' in x['command'])})
policy=json.loads((ROOT/'configs/sv1_3j/policy.json').read_text())
for k in ('candidates','stop_on_first_make_pass','MPI_fallback_max'):policy.pop(k,None)
policy.update(frozen_before_runs=now(),CUDA_version='12.3.2',MPI_version='4.1.6',MPI_Fortran_bindings='all',compiler_major=12,datatype_probe_repeats=3,CPU_GPU_same_executable=True)
write_json(ROOT/'configs/sv1_3l/policy.json',policy)
print(json.dumps({'historical_entries':len(files),'old_FEM_entries':len(old['files']),'reference_files':len(paths)}))
