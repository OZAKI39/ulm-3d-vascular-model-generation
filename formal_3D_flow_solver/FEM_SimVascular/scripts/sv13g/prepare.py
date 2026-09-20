#!/usr/bin/env python3
"""One-shot freeze of the pre-existing scientific and development evidence."""
import gzip,json,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import inventory,git_state,sha256,write_json,now
R=ROOT/'reports/sv1_3g'
assert not (R/'reference_manifest.json').exists(), 'Frozen manifest already exists'
for kind in ('reports','outputs','logs','configs','benchmarks'):(ROOT/kind/'sv1_3g').mkdir(parents=True,exist_ok=True)
shutil.copyfile('/home/lzy/.codex/attachments/5a09016e-9810-4616-9e5e-32af3cb755ce/pasted-text.txt',R/'USER_REQUEST.txt')
files={}; scopes=[]
for kind in ('reports','outputs','logs','configs','benchmarks'):
    for p in (ROOT/kind).iterdir():
        if p.is_dir() and p.name!='sv1_3g':scopes.append(str(p.relative_to(ROOT)))
scopes+=['inputs']
for name in scopes:files.update({name+'/'+k:v for k,v in inventory(ROOT/name)['files'].items()})
for name in ('configs','src','scripts','tests'):
    for k,v in inventory(ROOT/name)['files'].items():
        if 'sv13g' in k or 'sv1_3g' in k or '__pycache__' in k:continue
        files[name+'/'+k]=v
(R/'history_baseline.json.gz').write_bytes(gzip.compress(json.dumps({'timestamp':now(),'files':files,'scopes':scopes}).encode(),mtime=0))
old=inventory(ROOT.parent/'FEM');old['git']=git_state(ROOT.parent/'FEM')
(R/'old_fem_baseline.json.gz').write_bytes(gzip.compress(json.dumps(old).encode(),mtime=0))
ref=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
accepted=json.loads((ROOT/'reports/sv1_3/cpu_accepted_solution.json').read_text())
paths={f['path'] for f in ref['files']}
paths.update((accepted['path'],'reports/sv1_3/cpu_accepted_solution.json','configs/sv1_3/production_performance.yaml','configs/sv1_3/policy.json','configs/sv1_3/petsc_options.json','reports/sv1_3/gpu_environment.json','reports/sv1_3/cuda_petsc_build.json'))
for f in ref['files']:assert sha256(ROOT/f['path'])==f['sha256'], 'REFERENCE_CHANGED: '+f['path']
write_json(R/'reference_manifest.json',{'status':'PASS','timestamp':now(),'source':'Stage SV1.3 artifacts (read, not prompt constants)','files':[{'path':p,'sha256':sha256(ROOT/p)} for p in sorted(paths)],'SV1.2':ref['accepted_solution'],'SV1.3':accepted,'physics':ref['physics'],'dt':ref['dt'],'BC':ref['BC'],'build':ref['build'],'official_source_git':git_state(ROOT/'external/svMultiPhysics'),'CPU_production_read_only':True})
write_json(ROOT/'configs/sv1_3g/policy.json',{'frozen_before_runs':now(),'proof_steps':20,'initial_state':'t=0','GPU_ranks':1,'GPU_count':1,'OMP_NUM_THREADS':1,'velocity_relative_L2':1e-5,'pressure_relative_L2':1e-5,'Qin_relative':1e-6,'Qout_relative':1e-5,'mass_error_difference':1e-6,'MPI_timeout_s':10,'MPI_rank1_repetitions':5,'MPI_fallback_max':1,'PETSc_smoke_repetitions':3,'benchmark_min_repeats':2,'benchmark_third_repeat_difference':0.1,'GPU_speedup_threshold':1.25,'GPU_B_allowed':False,'production_changes_allowed':False})
print(json.dumps({'historical_entries':len(files),'old_FEM_entries':len(old['files']),'reference_files':len(paths)}))
