"""Verify protected stacks, all stage sources and fixed inputs; mirror new artifacts."""
import shutil
from runner_remote import *
from environment_remote import snapshot
before=load('pre_install_environment');after=snapshot();write('final_environment',after)
keys=('driver_files','cuda_default_link','cuda13_prefix','cuda13_files','cuda126_prefix','cuda126_files','default_compilers','ld_configuration','historical_gpu_stacks')
checks={k:before[k]==after[k] for k in keys};write('remote_preservation',dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,includes_Stage_O=True));assert all(checks.values())
active=[]
for proc in Path('/proc').iterdir():
 if not proc.name.isdigit():continue
 try:exe=(proc/'exe').resolve(strict=True)
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if str(exe).startswith(str(BASE)) and exe.name in ('svmultiphysics','gpu_sparse_hypre','gpu_sparse_amgx'):active.append(int(proc.name))
write('native_solver_stopped',dict(status='PASS' if not active else 'FAIL',active_solver_pids=active));assert not active
patch=json.loads((BASE/'configs/source_patch.json').read_text());source_verified={};files=[]
for name in ['svmp_reuse_build','svmp_hypre_build','svmp_amgx_build']:
 if not (R/(name+'.json')).exists():continue
 d=load(name)
 if d['status']!='PASS':continue
 source=Path(d['source']);expected=patch['after'] if name=='svmp_reuse_build' else patch['before']
 for n,h in expected.items():assert digest(source/n)==h,n
 source_verified[name]=len(expected);files.append(dict(remote_path=d['executable'],local_path='outputs/sv1_3p/native/'+name+'/svmultiphysics',sha256=digest(d['executable'])))
 for p in [Path(d['build'])/'CMakeCache.txt',Path(d['build'])/'svMultiPhysics-build/CMakeCache.txt']:
  files.append(dict(remote_path=str(p),local_path='outputs/sv1_3p/native/'+name+'/'+p.parent.name+'_CMakeCache.txt',sha256=digest(p)))
for kind in ('hypre','amgx'):
 path=BASE/('external/petsc325_cuda13_'+kind);p=path/(kind+'.tar.gz')
 if p.exists():files.append(dict(remote_path=str(p),local_path='external/petsc325_cuda13_'+kind+'/'+p.name,sha256=digest(p)))
 p=path/'install/lib/libpetsc.so'
 if p.exists():files.append(dict(remote_path=str(p),local_path='external/petsc325_cuda13_'+kind+'/native_libpetsc.so',sha256=digest(p)))
 if kind=='hypre':
  for lib in sorted((path/'install/lib').glob('libHYPRE-*.so')):
   files.append(dict(remote_path=str(lib),local_path='external/petsc325_cuda13_hypre/'+lib.name,sha256=digest(lib)))
N=BASE.parent/'sv1_3n';manifest=json.loads((N/'configs/baseline_L_flow_input_manifest.json').read_text());verified=[]
for f in manifest['files']:
 if not f['destination'].startswith('SV_MESH/'):continue
 p=BASE/'outputs'/f['destination'];assert digest(p)==f['sha256'];verified.append(dict(path=str(p),sha256=f['sha256']))
policy=json.loads((BASE/'configs/policy.json').read_text());assert digest(N/'outputs/REAL_VASCULAR_GPU/1-procs/stFile_060.bin')==policy['checkpoint_sha256']
write('final_input_integrity',dict(status='PASS',files=verified,source_files_verified=source_verified,checkpoint_preserved=True))
for f in files:f['bytes']=Path(f['remote_path']).stat().st_size
write('native_artifacts',dict(status='PASS',files=files));print('Native Stage P preservation and artifact inventory PASS.',flush=True)
