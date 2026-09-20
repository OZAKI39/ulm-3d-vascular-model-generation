"""Verify protected stacks, all stage sources and fixed inputs; mirror new artifacts."""
import shutil
from runner_remote import *
from environment_remote import snapshot
before=load('pre_install_environment');after=snapshot();write('final_environment',after)
keys=('driver_files','cuda_default_link','cuda13_prefix','cuda13_files','cuda126_prefix','cuda126_files','default_compilers','ld_configuration','historical_gpu_stacks')
checks={k:before[k]==after[k] for k in keys};write('remote_preservation',dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,includes_Stage_O=True,includes_Stage_P=True));assert all(checks.values())
active=[]
for proc in Path('/proc').iterdir():
 if not proc.name.isdigit():continue
 try:exe=(proc/'exe').resolve(strict=True)
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if str(exe).startswith(str(BASE)) and exe.name in ('svmultiphysics','recovery_probe'):active.append(int(proc.name))
write('native_solver_stopped',dict(status='PASS' if not active else 'FAIL',active_solver_pids=active));assert not active
patch=json.loads((BASE/'configs/source_patch.json').read_text());source_verified={};files=[]
for name in ['svmp_reuse_build']:
 if not (R/(name+'.json')).exists():continue
 d=load(name)
 if d['status']!='PASS':continue
 source=Path(d['source']);expected=patch['after'] if name=='svmp_reuse_build' else patch['before']
 for n,h in expected.items():assert digest(source/n)==h,n
 source_verified[name]=len(expected);files.append(dict(remote_path=d['executable'],local_path='outputs/sv1_3q/native/'+name+'/svmultiphysics',sha256=digest(d['executable'])))
 for p in [Path(d['build'])/'CMakeCache.txt',Path(d['build'])/'svMultiPhysics-build/CMakeCache.txt']:
  files.append(dict(remote_path=str(p),local_path='outputs/sv1_3q/native/'+name+'/'+p.parent.name+'_CMakeCache.txt',sha256=digest(p)))
N=BASE.parent/'sv1_3n';manifest=json.loads((N/'configs/baseline_L_flow_input_manifest.json').read_text());verified=[]
for f in manifest['files']:
 if not f['destination'].startswith('SV_MESH/'):continue
 p=BASE/'outputs'/f['destination'];assert digest(p)==f['sha256'];verified.append(dict(path=str(p),sha256=f['sha256']))
policy=json.loads((BASE/'configs/policy.json').read_text());assert digest(N/'outputs/REAL_VASCULAR_GPU/1-procs/stFile_060.bin')==policy['checkpoint_sha256']
assert digest(BASE.parent/policy['early_checkpoint_remote'])==policy['early_checkpoint_sha256']
write('final_input_integrity',dict(status='PASS',files=verified,source_files_verified=source_verified,checkpoint_preserved=True))
for remote_name,local_name in [('benchmarks/recovery_probe','recovery_probe'),('recovery_probe_profile.txt','recovery_probe_profile.txt')]:
 p=BASE/remote_name
 assert p.exists()
 files.append(dict(remote_path=str(p),local_path='benchmarks/sv1_3q/'+local_name,sha256=digest(p)))
for f in files:f['bytes']=Path(f['remote_path']).stat().st_size
write('native_artifacts',dict(status='PASS',files=files));print('Native Stage Q preservation and artifact inventory PASS.',flush=True)
