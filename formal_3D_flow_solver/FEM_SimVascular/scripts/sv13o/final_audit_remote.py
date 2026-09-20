"""Close native preservation and archive the adopted output-only solver build."""
import shutil
from runner_remote import *
from environment_remote import snapshot,inventory
before=load('pre_install_environment');after=snapshot();write('final_environment',after)
keys=('driver_files','cuda_default_link','cuda13_prefix','cuda13_files','cuda126_prefix','cuda126_files','default_compilers','ld_configuration','historical_gpu_stacks')
checks={k:before[k]==after[k] for k in keys};write('remote_preservation',dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,includes_Stage_N=True));assert all(checks.values())
patch=json.loads((BASE/'configs/source_patch.json').read_text());build=load('svmp_build');s=Path(build['source'])
active=[]
for proc in Path('/proc').iterdir():
 if not proc.name.isdigit():continue
 try:exe=(proc/'exe').resolve(strict=True)
 except (FileNotFoundError,PermissionError,ProcessLookupError):continue
 if str(exe)==build['executable']:active.append(int(proc.name))
write('native_solver_stopped',dict(status='PASS' if not active else 'FAIL',active_solver_pids=active,executable=build['executable']))
assert not active,'Stage O CFD must be finished before artifact closure'
for n,h in patch['after'].items():assert digest(s/n)==h,n
N=BASE.parent/'sv1_3n';manifest=json.loads((N/'configs/baseline_L_flow_input_manifest.json').read_text());verified=[]
for f in manifest['files']:
 if not f['destination'].startswith('SV_MESH/'):continue
 p=BASE/'outputs'/f['destination'];assert digest(p)==f['sha256'];verified.append(dict(path=str(p),sha256=f['sha256']))
write('final_input_integrity',dict(status='PASS',files=verified,source_files_verified=len(patch['after'])))
O=BASE/'outputs/native';O.mkdir();binary=O/'svmultiphysics';shutil.copyfile(build['executable'],binary)
for p in [Path(build['build'])/'CMakeCache.txt',Path(build['build'])/'svMultiPhysics-build/CMakeCache.txt']:
 shutil.copyfile(p,O/(p.parent.name+'_CMakeCache.txt'))
write('native_artifacts',dict(status='PASS',files=[dict(remote_path=str(p),sha256=digest(p),bytes=p.stat().st_size) for p in O.iterdir()],PETSc_reused_from_frozen_Stage_N=True))
print('Native Stage O preservation, source/input hashes and solver archive PASS.',flush=True)
