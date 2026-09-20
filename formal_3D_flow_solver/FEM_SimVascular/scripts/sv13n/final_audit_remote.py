"""Read-only verification of frozen stacks and archival of Stage N native builds."""
import subprocess,tarfile
from environment_remote import snapshot
from runner_remote import *
before=load('pre_install_environment');after=snapshot();write('final_environment',after)
keys=('driver_files','cuda_default_link','cuda13_prefix','cuda13_files','cuda126_prefix','cuda126_files','default_compilers','ld_configuration','historical_gpu_stacks')
checks={k:before[k]==after[k] for k in keys}
write('remote_preservation',dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,historical_stages=['J','L','M','G MPI']))
assert all(checks.values())
flow_manifest=json.loads((BASE/'configs/baseline_L_flow_input_manifest.json').read_text())
flow_files=[]
for f in flow_manifest['files']:
 p=BASE/'outputs'/f['destination'];actual=digest(p)
 assert actual==f['sha256'],('FROZEN_NATIVE_INPUT_CHANGED',f['destination'])
 flow_files.append(dict(path=str(p),sha256=actual))
write('final_native_input_integrity',dict(status='PASS',files=flow_files,scope='Actual native mesh/geometry/template inputs checked against frozen original manifest; generated run-control XML changes recorded separately.'))
integrity=[]
for p in R.glob('petsc_*_source_integrity.json'):
 d=json.loads(p.read_text());s=Path(d['source']);changes=[n for n,h in d['files'].items() if not (s/n).is_file() or digest(s/n)!=h]
 integrity.append(dict(source=str(s),files=len(d['files']),changes=changes));assert not changes
for p in R.glob('svmp_*_source_expected.json'):
 d=json.loads(p.read_text());before=load(p.stem.replace('_source_expected','_source_before'));s=Path(before['source'])
 changes=[n for n,h in d['files'].items() if not (s/n).is_file() or digest(s/n)!=h]
 integrity.append(dict(source=str(s),files=len(d['files']),expected_compatibility_files=d['modified'],unexpected=changes));assert not changes
write('final_source_integrity',dict(status='PASS',sources=integrity))
O=BASE/'outputs/native';O.mkdir(exist_ok=True);paths=[]
for candidate in ('cpu','gpu13','gpu123'):
 p=R/('petsc_'+candidate+'_build.json')
 if not p.exists():continue
 w=json.loads(p.read_text())
 if w['status']!='PASS':continue
 archive=O/('petsc325_'+candidate+'.tar.gz');assert not archive.exists()
 with tarfile.open(archive,'w:gz') as t:t.add(w['prefix'],arcname='petsc325_'+candidate)
 paths.append(archive)
 p=R/('svmp_'+candidate+'_build.json')
 if p.exists():
  s=json.loads(p.read_text())
  if s['status']=='PASS':
   dst=O/('svmultiphysics_'+candidate);import shutil;shutil.copyfile(s['executable'],dst);paths.append(dst)
for p in B.iterdir():
 if p.is_file() and p.suffix not in ('.c','.cpp','.py'):paths.append(p)
write('native_artifacts',dict(status='PASS',artifacts=[dict(remote_path=str(p),sha256=digest(p),bytes=p.stat().st_size) for p in paths],dependency_note='Linux build-host binaries and RPATH; MPI/CUDA/VTK identities retained, not installed into WSL production.'))
print('Native history and unmodified PETSc / approved svMP source integrity verified.',flush=True)
