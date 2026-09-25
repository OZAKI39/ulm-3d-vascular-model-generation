"""Compile/run only this stage on Vast using existing CPU/HDF5/MPI libraries."""
import hashlib
import io
import json
import subprocess
import tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
paths=json.loads((ROOT/'provenance/STAGE_PATHS.json').read_text())
assert json.loads((ROOT/'validation/CF2003_REFERENCE_VALIDATION.json').read_text())['REFERENCE_RESOLUTION_PASS']
ssh=['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=10','-o','StrictHostKeyChecking=yes',
     '-o','UpdateHostKeys=no','-o','ControlMaster=no','-o','ControlPath=none','vast4090']
files=[p for p in ROOT.rglob('*') if p.is_file() and (p.parts[-2] in ['src','tests','scripts','contracts']
       or p.name in ['CMakeLists.txt','PHASE1_INPUTS.csv','INPUT_COVERAGE.json'])]
payload=io.BytesIO();manifest={}
with tarfile.open(fileobj=payload,mode='w:gz') as tar:
    for p in sorted(files):
        rel=str(p.relative_to(ROOT));manifest[rel]=hashlib.sha256(p.read_bytes()).hexdigest();tar.add(p,arcname=rel)
r=subprocess.run(ssh+['tar -xzf - -C '+paths['remote_stage']],input=payload.getvalue(),capture_output=True,timeout=50)
assert r.returncode==0,r.stderr
code='''from pathlib import Path
import subprocess,json,hashlib,time
r=Path(__REMOTE__);b=Path(__WORK__)/'build';manifest=__MANIFEST__
assert all(hashlib.sha256((r/p).read_bytes()).hexdigest()==h for p,h in manifest.items())
(r/'provenance/REMOTE_BUILD_SOURCE_CHECK.json').write_text(json.dumps({'match':True,'sha256':manifest},indent=2)+'\\n')
steps=[('configure',['cmake','-S',str(r),'-B',str(b),'-DCMAKE_BUILD_TYPE=Release']),
       ('build',['cmake','--build',str(b),'-j','2']),
       ('standalone',['str_binary',__TABLE__,str(r/'raw/PHASE1_INPUTS.csv'),str(r/'raw/CPP_PHASE1_OUTPUT.csv'),str(r/'validation/CPP_PHASE1_VALIDATION.json')])]
steps[2][1][0]=str(b/'test_reduced_wall')
results=[]
for name,cmd in steps:
 t=time.monotonic();p=subprocess.run(cmd,text=True,capture_output=True);elapsed=time.monotonic()-t
 (r/f'provenance/{name}.log').write_text(p.stdout+p.stderr)
 results.append({'step':name,'returncode':p.returncode,'elapsed_seconds':elapsed,'command':cmd})
 print(name,p.returncode,(p.stdout+p.stderr)[-10000:],flush=True)
 if p.returncode:break
(r/'provenance/BUILD_RUN.json').write_text(json.dumps(results,indent=2)+'\\n')
'''.replace('__REMOTE__',repr(paths['remote_stage'])).replace('__WORK__',repr(paths['remote_work'])).replace('__MANIFEST__',repr(manifest)).replace('__TABLE__',repr('/workspace/microbubble_lammps/results/wall_hydrodynamics_v0_20260916_130547/tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5'))
r=subprocess.run(ssh+['python3 -B -'],input=code,text=True,capture_output=True,timeout=120)
print(r.stdout);assert r.returncode==0,r.stderr
pull='''from pathlib import Path
import sys,tarfile
r=Path(__REMOTE__)
files=['provenance/configure.log','provenance/build.log','provenance/standalone.log','provenance/BUILD_RUN.json','provenance/REMOTE_BUILD_SOURCE_CHECK.json','raw/CPP_PHASE1_OUTPUT.csv','validation/CPP_PHASE1_VALIDATION.json']
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for name in files:
  if (r/name).exists():tar.add(r/name,arcname=name)
'''.replace('__REMOTE__',repr(paths['remote_stage']))
p=subprocess.run(ssh+['python3 -B -'],input=pull.encode(),capture_output=True,timeout=60);assert p.returncode==0,p.stderr
with tarfile.open(fileobj=io.BytesIO(p.stdout),mode='r:gz') as tar:
    for member in tar:
        target=ROOT/member.name
        assert member.isfile() and target.resolve().is_relative_to(ROOT)
        target.write_bytes(tar.extractfile(member).read())
record=json.loads((ROOT/'provenance/BUILD_RUN.json').read_text())
assert len(record)==3 and all(x['returncode']==0 for x in record),'REMOTE_BUILD_OR_TEST_FAILED'
print('Remote CPU build/test complete; outputs mirrored locally')
