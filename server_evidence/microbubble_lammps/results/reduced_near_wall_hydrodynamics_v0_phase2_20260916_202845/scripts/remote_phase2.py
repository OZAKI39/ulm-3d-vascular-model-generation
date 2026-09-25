"""Upload only the new stage, build with the frozen CPU library, and run bounded synthetic cases."""
import argparse,hashlib,io,json,subprocess,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
P=json.loads((ROOT/'provenance/STAGE_PATHS.json').read_text())
SSH=['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=10','-o','StrictHostKeyChecking=yes','-o','UpdateHostKeys=no','-o','ControlMaster=no','-o','ControlPath=none','vast4090']
def remote(code,timeout=45):
    p=subprocess.run(SSH+['python3 -B -'],input=code,text=True,capture_output=True,timeout=timeout)
    print(p.stdout,flush=True);assert p.returncode==0,p.stderr
    return p.stdout

def upload():
    files=[p for p in ROOT.rglob('*') if p.is_file() and (p.relative_to(ROOT).parts[0] in ['src','tests','scripts','contracts','cases','reference'] or p.name=='CMakeLists.txt')]
    data=io.BytesIO();manifest={}
    with tarfile.open(fileobj=data,mode='w:gz') as tar:
        for p in sorted(files):
            rel=str(p.relative_to(ROOT));manifest[rel]=hashlib.sha256(p.read_bytes()).hexdigest();tar.add(p,arcname=rel)
    r=subprocess.run(SSH+['tar -xzf - -C '+P['remote_stage']],input=data.getvalue(),capture_output=True,timeout=45);assert r.returncode==0,r.stderr
    remote('from pathlib import Path\nimport hashlib,json\nr=Path('+repr(P['remote_stage'])+')\nm='+repr(manifest)+'\nassert all(hashlib.sha256((r/p).read_bytes()).hexdigest()==h for p,h in m.items()),"UPLOAD_HASH_MISMATCH"\n(r/"provenance").mkdir(exist_ok=True)\n(r/"provenance/COMPILED_INPUTS.json").write_text(json.dumps(m,indent=2)+"\\n")\nprint("UPLOAD_VERIFIED",len(m))')

def pull():
    code='''from pathlib import Path
import sys,tarfile
r=Path(ROOT)
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for d in ['raw','validation','provenance']:
  for p in (r/d).rglob('*'):
   if p.is_file():tar.add(p,arcname=str(p.relative_to(r)))
'''.replace('ROOT',repr(P['remote_stage']))
    p=subprocess.run(SSH+['python3 -B -'],input=code.encode(),capture_output=True,timeout=90);assert p.returncode==0,p.stderr
    with tarfile.open(fileobj=io.BytesIO(p.stdout),mode='r:gz') as tar:
        for m in tar:
            target=ROOT/m.name;assert m.isfile() and target.resolve().is_relative_to(ROOT)
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(tar.extractfile(m).read())
    print('RAW_OUTPUTS_MIRRORED',flush=True)

parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','run','pull']);parser.add_argument('--names',default='ALL');args=parser.parse_args()
if args.mode=='build':
    upload()
    remote('''from pathlib import Path
import subprocess,json,time,os
P=PATHS
r=Path(P['remote_stage']);b=Path(P['remote_work'])/'build'
(r/'validation').mkdir(exist_ok=True);(r/'raw').mkdir(exist_ok=True)
steps=[('configure',['cmake','-S',str(r),'-B',str(b),'-DCMAKE_BUILD_TYPE=Release','-DENGINE_ROOT='+P['engine_root']]),('build',['cmake','--build',str(b),'-j','2']),('contract',[str(b/'test_phase2_contract'),P['table_remote'],str(r/'validation/CPP_CONTRACT.csv')])]
record=[]
for name,cmd in steps:
 t=time.monotonic();p=subprocess.run(cmd,capture_output=True,text=True,env=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1'))
 (r/f'provenance/{name}.log').write_text(p.stdout+p.stderr);record.append({'step':name,'command':cmd,'returncode':p.returncode,'seconds':time.monotonic()-t});print(name,p.returncode,(p.stdout+p.stderr)[-9000:],flush=True)
 if p.returncode:break
(r/'provenance/BUILD_RUN.json').write_text(json.dumps(record,indent=2)+'\\n')
assert len(record)==3 and all(x['returncode']==0 for x in record),'BUILD_OR_CONTRACT_FAILED'
'''.replace('PATHS',repr(P)),timeout=120)
    pull()
elif args.mode=='run':
    names=args.names.split(',')
    remote('''from pathlib import Path
import subprocess,json,time,os
P=PATHS;names=NAMES
r=Path(P['remote_stage']);b=Path(P['remote_work'])/'build';cases=json.loads((r/'contracts/CASE_INDEX.json').read_text())
record=[]
for c in cases:
 if names!=['ALL'] and c['name'] not in names:continue
 out=r/'raw'/c['name'];out.mkdir(exist_ok=False)
 cmd=['mpirun','--allow-run-as-root','--bind-to','none','-np',str(c['mpi_ranks']),str(b/'rigid_lmp_cpu'),str(r/'cases'/c['name']/'case.cfg')]
 t=time.monotonic();p=subprocess.run(cmd,cwd=out,text=True,capture_output=True,timeout=50,env=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
 (out/'engine.log').write_text(p.stdout+p.stderr)
 state=json.loads((out/'RUN_STATE.json').read_text()) if (out/'RUN_STATE.json').exists() else {}
 okay=(p.returncode==0 and state.get('status')=='COMPLETE') if c['expected']=='COMPLETE' else (p.returncode!=0 and state.get('reason')==c['expected'])
 row={'case':c['name'],'command':cmd,'returncode':p.returncode,'seconds':time.monotonic()-t,'expected':c['expected'],'outcome_match':okay,'state':state};record.append(row)
 (out/'RUN_COMMAND.json').write_text(json.dumps(row,indent=2)+'\\n');print(c['name'],p.returncode,'MATCH' if okay else 'FAIL',flush=True)
 if not okay:
  print((p.stdout+p.stderr)[-5000:],flush=True);break
(r/'provenance/RUN_BATCH_LAST.json').write_text(json.dumps(record,indent=2)+'\\n')
assert record and all(c['outcome_match'] for c in record),'ENGINE_RUN_FAILED'
'''.replace('PATHS',repr(P)).replace('NAMES',repr(names)),timeout=600)
    pull()
else:pull()
