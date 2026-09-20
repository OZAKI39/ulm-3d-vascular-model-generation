#!/usr/bin/env python3
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import REPORT,load
from sv_validation.provenance import now,write_json
from remote import run_python
assert load('cpu_validation')['status']=='PASS'
assert (REPORT/'remote_agent_guide.txt').exists()
code=r'''
import subprocess,json,shutil,glob,os
from pathlib import Path
probes=[]
for name,args in [('nvidia-smi',['nvidia-smi']),('gpu-query',['nvidia-smi','--query-gpu=name,memory.total,driver_version,compute_cap,pcie.link.gen.current,pcie.link.gen.max,pcie.link.width.current,pcie.link.width.max','--format=csv,noheader']),('nvcc',['nvcc','--version']),('cpu',['lscpu','-J']),('mpi',['mpiexec','--version']),('cmake',['cmake','--version']),('gcc',['gcc','--version']),('memory',['free','-b'])]:
    exe=shutil.which(args[0])
    if not exe and args[0]=='nvcc' and Path('/usr/local/cuda/bin/nvcc').exists():exe='/usr/local/cuda/bin/nvcc'
    if not exe:probes.append({'name':name,'available':False,'command':args});continue
    args[0]=exe
    p=subprocess.run(args,capture_output=True,text=True,timeout=20)
    probes.append({'name':name,'available':True,'command':args,'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
cap=subprocess.run(['vast-capabilities'],capture_output=True,text=True,timeout=20)
capdata=json.loads(cap.stdout) if cap.returncode==0 else {}
print(json.dumps({'probes':probes,'hardware':capdata.get('hardware'),
    'workspace_is_volume':capdata.get('instance',{}).get('workspace_is_volume'),
    'cuda_paths':glob.glob('/usr/local/cuda*'),'profilers':{p:shutil.which(p) for p in ('nsys','ncu')},
    'source_root_exists':Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy').exists(),
    'stage_root_exists':Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3').exists()}))
'''
p=run_python(code,timeout=60)
write_json(REPORT/'gpu_remote_probe.json',{'timestamp':now(),'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
assert p.returncode==0,p.stderr
d=json.loads(p.stdout);query=next(r for r in d['probes'] if r['name']=='gpu-query')
assert query['exit_code']==0
fields=[s.strip() for s in query['stdout'].strip().split(',')]
result={'timestamp':now(),'status':'AVAILABLE' if '4090' in fields[0] else 'BLOCKED',
    'reason':None if '4090' in fields[0] else 'Target RTX4090 not available',
    'site':'existing remote native-compute instance','gpu_model':fields[0],
    'gpu_memory_bytes':int(fields[1].split()[0])*1024**2,'driver_version':fields[2],
    'compute_capability':fields[3],'PCIe':dict(zip(['current_gen','max_gen','current_width','max_width'],fields[4:])),
    'probes':d['probes'],'hardware':d['hardware'],'cuda_paths':d['cuda_paths'],
    'profilers':d['profilers'],'workspace_is_volume':d['workspace_is_volume'],
    'storage_policy':'WSL source of truth; all final evidence fetched to WSL; no instance lifecycle changes',
    'build_mode':'native; no Docker commands'}
write_json(REPORT/'gpu_environment.json',result)
print(json.dumps(result,indent=2))
