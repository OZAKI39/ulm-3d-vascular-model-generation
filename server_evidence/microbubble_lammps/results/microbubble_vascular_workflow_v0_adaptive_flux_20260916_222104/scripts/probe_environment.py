"""Read-only donor verification and one-time creation of this new remote stage."""
import json,hashlib,subprocess
from pathlib import Path
S=Path(__file__).resolve().parents[1]
p=json.loads((S/'provenance/STAGE_PATHS.json').read_text())
donor=Path(p['donor_local'])
old=json.loads((donor/'provenance/REMOTE_PREFLIGHT.json').read_text())
manifest={str(Path(k).relative_to(donor)):v for k,v in json.loads((S/'provenance/DONOR_HASHES_BEFORE.json').read_text()).items()}
payload={'paths':p,'original_remote_inputs':old['inputs'],'previous_stage_manifest':manifest,'hostname':old['hostname']}
code='''import json,hashlib,platform,subprocess
from pathlib import Path
payload=PAYLOAD
p=payload['paths'];errors=[]
for path,v in payload['original_remote_inputs'].items():
 if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=v['sha256']:errors.append(path)
for rel,v in payload['previous_stage_manifest'].items():
 path=Path(p['donor_remote'])/rel
 if hashlib.sha256(path.read_bytes()).hexdigest()!=v['sha256']:errors.append(str(path))
assert not errors,('STOP_BASELINE_DRIFT',errors)
assert platform.node()==payload['hostname'],'STOP_BASELINE_DRIFT: hostname'
for key in ['remote_stage','remote_work']:Path(p[key]).mkdir(exist_ok=False)
for folder in ['src','scripts','configs','contracts','geometry','fields','inputs','cases','raw','validation','visualization','provenance','report']:(Path(p['remote_stage'])/folder).mkdir()
print(json.dumps({'status':'PASS','hostname':platform.node(),'original_remote_inputs_verified':len(payload['original_remote_inputs']),'previous_stage_files_verified':len(payload['previous_stage_manifest']),'compiler':subprocess.check_output(['c++','--version'],text=True).splitlines()[0],'mpi':subprocess.check_output(['mpirun','--version'],text=True).splitlines()[0],'guide':Path('/etc/vast-agents-guide.md').read_text()},indent=2))
'''.replace('PAYLOAD',repr(payload))
(S/'provenance/REMOTE_PROBE_SCRIPT.py').write_text(code)
r=subprocess.run(['ssh',p['ssh_alias'],'python3','-'],input=code,text=True,capture_output=True,check=True)
(S/'provenance/REMOTE_PROBE_SSH.log').write_text(r.stderr)
res=json.loads(r.stdout);guide=res.pop('guide');(S/'provenance/VAST_AGENT_GUIDE.md').write_text(guide)
(S/'provenance/REMOTE_BASELINE_CHECK.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2));print('Guide matches previously read version:',hashlib.sha256(guide.encode()).hexdigest()==hashlib.sha256((donor/'provenance/VAST_AGENT_GUIDE.md').read_bytes()).hexdigest())
