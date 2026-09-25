from pathlib import Path
import subprocess,json
S=Path(__file__).resolve().parents[1];paths=json.loads((S/'provenance/STAGE_PATHS.json').read_text())
remote='''from pathlib import Path
import platform,subprocess,hashlib,json,os
s=Path(STAGE);w=Path(WORK)
assert platform.node()=='f7c62a262077','STATE_DRIFT_DETECTED'
for p in ['/workspace/microbubble_lammps','/workspace/hemocell_restore']:assert Path(p).is_dir(),'STATE_DRIFT_DETECTED'
s.mkdir(parents=True,exist_ok=False);w.mkdir(parents=True,exist_ok=False)
for d in ['src','scripts','configs','contracts','geometry','fields','inputs','cases','raw','validation','visualization','provenance','report']:(s/d).mkdir()
r=Path('/workspace/microbubble_lammps/results');base=r/'stable_rigid_sphere_near_field_v0_20260916_101617';wall=r/'wall_hydrodynamics_v0_20260916_130547';files=list((base/'src').glob('*'))+[wall/'src/wall_distance.cpp',wall/'src/wall_distance.hpp',wall/'fields/FROZEN_FLOW_FIELD_V0.h5',wall/'provenance/closed_geometry_m.stl',base/'CMakeLists.txt',base/'provenance/PROJECT_BUILD_PROVENANCE.json']
inputs={str(p):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in files if p.is_file()}
def cmd(a):
 r=subprocess.run(a,capture_output=True,text=True);return {'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
state={'hostname':platform.node(),'platform':platform.platform(),'inputs':inputs,'compiler':cmd(['c++','--version']),'mpi':cmd(['mpirun','--version']),'timestamp_utc':cmd(['date','-u','+%Y-%m-%dT%H:%M:%SZ']),'stage':str(s),'work':str(w),'historical_build':json.loads((base/'provenance/PROJECT_BUILD_PROVENANCE.json').read_text())}
(s/'provenance/REMOTE_PREFLIGHT.json').write_text(json.dumps(state,indent=2)+'\\n');(s/'provenance/VAST_AGENT_GUIDE.md').write_text(Path('/etc/vast-agents-guide.md').read_text());print(json.dumps(state))
'''.replace('STAGE',repr(paths['remote_stage'])).replace('WORK',repr(paths['remote_work']))
# Avoid substituting identifiers inside error strings; values are Python literals, not shell interpolation.
(S/'provenance/REMOTE_PREFLIGHT_SCRIPT.py').write_text(remote)
r=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','vast4090','python3','-'],input=remote,capture_output=True,text=True);(S/'provenance/REMOTE_PREFLIGHT_SSH.log').write_text(r.stderr);assert r.returncode==0,r.stderr
state=json.loads(r.stdout);(S/'provenance/REMOTE_PREFLIGHT.json').write_text(json.dumps(state,indent=2)+'\n');print(json.dumps({k:state[k] for k in ['hostname','stage','compiler','mpi']},indent=2))
