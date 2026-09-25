"""Read-only baseline/remote probes; write only this stage's provenance."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import platform
import socket
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = Path('/home/lzy/projects/github_sync/ulm_microbubble_20260916_173258')
REF = Path('/home/lzy/projects/compre_output/wall_hydrodynamics_reference_audit/20260916_115322')
WALL = Path('/home/lzy/projects/compre_output/wall_hydrodynamics_v0/20260916_130547')
STATE = REPO/'review/server_environment_audit_20260916/SERVER_ENVIRONMENT_STATE.json'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['before','after'],required=True);args=parser.parse_args()
    s=json.loads(STATE.read_text())
    env=dict(__import__('os').environ,GIT_OPTIONAL_LOCKS='0')
    def git(*a):return subprocess.check_output(['git',*a],cwd=REPO,text=True,env=env).strip()
    branch=git('branch','--show-current');commit=git('rev-parse','HEAD');status=git('status','--short','--untracked-files=all').splitlines()
    expected_status=['?? review/server_environment_audit_20260916/SERVER_ENVIRONMENT_REPORT.md','?? review/server_environment_audit_20260916/SERVER_ENVIRONMENT_STATE.json']
    assert branch=='codex/microbubble-wall-audits-handoff-20260916-173258' and commit=='a466b77be19559bc7d1d91b49f88ab78a976c529' and sorted(status)==sorted(expected_status),'STOP_STATE_DRIFT'
    selected=[x for x in s['source_identity']['files'] if ('stable_rigid_sphere_near_field' in x['path'] and x['path'].endswith('rigid_math.cpp')) or x['path'].endswith(('wall_model.cpp','wall_model.hpp'))]
    rows=[];remote=[]
    for x in selected:
        p=Path(x['local_archive']);assert sha(p)==x['sha256'],'STOP_STATE_DRIFT'
        classification='PRODUCTION' if 'stable_rigid' in str(p) else 'AUDIT_ONLY'
        rows.append({'hostname':socket.gethostname(),'absolute_path':str(p),'classification':classification,'sha256':sha(p),'size':p.stat().st_size})
        remote.append({'path':x['path'],'expected_sha256':x['sha256'],'classification':classification})
    for name in ['classical_wall_reference.py','compare_oneill_table.py','rmbw_reference_adapter.py']:
        p=REF/'src'/name;q=REPO/'review/microbubble_wall_handoff_20260916/wall_reference/src'/name
        assert sha(p)==sha(q),'STOP_STATE_DRIFT'
        rows.append({'hostname':socket.gethostname(),'absolute_path':str(p),'classification':'REFERENCE','sha256':sha(p),'size':p.stat().st_size})
    table=WALL/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5'
    expected='71e62ef6bc22edc66ae0c4d1dcc224378c2f6c5d3848578f3fbedd2cd9a43f98'
    assert sha(table)==expected,'STOP_STATE_DRIFT'
    rows.append({'hostname':socket.gethostname(),'absolute_path':str(table),'classification':'REFERENCE','sha256':sha(table),'size':table.stat().st_size})
    remote.append({'path':'/workspace/microbubble_lammps/results/wall_hydrodynamics_v0_20260916_130547/tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5','expected_sha256':expected,'classification':'REFERENCE'})
    for p in [WALL/'contracts/RMBW_WALL_TABLE_CONTRACT.json',WALL/'contracts/MICROBUBBLE_WALL_HYDRODYNAMICS_V0_CONTRACT.json',REF/'WALL_REFERENCE_CONVENTION.json',REF/'validation/THEORY_SUPPORT_AUDIT.json',STATE,STATE.with_name('SERVER_ENVIRONMENT_REPORT.md'),REPO/'review/microbubble_wall_handoff_20260916/CURRENT_STATE.json',REPO/'review/microbubble_wall_handoff_20260916/CURRENT_PROJECT_STATUS.md']:
        rows.append({'hostname':socket.gethostname(),'absolute_path':str(p),'classification':'IMMUTABLE_HANDOFF_EVIDENCE','sha256':sha(p),'size':p.stat().st_size})
    code='''import socket,platform,hashlib,json,subprocess
from pathlib import Path
items=ITEMS
rows=[]
for x in items:
 p=Path(x['path']);h=hashlib.sha256(p.read_bytes()).hexdigest()
 assert h==x['expected_sha256'],'STOP_STATE_DRIFT'
 rows.append({'hostname':socket.gethostname(),'absolute_path':str(p),'classification':x['classification'],'sha256':h,'size':p.stat().st_size})
def version(a):return subprocess.check_output(a,text=True).splitlines()[0]
print(json.dumps({'hostname':socket.gethostname(),'OS':Path('/etc/os-release').read_text(),'platform':platform.platform(),'compiler':version(['g++','--version']),'cmake':version(['cmake','--version']),'MPI':version(['mpirun','--version']),'hdf5':version(['pkg-config','--modversion','hdf5']),'upstream_inputs':rows}))
'''.replace('ITEMS',repr(remote))
    r=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=10','-o','StrictHostKeyChecking=yes','-o','UpdateHostKeys=no','-o','ControlMaster=no','-o','ControlPath=none','vast4090','python3 -B -'],input=code,text=True,capture_output=True,timeout=35)
    assert r.returncode==0,r.stderr
    data=json.loads(r.stdout);assert data['hostname']=='f7c62a262077','STOP_STATE_DRIFT'
    rows+=data.pop('upstream_inputs')
    def ver(a):return subprocess.check_output(a,text=True).splitlines()[0]
    record={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'git':{'branch':branch,'commit':commit,'status':status},'local':{'hostname':socket.gethostname(),'OS':Path('/etc/os-release').read_text(),'platform':platform.platform(),'compiler':ver(['g++','--version']),'cmake':ver(['cmake','--version']),'MPI':ver(['/usr/bin/mpirun','--version'])},'remote':data,'upstream_inputs':rows}
    target=ROOT/'provenance'/('BASELINE_'+args.phase.upper()+'.json')
    assert not target.exists(),'RECORD_ALREADY_EXISTS'
    target.write_text(json.dumps(record,indent=2)+'\n')
    if args.phase=='before':
        import csv
        with (ROOT/'provenance/UPSTREAM_INPUTS.tsv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=['hostname','absolute_path','classification','sha256','size'],delimiter='\t');writer.writeheader();writer.writerows(rows)
    else:
        before=json.loads((ROOT/'provenance/BASELINE_BEFORE.json').read_text())
        assert rows==before['upstream_inputs'] and record['git']==before['git'],'STOP_STATE_DRIFT'
    print(json.dumps({'phase':args.phase,'verified_inputs':len(rows),'git_state_match':True,'remote_hostname':data['hostname']}))


if __name__=='__main__':main()
