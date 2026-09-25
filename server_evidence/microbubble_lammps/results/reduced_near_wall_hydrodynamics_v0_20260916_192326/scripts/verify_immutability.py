"""Read-only comparison of previous stages, supplied PDF, and compiled inputs."""
import hashlib
import json
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before=json.loads((ROOT/'provenance/PREVIOUS_STAGE_BEFORE.json').read_text());old=Path(before['root'])
files=[{'path':str(p.relative_to(old)),'sha256':sha(p),'size':p.stat().st_size} for p in sorted(old.rglob('*')) if p.is_file()]
assert files==before['files'],'STOP_STATE_DRIFT_PREVIOUS_LOCAL_STAGE'
(ROOT/'provenance/PREVIOUS_STAGE_AFTER.json').write_text(json.dumps({'root':str(old),'files':files},indent=2)+'\n')
remote_before=json.loads((ROOT/'provenance/PREVIOUS_REMOTE_STAGE_BEFORE.json').read_text())
code='''from pathlib import Path
import hashlib,json
old=Path(__OLD__)
files=[{'path':str(p.relative_to(old)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'size':p.stat().st_size} for p in sorted(old.rglob('*')) if p.is_file()]
print(json.dumps({'root':str(old),'files':files}))
'''.replace('__OLD__',repr(remote_before['root']))
p=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=10','-o','StrictHostKeyChecking=yes','-o','UpdateHostKeys=no','-o','ControlMaster=no','-o','ControlPath=none','vast4090','python3 -B -'],input=code,text=True,capture_output=True,timeout=35)
assert p.returncode==0,p.stderr;remote_after=json.loads(p.stdout)
assert remote_before==remote_after,'STOP_STATE_DRIFT_PREVIOUS_REMOTE_STAGE'
(ROOT/'provenance/PREVIOUS_REMOTE_STAGE_AFTER.json').write_text(json.dumps(remote_after,indent=2)+'\n')
source=json.loads((ROOT/'reference/materials/user_supplied_paper.json').read_text())
assert sha(Path(source['original_path']))==source['sha256'],'USER_SUPPLIED_PDF_CHANGED'
build=json.loads((ROOT/'provenance/REMOTE_BUILD_SOURCE_CHECK.json').read_text())['sha256']
compiled={p:h for p,h in build.items() if p.startswith(('src/','tests/')) or p in ['CMakeLists.txt','raw/PHASE1_INPUTS.csv']}
assert all(sha(ROOT/p)==h for p,h in compiled.items()),'COMPILED_SOURCE_DRIFT'
before_baseline=json.loads((ROOT/'provenance/BASELINE_BEFORE.json').read_text())
after_baseline=json.loads((ROOT/'provenance/BASELINE_AFTER.json').read_text())
assert before_baseline['git']==after_baseline['git'] and before_baseline['upstream_inputs']==after_baseline['upstream_inputs']
for row in (ROOT/'reference/CF2003_FREE_SHEAR_REFERENCE_V0.sha256').read_text().splitlines():
    checksum,name=row.split();assert sha(ROOT/'reference'/name)==checksum,'FROZEN_REFERENCE_DRIFT'
record={'previous_local_stage_unchanged':True,'previous_local_file_count':len(files),
        'previous_remote_stage_unchanged':True,'previous_remote_file_count':len(remote_after['files']),
        'old_stage_status':json.loads((old/'CURRENT_STATE.json').read_text())['status'],
        'user_supplied_pdf_unchanged':True,'compiled_inputs_match_final_sources':True,
        'compiled_input_count':len(compiled),'baseline_inputs_unchanged':True,'git_state_unchanged':True,'frozen_reference_sha256_match':True}
assert record['old_stage_status']=='BLOCKED'
(ROOT/'provenance/IMMUTABILITY_CHECK.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
