"""Freeze the completed Stage Q delivery; no CFD occurs here."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,now
R=ROOT/'reports/sv1_3q'
def read(n):return json.loads((R/(n+'.json')).read_text())
assert read('stage_status')['completed']
for n in ['preservation_audit','remote/remote_preservation','native_artifact_mirror','visuals']:assert read(n)['status']=='PASS'
assert read('pytest_summary')['stage_failures']==0
if read('winner')['full_required']:assert read('winner_steady_candidate')['status']=='PASS'
paths=[]
for folder in ['reports/sv1_3q','outputs/sv1_3q','logs/sv1_3q','configs/sv1_3q','patches/sv1_3q','benchmarks/sv1_3q','scripts/sv13q']:
 paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.name!='delivery_manifest.json' and '__pycache__' not in p.parts and 'plot_cache' not in p.parts)
paths.extend([ROOT/'src/sv_validation/sv13q.py',*sorted((ROOT/'tests').glob('*sv13q*.py'))])
files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha256(p)) for p in sorted(set(paths))]
d=dict(status='PASS',timestamp=now(),stage_status=read('stage_status')['status'],scientific_equivalence='DEFERRED',production='CPU_EARLY_STOP_PRODUCTION',production_unchanged=True,WSL_is_source_of_truth=True,files=files,source_manifest='source_patch.json',native_build_manifest='native_artifact_mirror.json',expensive_work_finished=True)
(R/'delivery_manifest.json').write_text(json.dumps(d,indent=2)+'\n');print('Stage Q delivery frozen:',len(files),'files')
