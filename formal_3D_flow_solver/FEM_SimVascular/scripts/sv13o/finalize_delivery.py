import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,now
R=ROOT/'reports/sv1_3o'
def read(n):return json.loads((R/(n+'.json')).read_text())
assert read('stage_result')['completed']
for n in ('optimized_steady_candidate','preservation_audit','remote_preservation','native_artifact_mirror','visuals'):assert read(n)['status']=='PASS'
paths=[]
for folder in ('reports/sv1_3o','outputs/sv1_3o','logs/sv1_3o','configs/sv1_3o','patches/sv1_3o','benchmarks/sv1_3o','scripts/sv13o'):
 paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.name!='delivery_manifest.json' and '__pycache__' not in p.parts and 'plot_cache' not in p.parts)
paths.extend([ROOT/'src/sv_validation/sv13o.py',*sorted((ROOT/'tests').glob('*sv13o*.py'))])
files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha256(p)) for p in sorted(set(paths))]
d=dict(status='PASS',timestamp=now(),stage_status=read('stage_result')['status'],classification=read('optimized_steady_candidate')['classification'],scientific_equivalence='DEFERRED',production='CPU_EARLY_STOP_PRODUCTION',production_unchanged=True,WSL_is_source_of_truth=True,files=files,adopted_source_manifest='source_patch.json',native_build_manifest='native_artifact_mirror.json',expensive_work_finished=True)
(R/'delivery_manifest.json').write_text(json.dumps(d,indent=2)+'\n');print(f'Stage O delivery frozen: {len(files)} files. No further CFD.')
