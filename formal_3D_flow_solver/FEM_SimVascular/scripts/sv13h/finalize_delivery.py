"""Hash only this stage's deliverables, then verify the report's local links."""
import ast, json, re, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,write_json,now
R=ROOT/'reports/sv1_3h'
paths=set()
for folder in ('reports/sv1_3h','outputs/sv1_3h','logs/sv1_3h','configs/sv1_3h',
               'benchmarks/sv1_3h','scripts/sv13h'):
    paths.update(p for p in (ROOT/folder).rglob('*') if p.is_file()
                 and '__pycache__' not in p.parts and 'plot_cache' not in p.parts)
paths.update((ROOT/'tests').glob('test_sv13h_*.py'))
paths.update(ROOT/p for p in ('tests/sv13h_support.py','src/sv_validation/sv13h.py','scripts/use_cuda12_gpu_env.sh'))
paths.discard(R/'delivery_manifest.json')
for p in paths:
    if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
write_json(R/'delivery_manifest.json',{'timestamp':now(),'status':'EVIDENCE_COMPLETE_STAGE_FAILED',
    'stage_status':'FAIL','reason':'PETSC_CUDA12_BUILD_FAIL','entries':len(paths),
    'files':[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'size':p.stat().st_size} for p in sorted(paths)],
    'scope':'Only SV1.3H additions; historical artifacts are covered by preservation_audit.json',
    'exclusions':['delivery_manifest.json itself','plot_cache','__pycache__']})
text=(R/'REPORT.md').read_text()
assert len(re.findall(r'^## ',text,re.M))==13
for link in re.findall(r'\]\(([^)]+)\)',text):
    if not link.startswith('https://'):assert (R/link).exists(),link
for e in json.loads((R/'delivery_manifest.json').read_text())['files']:
    assert sha256(ROOT/e['path'])==e['sha256']
print(f'Verified {len(paths)} deliverable hashes, 13 report sections, all local report links.')
