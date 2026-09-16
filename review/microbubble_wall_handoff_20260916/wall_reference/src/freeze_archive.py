"""Freeze all deliverables, validate the manifest, write an external receipt."""
from pathlib import Path
import ast,hashlib,json,re,subprocess,datetime
R=Path(__file__).resolve().parents[1]
W=Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['work'])
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
required=['RIGIDMULTIBLOBSWALL_PROVENANCE.json','PYSTOKES_PROVENANCE.json',
    'RIGIDMULTIBLOBSWALL_CAPABILITY_AUDIT.json','PYSTOKES_CAPABILITY_AUDIT.json',
    'PYSTOKES_WALL_API_AUDIT.json','CLASSICAL_SPHERE_WALL_THEORY_AUDIT.md',
    'WALL_REFERENCE_CONVENTION.json','WALL_REFERENCE_RAW_RESULTS.csv',
    'WALL_REFERENCE_MOBILITY_MATRICES.h5','WALL_REFERENCE_RESISTANCE_MATRICES.h5',
    'WALL_MOBILITY_MATRIX.csv','WALL_RESISTANCE_MATRIX.csv',
    'WALL_REFERENCE_CROSS_COMPARISON.csv','NEAR_WALL_COMPARISON.csv',
    'WALL_REFERENCE_SCALING_AUDIT.json','WALL_HYDRODYNAMICS_REFERENCE_AUDIT_REPORT.md',
    'FINAL_SUMMARY.json','FINAL_TERMINAL_SUMMARY.txt','README_REPRODUCE.md']
vis=json.loads((R/'VISUALIZATION_PROVENANCE.json').read_text())
required += [v['file'] for v in vis['artifacts']]
assert len(vis['artifacts'])==8
for name in required:assert (R/name).is_file() and (R/name).stat().st_size>0,name
assert json.loads((R/'validation/INDEPENDENT_FINALIZER.json').read_text())['status']=='PASS'
for p in (R/'src').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
for p in R.rglob('*.json'):json.loads(p.read_text())
for p in R.glob('*.md'):
    for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
        if not target.startswith(('http:','https:','#')):assert (p.parent/target).is_file(),(p,target)
raw=json.loads((R/'provenance/FIRST_OUTPUT_IDENTITY.json').read_text())['sha256']
assert all(sha(R/p)==h for p,h in raw.items())
codes=sorted((R/'src').glob('*.py'))
(R/'CODE_SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(R)}\n' for p in codes))
files=sorted(p for p in R.rglob('*') if p.is_file() and p!=R/'SHA256SUMS')
assert not any(p.is_symlink() for p in files),'Archive must not depend on external symlinks'
(R/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(R)}\n' for p in files))
result=subprocess.run(['sha256sum','--check','--quiet','SHA256SUMS'],cwd=R,text=True,capture_output=True)
assert result.returncode==0,result.stderr+result.stdout
receipt={'status':'PASS','manifest_sha256':sha(R/'SHA256SUMS'),'manifest':str(R/'SHA256SUMS'),
    'files_verified':len(files),'code_files':len(codes),'bytes_covered':sum(p.stat().st_size for p in files),
    'required_outputs':'PASS','JSON_parse':'PASS','Python_syntax':'PASS','report_local_links':'PASS',
    'first_raw_output_identity':'PASS','generated_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'receipt_placement':'Outside result directory to avoid circular manifest dependency',
    'scientific_status':'PASS_WITH_LIMITATIONS; see report'}
(W/'RESULT_SHA256_VERIFICATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
