"""Freeze this stage's delivery, excluding no historical writes and no new CFD."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256,now
R=ROOT/'reports/sv1_3n'
def read(n):return json.loads((R/(n+'.json')).read_text())
assert read('stage_result')['completed'] and read('stage_result')['status']=='GPU_DEVELOPMENT_PASS'
for n in ('preservation_audit','remote_preservation','WSL_source_mirror','native_artifact_mirror','gpu_steady_candidate','visuals'):
 assert read(n)['status']=='PASS',n
files=[]
for folder in ('reports/sv1_3n','outputs/sv1_3n','logs/sv1_3n','configs/sv1_3n','benchmarks/sv1_3n','patches/sv1_3n','scripts/sv13n'):
 for p in sorted((ROOT/folder).rglob('*')):
  if not p.is_file() or p.name=='delivery_manifest.json' or '__pycache__' in p.parts or 'plot_cache' in p.parts:continue
  files.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha256(p)))
for p in [ROOT/'src/sv_validation/sv13n.py',*sorted((ROOT/'tests').glob('*sv13n*.py'))]:
 files.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha256(p)))
data=dict(status='PASS',timestamp=now(),WSL_is_source_of_truth=True,classification='GPU_STEADY_CANDIDATE',scientific_validation='DEFERRED_BY_USER',supported_MPI_ranks=1,MPI_CUDA_2R='UNSUPPORTED: raw numeric failure retained',production='CPU_EARLY_STOP_PRODUCTION',production_unchanged=True,files=files,adopted_source_manifest='WSL_source_mirror.json',official_source_manifest='petsc325_source.json',native_build_archive_manifest='native_artifact_mirror.json',expensive_work_finished=True)
(R/'delivery_manifest.json').write_text(json.dumps(data,indent=2)+'\n')
print(f'Stage N delivery frozen: {len(files)} files; source/build manifests linked; no further CFD.')
