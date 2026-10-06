from pathlib import Path
import os,shutil,subprocess,json,hashlib,datetime
SRC=Path('/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin');DST=Path('/home/lzy/projects/computation_examples/_archives/github_sync_brava_20261006');META=DST/'sync_metadata/brava_p2_20261006';META.mkdir(parents=True,exist_ok=True)
def sha(p):
 if p.is_symlink():return hashlib.sha256(os.fsencode(os.readlink(p))).hexdigest()
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def copy(source,target):
 target.parent.mkdir(parents=True,exist_ok=True)
 if target.exists() or target.is_symlink():
  if target.is_symlink()==source.is_symlink() and sha(source)==sha(target):return
  if target.is_symlink():target.unlink()
 if source.is_symlink():target.symlink_to(os.readlink(source))
 else:shutil.copy2(source,target)
excluded_csv='population_v1/analysis/stabilization_local_conservation_review_20261006/STABILIZATION_TERMS.csv'
tracked=[p.decode() for p in subprocess.check_output(['git','ls-files','-z'],cwd=SRC).split(b'\0') if p]
files=[]
for name in tracked:
 if name==excluded_csv:continue
 s=SRC/name;dest='brava_flow_roi_18mlmin/'+name;copy(s,DST/dest);files.append(dict(path=dest,source=str(s),sha256=sha(s),bytes=s.lstat().st_size if not s.is_symlink() else len(os.fsencode(os.readlink(s))),kind='symlink' if s.is_symlink() else 'file',role='current BraVa tracked source/evidence'))
shared=Path('/home/lzy/projects/computation_examples/_shared/particle_workflow/particle_3d');shared_records=[]
for group in ['src','tests']:
 for s in sorted((shared/group).rglob('*')):
  if not s.is_file() or any(t in {'__pycache__','.pytest_cache'} for t in s.parts) or s.suffix=='.pyc':continue
  dest='ulm_particle_formal_p9a5/particle_3d/'+str(s.relative_to(shared));copy(s,DST/dest);record=dict(path=dest,source=str(s),sha256=sha(s),bytes=s.stat().st_size,kind='file',role='production particle source/test');files.append(record);shared_records.append(record)
(META/'INCLUDED_FILES.json').write_text(json.dumps(dict(source_BraVa_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=SRC,text=True).strip(),base_GitHub_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=DST,text=True).strip(),branch='sync/brava-p2-algebra-20261006',files=files),indent=2)+'\n')
large=json.loads((SRC/'EXCLUDED_LARGE_ARTIFACTS.json').read_text())
(META/'EXCLUDED_FILES.json').write_text(json.dumps(dict(status='KEEP_LOCAL_AND_SERVER_NO_DELETION',new_exclusions=[dict(path='brava_flow_roi_18mlmin/'+excluded_csv,bytes=(SRC/excluded_csv).stat().st_size,sha256=sha(SRC/excluded_csv),reason='Redundant CSV; lossless NPZ and exporter retained; all rows validated',replacement='STABILIZATION_TERMS.npz')],inherited_exclusions_manifest='brava_flow_roi_18mlmin/GITHUB_SYNC_EXCLUDED.json',inherited_large_manifest='brava_flow_roi_18mlmin/EXCLUDED_LARGE_ARTIFACTS.json',inherited_large_count=len(large['files']),excluded_categories=['Environments/caches and compiled binaries/object libraries','Third-party complete build/source trees; source identities and patches retained','Full PETSc first-system binary dumps and repeated transient solver outputs; summaries, receipts, final state NPZ retained','Old large field/audit tables and obsolete visualization copies as listed in inherited manifest'],no_LFS=True,no_chunking=True),indent=2)+'\n')
copy(Path('/tmp/brava_github_sync_20261006/CSV_NPZ_EQUIVALENCE.json'),META/'CSV_NPZ_EQUIVALENCE.json')
for n in ['prepare_snapshot.py','verify_redundant_csv.py']:copy(Path('/tmp/brava_github_sync_20261006')/n,META/n)
print('Snapshot files',len(files),'shared',len(shared_records),'new destination changes follow separate audit')
