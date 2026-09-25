from pathlib import Path
import os,json,hashlib,shutil
A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination']);M=D/'sync_metadata'/C['snapshot'];V=Path('/home/lzy/projects/ulm_3D_vascular')
SKIP={'.git','.venv','venv','__pycache__','.pytest_cache','.ruff_cache','.mypy_cache','.codex_tmp','node_modules','build','dist','.cache','external','third_party','Ultraliser','external_reference','references','tmp','vessel_model','archive','outputs'}
EXT={'.py','.sh','.md','.txt','.toml','.yaml','.yml','.json','.xml','.csv','.cpp','.h','.f90','.patch','.ini'}
rows=[];changed=[]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
for base,dirs,files in os.walk(V,followlinks=False):
 dirs[:]=[n for n in dirs if n not in SKIP and not n.startswith('.') and not n.endswith('.egg-info')]
 for name in files:
  p=Path(base)/name;rel=p.relative_to(V);target=D/'vascular_network'/rel
  if p.is_symlink() or p.suffix not in EXT or name.startswith('.') or p.stat().st_size>25*1024**2:continue
  if str(rel)=='README.md':continue
  before=p.stat();h=sha(p);state='UNCHANGED'
  if not target.exists() or sha(target)!=h:
   target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);assert sha(target)==h;state='ADDED_OR_UPDATED';changed.append(str(rel))
  assert p.stat().st_mtime_ns==before.st_mtime_ns,'CONCURRENT_CHANGE_DETECTED '+str(p)
  rows.append(dict(source=str(p),target=str(target.relative_to(D)),sha256=h,status=state))
# A removed active module is retained in the parent commit, not as stale current code.
removed=[]
for row in json.loads((M/'inherited_source_recheck.json').read_text()):
 if row['status']=='SOURCE_NO_LONGER_PRESENT_RETAINED_HISTORICAL_SNAPSHOT' and row['target'].startswith('vascular_network/vascular_processing/'):
  p=D/row['target'];assert not Path(row['source']).exists()
  if p.exists():p.unlink();removed.append(row['target'])
record=dict(files=rows,added_or_updated=changed,removed_in_source_and_export=removed,bulk_topbrain_data_excluded=True,scope='Latest first-party vascular code/config/tests/docs snapshot; independent TopBrain branch not covered by the FEM/Particle 165-test gate')
(M/'vascular_source_refresh.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='files'},indent=2))
