"""One-time, audited relocation. No dataset/model/record content rewriting."""
from pathlib import Path
import os,json,hashlib,csv,shutil,time
from concurrent.futures import ThreadPoolExecutor
OLD=Path('/home/lzy/projects/ulm_3D_vascular');NEW=Path('/home/lzy/projects/vascular_printing')
AUDIT=NEW/'migration/20260926';AUDIT.mkdir(parents=True,exist_ok=True)
CACHE={'__pycache__','.pytest_cache','.ruff_cache','.mypy_cache'}
whole=['vascular_processing','tools','third_party','Ultraliser','vessel_model','data','config','test_data']
whole += [p.name for p in OLD.glob('s1-*.py')]+['s2_swc_stl_model_generate.py','main.py','cfd_lumen_config.yaml','TOPBRAIN_BRAVA_LABEL_TRANSFER.md','requirements-fixture.txt','requirements-vascularmd.txt']
whole += ['utils/'+p.name for p in (OLD/'utils').iterdir() if p.name not in CACHE|{'__init__.py','cfd_flow','cfd_preprocess','cfd_surface_prepare'}]
whole += ['configs/'+p.name for p in (OLD/'configs').iterdir() if not p.name.startswith('cfd_')]
whole += ['outputs/'+p.name for p in (OLD/'outputs').iterdir() if p.is_dir() and not p.name.startswith('cfd_')]
whole += ['docs/ULTRALISER_PIPELINE.md','docs/research_history']
whole += ['references/'+p.name for p in (OLD/'references').iterdir() if p.name=='image' or (p.suffix=='.md' and (any(w in p.name for w in ['BraVa','VascularMD','TopBrain','SWC']) and 'CFD' not in p.name))]
keywords=('all_port_attachment','boundary_review_export','brava_input','compact_manufacturing','manufacturing','mevo','nn_production','o3_smooth','roi_display','s1_','swc_export','swc_roi','swc_stl','topbrain','vascularmd','rodent_vasculature','sampling_pipeline','schmid_pkl','ultraliser_backend','pipeline','mesh_cleanup','voxelization','hierarchical_graph','nne2','cli_smoke')
whole += ['tests/'+p.name for p in (OLD/'tests').iterdir() if p.name in {'data','swc_roi_inspection_output','README_swc_roi_inspection.md','inspect_swc_rois.py','schmid_test_data.py'} or (p.is_file() and p.name.startswith('test_') and any(k in p.name for k in keywords) and not p.name.startswith('test_nne2_staged'))]
# Include the NNE2 pipeline regression as part of the same model-generation family.
whole += ['tests/test_nne2_staged_pipeline.py']
whole=sorted(set(whole));copies=['utils/__init__.py','tests/conftest.py','requirements.txt']
assert all((OLD/r).exists() for r in whole+copies)
assert all(not (NEW/r).exists() for r in whole+copies)
plan=dict(source=str(OLD),destination=str(NEW),move=whole,copy_shared_boilerplate=copies,
    preserve_in_old=['CFD solvers/outputs/configs','network_1d0d','CFD reports/viewers','Git repository/history'],
    compatibility='Each moved source entry becomes a relative symlink to the new physical location. Frozen absolute-path records remain byte-identical.',
    inventory_excludes=['.git metadata','regenerable Python/test caches'])
(AUDIT/'migration_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2))
def entries(path):
 if path.is_symlink():yield path;return
 if path.is_file():yield path;return
 for root,dirs,files in os.walk(path,followlinks=False):
  dirs[:]=[d for d in dirs if d not in CACHE and d!='.git']
  for d in list(dirs):
   p=Path(root)/d
   if p.is_symlink():yield p;dirs.remove(d)
  for name in files:
   if name.endswith(('.pyc','.pyo')):continue
   yield Path(root)/name
paths=list(dict.fromkeys(p for r in whole+copies for p in entries(OLD/r)))
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while block:=f.read(8*1024*1024):h.update(block)
 return h.hexdigest()
def record(p):
 rel=str(p.relative_to(OLD))
 if p.is_symlink():return dict(path=rel,kind='symlink',target=os.readlink(p),bytes=0)
 return dict(path=rel,kind='file',bytes=p.stat().st_size,sha256=digest(p))
print('INVENTORY_START',len(paths),flush=True);start=time.monotonic();rows=[]
with ThreadPoolExecutor(max_workers=4) as pool,(AUDIT/'inventory_before.jsonl').open('w') as f:
 for i,row in enumerate(pool.map(record,paths),1):
  rows.append(row);f.write(json.dumps(row,ensure_ascii=False)+'\n')
  if i%2000==0:print('HASHED',i,'/',len(paths),flush=True)
print('INVENTORY_DONE',sum(r['bytes'] for r in rows),round(time.monotonic()-start,1),'seconds',flush=True)
links=[]
for rel in whole:
 src=OLD/rel;dst=NEW/rel;dst.parent.mkdir(parents=True,exist_ok=True)
 src.rename(dst);src.symlink_to(os.path.relpath(dst,src.parent),target_is_directory=dst.is_dir());links.append(dict(old=str(src),new=str(dst),target=os.readlink(src)))
for rel in copies:
 dst=NEW/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(OLD/rel,dst)
(AUDIT/'compatibility_links.json').write_text(json.dumps(links,ensure_ascii=False,indent=2))
# Existing TopBrain symlink originally targets the old location. Give the new
# project a direct internal relative link, keeping the same dataset target.
alias=NEW/'data/TopBrain';old_target=os.readlink(alias);resolved=alias.resolve()
assert resolved.is_relative_to(NEW) and resolved.is_dir()
alias.unlink();alias.symlink_to(os.path.relpath(resolved,alias.parent),target_is_directory=True)
(AUDIT/'symlink_updates.json').write_text(json.dumps([dict(path='data/TopBrain',before=old_target,after=os.readlink(alias),same_dataset=True)],indent=2))
issues=[]
def verify(row):
 p=NEW/row['path']
 if row['kind']=='file':return None if p.is_file() and digest(p)==row['sha256'] else row['path']
 if row['path']=='data/TopBrain':return None if p.resolve()==resolved else row['path']
 return None if p.is_symlink() and os.readlink(p)==row['target'] else row['path']
print('VERIFY_START',flush=True)
with ThreadPoolExecutor(max_workers=4) as pool:
 for issue in pool.map(verify,rows):
  if issue:issues.append(issue)
result=dict(status='PASS' if not issues else 'FAIL',file_count=len(rows),total_bytes=sum(r['bytes'] for r in rows),moved_entries=len(whole),shared_files_copied=copies,hash_mismatches=issues,elapsed_seconds=time.monotonic()-start)
(AUDIT/'relocation_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False),flush=True)
if issues:raise RuntimeError('Migration hash check failed')
