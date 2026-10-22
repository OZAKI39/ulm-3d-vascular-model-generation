from pathlib import Path
import os,hashlib,json,shutil
new=Path('/home/lzy/projects/vascular_printing');old=Path('/home/lzy/projects/ulm_3D_vascular');audit=new/'migration/20260926'
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(4194304):h.update(b)
 return h.hexdigest()
paths=set();cache={'__pycache__','.pytest_cache','.ruff_cache','.mypy_cache'}
# Only caches in source trees; never traverse shared runtime, raw datasets or Git.
for root in [new/'vascular_processing',new/'tools',new/'tests',new/'utils',new/'third_party',old/'__pycache__',old/'.pytest_cache',old/'.ruff_cache']:
 if not root.exists():continue
 if root.name in cache:paths.update(p for p in root.rglob('*') if p.is_file() and not p.is_symlink());continue
 for folder,dirs,files in os.walk(root,followlinks=False):
  dirs[:]=[d for d in dirs if d not in {'.git','.venv'} and not (Path(folder)/d).is_symlink()]
  iscache=bool(set(Path(folder).parts)&cache)
  for name in files:
   p=Path(folder)/name
   if not p.is_symlink() and (iscache or name.endswith(('.pyc','.pyo'))):paths.add(p)
# Superseded axis/layout prototypes, excluding the unrelated CFD presentation.
scratch=old/'.codex_tmp'
if scratch.exists():
 for p in scratch.rglob('*'):
  if p.is_file() and not p.is_symlink() and 'cfd_ppt_boundary_review' not in p.parts:
   if p.name=='union_printer_probe.json':
    target=audit/'legacy_diagnostics'/p.name;target.parent.mkdir(exist_ok=True);shutil.copy2(p,target);assert digest(p)==digest(target)
   paths.add(p)
rows=[]
for p in sorted(paths):
 rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p),reason='Superseded visual-layout scratch' if scratch in p.parents else 'Regenerable Python/tool cache'))
(audit/'cleanup_deleted_files.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
for p in sorted(paths):p.unlink()
# Remove emptied directories only, deepest first; never remove a populated tree.
roots=[new/'vascular_processing',new/'tools',new/'tests',new/'utils',new/'third_party',old/'__pycache__',old/'.pytest_cache',old/'.ruff_cache',scratch]
for root in roots:
 if not root.exists():continue
 for folder,dirs,files in os.walk(root,topdown=False,followlinks=False):
  p=Path(folder)
  if not p.is_symlink() and not any(p.iterdir()):p.rmdir()
result=dict(deleted_files=len(rows),deleted_bytes=sum(r['bytes'] for r in rows),production_datasets_deleted=False,accepted_models_deleted=False)
(audit/'cleanup_summary.json').write_text(json.dumps(result,indent=2));print(result)
