import json, os, shutil
from pathlib import Path
from collections import Counter

root=Path(__file__).resolve().parents[1];stage=root/'stage'
rows=[json.loads(x) for x in (root/'inventory/SELECTED_MANIFEST.jsonl').open()]
kept=[];changes=[];removed=[]
for row in rows:
    path=stage/row['destination']
    if 'pytest_tmp' in Path(row['destination']).parts:
        removed.append(row['destination']);path.unlink();continue
    if row['kind']=='symlink':
        target=row['target']
        for prefix in ['/workspace/brava_flow_roi_18mlmin_20260928/','/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/']:
            if target.startswith(prefix):
                new=os.path.relpath(stage/'project'/target[len(prefix):],path.parent)
                path.unlink();path.symlink_to(new)
                changes.append(dict(path=row['destination'],original_target=target,new_target=new,reason='Relocate link inside dedicated project root; target file contents unchanged'))
                row['original_target']=target;row['target']=new
                break
    kept.append(row)
for path in sorted(stage.rglob('pytest_tmp'),reverse=True):
    if path.is_dir():shutil.rmtree(path)
local=Path('/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/flow')
if local.is_symlink():
    path=stage/'project/flow';target=os.readlink(local)
    if not path.is_symlink():path.symlink_to(target)
    kept.append(dict(destination='project/flow',kind='symlink',target=target,source=str(local),source_path=str(local),origin='current_local_git_worktree'))
broken=[r['destination'] for r in kept if r['kind']=='symlink' and not (stage/r['destination']).exists()]
assert not broken,broken
summary=json.loads((root/'inventory/SELECTION_SUMMARY.json').read_text())
counts=Counter()
for r in kept:
    if r['kind']=='file':counts[r['destination'].split('/')[0]]+=r['size_bytes']
summary.update(files=sum(r['kind']=='file' for r in kept),symlinks=sum(r['kind']=='symlink' for r in kept),logical_bytes=sum(counts.values()),bytes_by_category=dict(counts),relocated_links=len(changes),broken_links=0,removed_test_temporary_entries=len(removed))
(root/'inventory/SELECTION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
(root/'inventory/RELOCATED_LINKS.json').write_text(json.dumps(changes,indent=2)+'\n')
(root/'inventory/SELECTED_MANIFEST.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in kept))
bootstrap=[];project=[]
for r in kept:
    d=r['destination']
    (bootstrap if d.startswith(('deps/','sources/','provenance/','project/gpu_solver_fix/')) else project).append(d)
for name,paths in [('BOOTSTRAP_FILES.txt',bootstrap),('PROJECT_FILES.txt',project)]:
    (root/'inventory'/name).write_text('\n'.join(paths)+'\n')
(stage/'provenance/migration').mkdir(parents=True,exist_ok=True)
for name in ['SELECTED_MANIFEST.jsonl','SELECTION_SUMMARY.json','RELOCATED_LINKS.json','SOURCE_IDENTITY.json','LOCAL_SERVER_DIFFERENCES.json','NEW_SERVER_INITIAL.json','NEW_SERVER_AGENT_GUIDE.md']:
    shutil.copy2(root/'inventory'/name,stage/'provenance/migration'/name)
print(json.dumps(summary,indent=2))
