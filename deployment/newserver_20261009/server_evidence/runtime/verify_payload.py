import argparse,hashlib,json,os,stat
from pathlib import Path
from datetime import datetime,timezone

ROOT=Path('/workspace/brava_cfd')
parser=argparse.ArgumentParser();parser.add_argument('--bootstrap',action='store_true');args=parser.parse_args()
rows=[json.loads(x) for x in (ROOT/'provenance/migration/SELECTED_MANIFEST.jsonl').open()]
if args.bootstrap:rows=[r for r in rows if r['destination'].startswith(('deps/','sources/','provenance/','project/gpu_solver_fix/'))]
changes_path=ROOT/'provenance/migration/EXPECTED_BUILD_CHANGES.json'
changes=json.loads(changes_path.read_text()) if changes_path.exists() and not args.bootstrap else {}
failures=[];count=0;size=0;regenerated=[]
for row in rows:
    p=ROOT/row['destination'];reason=None
    if row['kind']=='symlink':
        if not p.is_symlink() or os.readlink(p)!=row['target'] or not p.exists():reason='Missing/different/unresolved link'
    else:
        if not p.is_file() or p.is_symlink():reason='Missing/wrong file type'
        else:
            info=p.stat()
            with p.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
            if info.st_size!=row['size_bytes'] or digest!=row['sha256']:
                change=changes.get(row['destination'])
                if change and change['before_sha256']==row['sha256'] and digest==change['after_sha256'] and info.st_size==change['after_size_bytes']:
                    original=ROOT/change['original_preserved']
                    with original.open('rb') as f:original_sha=hashlib.file_digest(f,'sha256').hexdigest()
                    if original_sha!=row['sha256']:reason='Preserved original build configuration differs'
                    else:regenerated.append(dict(path=row['destination'],**change))
                else:reason='Size/SHA256 differs'
            if stat.S_IMODE(info.st_mode)!=row['mode']:reason='Mode differs'
            size+=info.st_size;count+=1
    if reason:failures.append(dict(path=row['destination'],reason=reason))
report=dict(status='FAIL' if failures else 'PASS',timestamp_utc=datetime.now(timezone.utc).isoformat(),scope='bootstrap' if args.bootstrap else 'all_selected_payload',files=count,bytes=size,symlinks=sum(r['kind']=='symlink' for r in rows),expected_build_regenerations=regenerated,byte_identical_files=count-len(regenerated),failures=failures)
name='BOOTSTRAP_INTEGRITY.json' if args.bootstrap else 'PAYLOAD_INTEGRITY.json'
(ROOT/'reports'/name).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));raise SystemExit(bool(failures))
