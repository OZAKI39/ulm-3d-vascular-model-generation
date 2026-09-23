#!/usr/bin/env python3
"""Collect measured Stage 3 outputs and logs without replacing completed evidence."""
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
records={}
for category in ('outputs','logs'):
    source=ROOT/'outputs/stage03/remote_return'/category/'stage03'
    for p in sorted(source.rglob('*')):
        if not p.is_file():continue
        assert not p.is_symlink()
        target=ROOT/category/'stage03'/p.relative_to(source)
        if target.exists():assert sha256(target)==sha256(p),target
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,target)
        records[str(target.relative_to(ROOT))]=sha256(target)
write_json(ROOT/'reports/stage03/artifact_manifest.json',{'timestamp':timestamp(),'files':records})
print('Collected',len(records),'verified files.')
