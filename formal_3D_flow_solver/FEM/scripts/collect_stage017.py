#!/usr/bin/env python3
"""Collect only the new stage's compute artifacts; preserve all earlier stages."""
import shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
source=ROOT/'outputs/stage01_7/remote_return/outputs/stage01_7';records={}
for name in ['.']:
    for p in source.rglob('*'):
        if not p.is_file(): continue
        target=ROOT/'outputs/stage01_7'/p.relative_to(source)
        if target.exists() and sha256(target)!=sha256(p): raise RuntimeError('Refusing to replace existing evidence')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
        assert sha256(p)==sha256(target)
        records[str(target.relative_to(ROOT))]=sha256(target)
write_json(ROOT/'reports/stage01_7/artifact_manifest.json',{'timestamp':timestamp(),'source':'Vast compute; WSL-owned source','files':records})
print('Collected',len(records),'new-stage artifacts')
