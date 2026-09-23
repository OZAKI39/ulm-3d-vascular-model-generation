#!/usr/bin/env python3
"""Bring compute artifacts into WSL Stage 2; never copy remote source code."""
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import sha256,timestamp,write_json
returned=ROOT/"outputs/stage02/remote_return"
records={}
for category in ("api","meshes","cases"):
    source=returned/"outputs/stage02"/category
    for path in sorted(source.rglob("*")):
        if path.is_file():
            target=ROOT/"outputs/stage02"/category/path.relative_to(source)
            target.parent.mkdir(exist_ok=True,parents=True)
            if target.exists() and sha256(target)!=sha256(path):
                archive=ROOT/"outputs/stage02/previous_local_versions"/sha256(target)/target.name
                archive.parent.mkdir(exist_ok=True,parents=True)
                shutil.copy2(target,archive)
            shutil.copy2(path,target)
            assert sha256(path)==sha256(target)
            records[str(target.relative_to(ROOT))]=sha256(target)
for path in (returned/"reports/stage02").glob("*.json"):
    target=ROOT/"reports/stage02"/path.name
    shutil.copy2(path,target)
    records[str(target.relative_to(ROOT))]=sha256(target)
write_json(ROOT/"reports/stage02/artifact_manifest.json",{"timestamp":timestamp(),"files":records})
print(f"Collected {len(records)} Stage 2 artifacts")
