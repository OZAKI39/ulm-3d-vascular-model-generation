#!/usr/bin/env python3
"""Promote byte-checked returned artifacts into local review directories."""
import argparse
import json
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import sha256, timestamp, write_json
p=argparse.ArgumentParser()
p.add_argument("--profile",required=True,choices=("coarse","medium"))
a=p.parse_args()
source=ROOT/"outputs/stage01/remote_return/outputs/stage01"/a.profile
dest=ROOT/"outputs/stage01"/a.profile
records={}
for path in sorted(source.rglob("*")):
    if path.is_file():
        target=dest/path.relative_to(source)
        if target.exists() and sha256(target)!=sha256(path):
            raise RuntimeError(f"Refusing to replace different artifact {target}")
        target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():
            shutil.copy2(path,target)
        records[str(target.relative_to(ROOT))]=sha256(target)
for name in ("geometry_qc.json","reload_r1.json","reload_r2.json"):
    result=json.loads((dest/"qc"/name).read_text())
    assert result.get("hard_gate_status",result.get("status"))=="PASS"
assert records
write_json(ROOT/f"reports/stage01/{a.profile}_artifact_manifest.json",
           {"timestamp":timestamp(),"origin":"Vast compute; WSL-owned code and SI input", "files":records})
print(f"Collected {len(records)} verified artifacts for {a.profile}")
