#!/usr/bin/env python3
"""Full Stage 0 reference mechanism, writing ONLY Stage 1 evidence."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import reference_snapshot, sha256, timestamp, write_json

baseline_path=ROOT/"reports/stage00/reference_baseline.json"
baseline=json.loads(baseline_path.read_text())
current=reference_snapshot(baseline["roots"])
write_json(ROOT/"reports/stage01/reference_final.json",current)
old,new=baseline["files"],current["files"]
modified=sorted(k for k in old.keys()&new.keys() if old[k]!=new[k])
deleted=sorted(old.keys()-new.keys())
added=sorted(new.keys()-old.keys())
per_root={}
for root in baseline["roots"]:
    def in_root(paths):
        return [p for p in paths if Path(p).is_relative_to(root)]
    per_root[root]={"modified":len(in_root(modified)),"deleted":len(in_root(deleted)),"added":len(in_root(added))}
result={"status":"FAIL" if modified or deleted or added else "PASS", "timestamp":timestamp(),
    "reference_files_modified":len(modified),"reference_files_deleted":len(deleted),"reference_files_added":len(added),
    "modified":modified,"deleted":deleted,"added":added,"per_root":per_root,"entries":len(new),
    "baseline_sha256":sha256(baseline_path),"final_sha256":sha256(ROOT/"reports/stage01/reference_final.json"),
    "method":"Full SHA256 + size, mode, mtime recomputation; directories and symlinks inventoried; same mechanism as Stage 0"}
write_json(ROOT/"reports/stage01/reference_integrity.json",result)
previous=json.loads((ROOT/"reports/stage01/stage00_artifact_baseline.json").read_text())
old_artifacts=previous
now={}
for folder in ("reports/stage00","inputs/stage00","outputs/stage00","logs/stage00"):
    for path in (ROOT/folder).rglob("*"):
        if path.is_file():
            now[str(path.relative_to(ROOT))]={"sha256":sha256(path),"size":path.stat().st_size}
changes=sorted(k for k in old_artifacts.keys()|now.keys() if old_artifacts.get(k)!=now.get(k))
preservation={"status":"FAIL" if changes else "PASS", "timestamp":timestamp(),"entries":len(now),"changes":changes,
    "method":"All Stage 0 report, input, output and log file contents and sizes compared with Stage 1 start snapshot"}
write_json(ROOT/"reports/stage01/stage00_preservation.json",preservation)
print(json.dumps({"reference_integrity":result,"stage00_preservation":preservation},indent=2))
sys.exit(int(result["status"]!="PASS" or preservation["status"]!="PASS"))
