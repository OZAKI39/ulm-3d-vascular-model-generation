#!/usr/bin/env python3
"""Read exact input copies and return a hash receipt; no vascular compute."""
import json
import platform
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import sha256,timestamp,write_json

manifest=json.loads((ROOT/"inputs/stage00/input_manifest.json").read_text())
files=[]
for row in manifest["copies"]:
    digest=sha256(ROOT/row["project_path"])
    assert digest==row["sha256"],row["project_path"]
    files.append({"project_path":row["project_path"],"sha256":digest})
write_json(ROOT/"outputs/stage00/roundtrip_result.json",{
    "status":"PASS","timestamp":timestamp(),"hostname":platform.node(),
    "remote_work_dir":str(ROOT),"source_of_truth":manifest["source_of_truth"],"files":files,
    "vascular_geometry_meshed":False,"vascular_flow_solved":False})
print("Input byte checks passed on",platform.node(),"in",ROOT)
