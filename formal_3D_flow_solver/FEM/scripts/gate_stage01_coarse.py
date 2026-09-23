#!/usr/bin/env python3
"""Record the hard-gate and visualization checkpoint before medium meshing."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import sha256, timestamp, write_json

base=ROOT/"outputs/stage01/coarse"
qc=json.loads((base/"qc/geometry_qc.json").read_text())
assert qc["hard_gate_status"]=="PASS"
for rank in (1,2):
    assert json.loads((base/f"qc/reload_r{rank}.json").read_text())["status"]=="PASS"
tests=ET.parse(ROOT/"reports/stage01/coarse_pytest.xml").getroot().find("testsuite")
assert int(tests.attrib["failures"])==int(tests.attrib["errors"])==0
visual=ROOT/"reports/stage01/coarse_review/visualization_manifest.json"
manifest=json.loads(visual.read_text())
assert len(manifest["images"])==5
for name,entry in manifest["images"].items():
    assert sha256(visual.parent/name)==entry["sha256"]
gate={"timestamp":timestamp(),"hard_gate_status":"PASS", "coarse_mesh_sha256":sha256(base/"mesh/volume_mesh.npz"),
    "coarse_quality_status":qc["quality_status"], "coarse_pytest":tests.attrib,
    "coarse_visualization_manifest_sha256":sha256(visual),
    "human_review":"PENDING; automated visual inspection does not replace user review",
    "medium_authorization":"User permits nondegenerate slivers as advisory; all coarse HARD gates passed and five review figures generated before medium"}
write_json(ROOT/"inputs/stage01/coarse_gate.json",gate)
write_json(ROOT/"reports/stage01/coarse_gate.json",gate)
print(json.dumps(gate,indent=2))
