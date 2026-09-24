#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
from fem3d.audit import sha256, timestamp, write_json
from fem3d.mesh_qc import audit_mesh

p=argparse.ArgumentParser()
p.add_argument("--profile",required=True,choices=("coarse","medium"))
a=p.parse_args()
out=ROOT/"outputs/stage01"/a.profile
data=np.load(out/"mesh/volume_mesh.npz")
source=np.load(ROOT/"inputs/stage01/tagged_surface_si.npz")
contract=json.loads((ROOT/"inputs/stage01/source_contract.json").read_text())
config=json.loads((ROOT/"configs/meshing.json").read_text())
result=audit_mesh(data,source,contract,config)
result.update(timestamp=timestamp(), mesh_sha256=sha256(out/"mesh/volume_mesh.npz"))
write_json(out/"qc/geometry_qc.json",result)
print(json.dumps(result,indent=2))
