#!/usr/bin/env python3
import json
import shutil
from pathlib import Path
import numpy as np
from fem3d.mesh_input import adapt_surface
from fem3d.audit import write_json,sha256

root=Path(__file__).resolve().parents[1]
points,triangles,tags,metadata,contract=adapt_surface(root)
out=root/"inputs/stage01"
out.mkdir(parents=True,exist_ok=True)
np.savez_compressed(out/"tagged_surface_si.npz",points_m=points,triangles=triangles,facet_tags=tags)
metadata["adapted_surface_sha256"]=sha256(out/"tagged_surface_si.npz")
write_json(out/"surface_adapter.json",metadata)
shutil.copyfile(root/"reports/stage00/source_contract.json",out/"source_contract.json")
write_json(root/"reports/stage01/source_adapter.json",metadata)
print(json.dumps(metadata,indent=2))
