import json
from pathlib import Path
import numpy as np
import pytest
from fem3d.mesh_qc import audit_mesh

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session",params=["coarse","medium"])
def stage01_mesh(request):
    profile=request.param
    base=ROOT/"outputs/stage01"/profile
    with np.load(base/"mesh/volume_mesh.npz") as archive:
        data=dict(archive)
    with np.load(ROOT/"inputs/stage01/tagged_surface_si.npz") as archive:
        source=dict(archive)
    contract=json.loads((ROOT/"reports/stage00/source_contract.json").read_text())
    config=json.loads((ROOT/"configs/meshing.json").read_text())
    # Recompute from actual artifacts instead of trusting a previously saved PASS.
    qc=audit_mesh(data,source,contract,config)
    return {"root":ROOT,"base":base,"profile":profile,"data":data,"source":source,"contract":contract,"config":config,"qc":qc}
