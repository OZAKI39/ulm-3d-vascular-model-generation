import json,shutil
import numpy as np
import pytest
from particle_3d.particle65_checkpoint import read_checkpoint
from particle_3d.audit import sha256

def test_restart_contract_cannot_silently_change(report,tmp_path):
    folder=tmp_path/'copy';shutil.copytree(report/'checkpoints/v1',folder)
    side=json.loads((folder/'particle_sidecar.json').read_text());provenance=side['provenance']
    side['contract']['h_molecular_floor_m']=3e-9
    (folder/'particle_sidecar.json').write_text(json.dumps(side))
    manifest=json.loads((folder/'checkpoint_manifest.json').read_text())
    manifest['files']['particle_sidecar.json']=sha256(folder/'particle_sidecar.json')
    (folder/'checkpoint_manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='CHECKPOINT_V1_PROVENANCE_MISMATCH'):
        read_checkpoint(folder,provenance,lambda i,s,t:(np.zeros(3),np.zeros(3)),.00345312)

def test_restart_provenance_cannot_silently_change(report):
    with pytest.raises(ValueError,match='CHECKPOINT_V1_PROVENANCE_MISMATCH'):
        read_checkpoint(report/'checkpoints/v1',{'different_frozen_input':True},lambda i,s,t:(np.zeros(3),np.zeros(3)),.00345312)
