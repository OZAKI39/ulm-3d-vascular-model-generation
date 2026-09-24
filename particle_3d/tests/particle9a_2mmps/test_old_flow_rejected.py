import json,pytest
from particle_3d.particle9a_provenance import FEM,OLD_FLOW_SHA,require_current_flow

def test_old_flow_rejected(tmp_path):
    manifest=json.loads((FEM/'frozen_reference/flow/flow_field_manifest.json').read_text())
    folder=tmp_path/'frozen_reference/flow';folder.mkdir(parents=True)
    for key,value in [('sha256',OLD_FLOW_SHA),('path','frozen_reference/flow/steady_flow_stage_sv1_3q.vtu'),('case_role','LEGACY_0P352841_MMPS')]:
        bad=dict(manifest);bad[key]=value
        (folder/'flow_field_manifest.json').write_text(json.dumps(bad))
        with pytest.raises(ValueError,match='REJECTS_LEGACY'):require_current_flow(tmp_path)
