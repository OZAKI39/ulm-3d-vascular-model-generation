#!/usr/bin/env python3
"""Recompute production tetra QC from the frozen Stage 1.7 artifact."""
import json
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
from fem3d.mature_quality import evaluate_gmsh,summarize
from fem3d.mature_boundary import verify_boundary

B=ROOT/'outputs/stage01_7/selected';O=ROOT/'outputs/stage03/mesh'
assert not (O/'qc.json').exists()
lock=json.loads((ROOT/'inputs/stage03/mesh_input_verification.json').read_text())
assert all(sha256(B/p)==h for p,h in lock['source_artifact_sha256'].items())
raw=dict(np.load(B/'mesh/volume_mesh.npz'));old=raw['min_sicn'].copy()
boundary=verify_boundary(raw,dict(np.load(B/'surface/tagged_surface_si.npz')))
assert boundary['status']=='PASS' and boundary['maximum_boundary_displacement_m']==0
evaluator=evaluate_gmsh(raw)
np.testing.assert_allclose(raw['min_sicn'],old,rtol=0,atol=1e-12)
measured=summarize(raw,boundary)
baseline=json.loads((B/'qc/volume_quality.json').read_text())
assert measured['proxy']==baseline['proxy']
for key in ('min_sicn','total_below_0_1','cap_adjacent_below_0_1','low_quality_nearest_boundary_counts'):
    assert measured['quality'][key]==baseline['quality'][key]
assert all(v==0 for v in measured['validity'].values())
assert measured['topology']['connected_fluid_components']==1
assert set(map(int,np.unique(raw['facet_tags'])))=={1,2,3,4,5}
measured.update(status='PASS',timestamp=timestamp(),evaluator=evaluator,
                source_mesh_sha256=sha256(B/'mesh/volume_mesh.npz'),mesh_modified=False)
write_json(O/'qc.json',measured)
np.save(O/'min_sicn.npy',raw['min_sicn'])
print(json.dumps({k:measured[k] for k in ('status','proxy','quality','validity')},indent=2))
