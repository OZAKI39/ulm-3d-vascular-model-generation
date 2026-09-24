import json
import numpy as np
import pyvista as pv
def test_old_new_mesh_identity(report):
    c=json.loads((report/'data/old_new_flow_contract.json').read_text())
    a,b=[pv.read(c['paths'][k]) for k in ['OLD','NEW']]
    for x,y in [(a.points,b.points),(a.cells,b.cells),(a.celltypes,b.celltypes)]:
        assert x.shape==y.shape and x.dtype==y.dtype and x.tobytes()==y.tobytes()
    for spec in c['boundaries'].values():
        a,b=map(pv.read,spec['paths'])
        for x,y in [(a.points,b.points),(a.faces,b.faces)]:assert x.tobytes()==y.tobytes()
        for attr in ['point_data','cell_data']:
            assert set(getattr(a,attr))==set(getattr(b,attr))
            for k in getattr(a,attr):assert np.array_equal(getattr(a,attr)[k],getattr(b,attr)[k])
