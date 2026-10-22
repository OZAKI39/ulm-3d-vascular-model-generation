"""NumPy-only reproduction of actual local WSS subsets in the reviewer bundle."""
from pathlib import Path
import importlib.util,json,hashlib
import numpy as np
V=Path(__file__).resolve().parents[1]
module=V/'evidence/production_source/wss.py'
if not module.exists():module=V.parent/'solver_support/src/flow_solver_support/wss.py'
spec=importlib.util.spec_from_file_location('audited_production_wss',module);wss=importlib.util.module_from_spec(spec);spec.loader.exec_module(wss)
rows=[]
for path in sorted((V/'evidence/local_wall_snapshots').glob('*.npz')):
    meta=json.loads(path.with_suffix('.json').read_text());assert hashlib.sha256(path.read_bytes()).hexdigest()==meta['subset_sha256']
    assert hashlib.sha256(module.read_bytes()).hexdigest()==meta['core_module_sha256']
    f=np.load(path);x=f['points_m'];t=f['tetra'];wall=f['wall_triangles'];u=f['velocity_m_s'];mu=float(f['mu_Pa_s'])
    own=wss.boundary_owners(t,wall);_,area,n=wss.wall_geometry(x,t,wall,own)
    traction=wss.tangential_traction(wss.p1_gradients(x,t[own],u),n,mu);magnitude=np.linalg.norm(traction,axis=1)
    de=float(np.max(abs(magnitude-f['expected_WSS_raw_Pa'])));dv=float(np.max(abs(traction-f['expected_traction_Pa'])))
    assert max(de,dv)<1e-9
    rows.append(dict(case=path.stem,facets=len(wall),maximum_WSS_difference_Pa=de,maximum_vector_component_difference_Pa=dv,tolerance_Pa=1e-9,meaning='Recovery regression on actual solved-field subset, not a new CFD or physical accuracy test'))
assert rows
out=V/'data/local_subset_reproduction.json';out.parent.mkdir(exist_ok=True);out.write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
