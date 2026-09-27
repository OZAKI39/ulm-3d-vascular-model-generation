"""Portable review using only files under wss_audit; no project checkout needed."""
from pathlib import Path
import sys,importlib.util,json,hashlib
import numpy as np
import pyvista as pv
sys.dont_write_bytecode=True
A=Path(__file__).resolve().parents[1];sys.path.insert(0,str(A/'evidence/source'))
spec=importlib.util.spec_from_file_location('production',A/'evidence/source/compute_field_diagnostics.py');p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
assert hashlib.sha256(Path(p.__file__).read_bytes()).hexdigest()=='efa12e218bcaf99ba07db7194d3b8da79251da709c0a54a97ad6f6179be6e438'
d=np.load(A/'inputs/H0/flow_arrays_si.npz');x=d['points_m'];t=d['tetra'];wall=d['boundary_triangles'][d['facet_tags']==1]
owner=p.boundary_owners(t,wall);_,area,n=p.wall_geometry(x,t,wall,owner)
w=np.linalg.norm(p.tangential_traction(p.p1_gradients(x,t,d['velocity_m_s'])[owner],n,.00345312),axis=1)
v,weight=p.nodal_average(wall,w,area,len(x));m=pv.read(A/'inputs/display/wall_wss_si.vtp')
raw=float(abs(w-m['WSS_raw_Pa']).max());display=float(abs(v[m['GlobalNodeID_zero_based']]-m['WSS_display_Pa']).max())
assert raw==display==0
result=dict(raw_error_Pa=raw,display_error_Pa=display,area_mean_Pa=float(np.average(w,weights=area)),min_Pa=float(w.min()),max_Pa=float(w.max()),portable_bundle_verified=True)
(A/'data/portable_review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
